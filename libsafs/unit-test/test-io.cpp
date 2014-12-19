#include "safs_file.h"
#include "io_interface.h"
#include "cache.h"

using namespace safs;

static const size_t FILE_SIZE = 16 * 1024 * 1024;
static const size_t IO_SIZE = 256 * 1024;

/*
 * This generate a random request that is inside a RAID block.
 */
std::pair<off_t, size_t> get_rand_req_in_block()
{
	size_t size = 0;
	// Get a non-zero size
	while (size == 0)
		size = ROUND(random() % params.get_RAID_block_size(), sizeof(long));
	// Get an offset that is no larger than FILE_SIZE - size
	off_t off_in_block = ROUND(
			random() % (params.get_RAID_block_size() - size), sizeof(long));
	assert(off_in_block <= params.get_RAID_block_size());
	size_t num_blocks = FILE_SIZE / params.get_RAID_block_size();
	off_t off = off_in_block
		+ (random() % num_blocks) * params.get_RAID_block_size();
	assert(off + size <= FILE_SIZE);
	return std::pair<off_t, size_t>(off, size);
}

std::pair<off_t, size_t> get_rand_req()
{
	size_t size = 0;
	// Get a non-zero size
	while (size == 0)
		size = ROUND(random() % IO_SIZE, sizeof(long));
	// Get an offset that is no larger than FILE_SIZE - size
	off_t off = ROUND(random() % (FILE_SIZE - size), sizeof(long));
	return std::pair<off_t, size_t>(off, size);
}

std::pair<off_t, size_t> get_rand_align_req()
{
	auto p = get_rand_req();
	p.first = ROUND(p.first, 512);
	p.second = ROUNDUP(p.second, 512);
	assert(p.first + p.second <= FILE_SIZE);
	return p;
}

//////////////////////////////// Test direct compute //////////////////////////

/*
 * This test user compute runs for two I/O requests.
 * The first request is issued outside the user compute and the second request
 * is issued by the user compute.
 */
class test_compute: public user_compute
{
	int file_id;
	off_t first_off;
	size_t first_size;
	page_byte_array *first_arr;

	off_t second_off;
	size_t second_size;

	int num;
public:
	test_compute(compute_allocator *alloc): user_compute(alloc) {
		file_id = -1;
		this->first_off = 0;
		this->first_size = 0;
		this->first_arr = NULL;
		second_off = 0;
		second_size = 0;
		num = 0;
	}

	void set_first(int file_id, off_t off, size_t size) {
		this->file_id = file_id;
		this->first_off = off;
		this->first_size = size;
	}

	void verify(const page_byte_array &arr, off_t off, size_t size);

	virtual int serialize(char *buf, int size) const {
		return 0;
	}

	virtual int get_serialized_size() const {
		return 0;
	}

	virtual void run(page_byte_array &arr);

	virtual bool has_completed() {
		return num == 2;
	}

	virtual int has_requests() {
		// We haven't issed the second request.
		return second_size == 0;
	}

	virtual request_range get_next_request() {
		std::pair<off_t, size_t> p = get_rand_req();
		data_loc_t loc(file_id, p.first);
		size_t size = p.second;
		second_off = p.first;
		second_size = p.second;
		return request_range(loc, size, READ, this);
	}
};

void test_compute::verify(const page_byte_array &arr, off_t off, size_t size)
{
	assert(arr.get_offset() == off);
	assert(arr.get_size() == size);
	long expected = off / sizeof(long);
	page_byte_array::seq_const_iterator<long> it = arr.get_seq_iterator<long>(
			0, arr.get_size());
	PAGE_FOREACH(long, v, it) {
		assert(v == expected);
		expected++;
	} PAGE_FOREACH_END
}

void test_compute::run(page_byte_array &arr)
{
	num++;

	if (first_arr == NULL) {
		first_arr = arr.clone();
		return;
	}

	verify(*first_arr, first_off, first_size);
	verify(arr, second_off, second_size);
	byte_array_allocator &alloc = first_arr->get_allocator();
	alloc.free(first_arr);
}

class test_compute_allocator: public compute_allocator
{
public:
	virtual user_compute *alloc() {
		return new test_compute(this);
	}

	virtual void free(user_compute *compute) {
		delete compute;
	}
};

void test_direct_comp(const std::string &data_file)
{
	file_io_factory::shared_ptr factory = create_io_factory(data_file,
			DIRECT_COMP_ACCESS);
	io_interface::ptr io = factory->create_io(thread::get_curr_thread());
	test_compute_allocator alloc;
	for (int i = 0; i < 1000; i++) {
		test_compute *compute = (test_compute *) alloc.alloc();
		std::pair<off_t, size_t> p = get_rand_req();
		compute->set_first(io->get_file_id(), p.first, p.second);
		data_loc_t loc(io->get_file_id(), p.first);
		io_request req(compute, loc, p.second, READ);
		io->access(&req, 1);
		while (io->num_pending_ios() > 32)
			io->wait4complete(1);
	}
	while (io->num_pending_ios() > 0)
		io->wait4complete(io->num_pending_ios());
	printf("direct compute I/O passed the test.\n");
}

//////////////////////////////// Test remote IO ///////////////////////////////

class test_callback: public callback
{
	void verify(const io_request &req);
public:
	virtual int invoke(io_request *reqs[], int num) {
		for (int i = 0; i < num; i++)
			verify(*reqs[i]);
		return 0;
	}
};

void test_callback::verify(const io_request &req)
{
	assert(req.get_offset() % sizeof(long) == 0);
	assert(req.get_size() % sizeof(long) == 0);
	long expected = req.get_offset() / sizeof(long);
	long *vs = (long *) req.get_buf();
	int num_longs = req.get_size() / sizeof(long);
	for (int i = 0; i < num_longs; i++) {
		assert(vs[i] == expected);
		expected++;
	}
	free(req.get_buf());
}

void test_remote_io(const std::string &data_file)
{
	file_io_factory::shared_ptr factory = create_io_factory(data_file,
			REMOTE_ACCESS);
	io_interface::ptr io = factory->create_io(thread::get_curr_thread());
	io->set_callback(callback::ptr(new test_callback()));
	for (int i = 0; i < 1000; i++) {
		std::pair<off_t, size_t> p = get_rand_align_req();
		data_loc_t loc(io->get_file_id(), p.first);
		char *buf = NULL;
		int ret = posix_memalign((void **) &buf, 512, p.second);
		assert(ret == 0);
		io_request req(buf, loc, p.second, READ);
		io->access(&req, 1);
		while (io->num_pending_ios() > 32)
			io->wait4complete(1);
	}
	while (io->num_pending_ios() > 0)
		io->wait4complete(io->num_pending_ios());
	printf("remote I/O passed the test.\n");
}

std::string prepare_file()
{
	std::string data_file_name = basename(tempnam(".", "test"));;
	safs_file f(get_sys_RAID_conf(), data_file_name);
	f.create_file(FILE_SIZE);

	file_io_factory::shared_ptr factory = create_io_factory(data_file_name,
			REMOTE_ACCESS);
	io_interface::ptr io = factory->create_io(thread::get_curr_thread());
	long *buf = NULL;
	int ret = posix_memalign((void **) &buf, 4096, FILE_SIZE);
	assert(ret == 0);
	size_t num_longs = FILE_SIZE / sizeof(buf[0]);
	for (size_t i = 0; i < num_longs; i++)
		buf[i] = i;
	data_loc_t loc(io->get_file_id(), 0);
	io_request req((char *) buf, loc, FILE_SIZE, WRITE);
	io->access(&req, 1);
	io->wait4complete(1);
	printf("finish preparing the file (%s)\n", data_file_name.c_str());
	return data_file_name;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_io_system(configs);

	std::string data_file = prepare_file();
	test_remote_io(data_file);
	test_direct_comp(data_file);

	safs_file f(get_sys_RAID_conf(), data_file);
	f.delete_file();
	destroy_io_system();
}
