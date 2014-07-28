#include "simple_KV_store.h"

size_t arr_len;
std::unique_ptr<int32_t[]> int_arr;

size_t num_runs;

class data_check
{
	size_t idx;
public:
	data_check() {
		idx = 0;
	}

	data_check(size_t idx) {
		this->idx = idx;
	}

	size_t get_idx() const {
		return idx;
	}

	size_t get_num_entries() const {
		return 2;
	}

	bool merge(const data_check &task) {
		return false;
	}

	void run(page_byte_array::seq_const_iterator<int32_t> &it) {
		assert(it.get_num_tot_entries() == 2);
		assert(int_arr[idx] == it.next());
		assert(int_arr[idx + 1] == it.next());
		num_runs++;
	}

	bool operator<(const data_check &check) const {
		return this->idx < check.idx;
	}
};

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "test conf_file file_name\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string file_name = argv[2];

	config_map configs(conf_file);
	init_io_system(configs);
	file_io_factory::shared_ptr factory = create_io_factory(file_name,
			GLOBAL_CACHE_ACCESS);
	size_t file_size = factory->get_file_size();
	thread *curr = thread::get_curr_thread();
	io_interface::ptr io = factory->create_io(curr);

	arr_len = file_size / sizeof(int32_t);
	int_arr = std::unique_ptr<int32_t[]>(new int32_t[arr_len]);
	io->access((char *) int_arr.get(), 0, arr_len * sizeof(int32_t), READ);

	simple_KV_store<int32_t, data_check>::ptr store
		= simple_KV_store<int32_t, data_check>::create(io);
	int num_checks = 1000;
	for (int test = 0; test < 5; test++) {
		for (int i = 0; i < num_checks; i++) {
			size_t idx = random() % arr_len;
			data_check check(idx);
			store->async_request(check);
		}
		store->flush_requests();
		printf("There are %d pending IOs\n", io->num_pending_ios());
		while (io->num_pending_ios() > 0) {
			io->wait4complete(1);
			printf("There are %d pending IOs\n", io->num_pending_ios());
		}
		printf("There are %ld runs\n", num_runs);
	}
}
