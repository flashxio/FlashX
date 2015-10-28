#include "in_mem_io.h"

using namespace safs;

size_t range_size_log = 20;
size_t range_size = 1 << range_size_log;

NUMA_buffer::ptr create_buf(size_t length, const NUMA_mapper &mapper)
{
	NUMA_buffer::ptr buf = NUMA_buffer::create(length, mapper);
	off_t off = 0;
	while (length > 0) {
		NUMA_buffer::data_info data = buf->get_data(off, length);
		assert(data.first);
		assert(data.second <= length);
		long *lptr = (long *) data.first;
		size_t lsize = data.second / sizeof(long);
		for (size_t i = 0; i < lsize; i++)
			lptr[i] = i + off / sizeof(long);
		off += data.second;
		length -= data.second;
	}
	assert((size_t) off == buf->get_length());
	return buf;
}

void test_in_mem(size_t length, const NUMA_mapper &mapper)
{
	printf("test a NUMA buffer with %ld bytes\n", length);
	NUMA_buffer::ptr buf = create_buf(length, mapper);
	std::unique_ptr<long[]> lbuf(new long[length / sizeof(long)]);
	char *raw_ptr = (char *) lbuf.get();
	buf->copy_to((char *) lbuf.get(), length, 0);
	for (size_t i = 0; i < length / sizeof(long); i++)
		assert(lbuf[i] == i);

	// Test copy from
	for (size_t i = 0; i < length / sizeof(long); i++)
		lbuf[i] = random();
	buf->copy_from((const char *) lbuf.get(), length, 0);
	
	// Test copy to.
	std::unique_ptr<long[]> lbuf1(new long[length / sizeof(long)]);
	buf->copy_to((char *) lbuf1.get(), length, 0);
	for (size_t i = 0; i < length / sizeof(long); i++)
		assert(lbuf1[i] == lbuf[i]);

	// Test with random locations.
	for (size_t i = 0; i < 1000; i++) {
		off_t off = random() % buf->get_length();
		size_t size = random() % (buf->get_length() - off);
		auto data = buf->get_data(off, size);
		assert(data.first);
		// There are two cases when getting part of data in the array.
		// We get all requested data.
		assert(data.second == size
				// The requested data reaches the end of a range.
				|| (off + data.second) % mapper.get_range_size() == 0);
		for (size_t i = 0; i < data.second; i++)
			assert(raw_ptr[off + i] == data.first[i]);
	}
}

void test_in_mem()
{
	size_t num_nodes = numa_num_configured_nodes();
	NUMA_mapper mapper(num_nodes, range_size_log);
	for (size_t i = 0; i < num_nodes; i++)
		test_in_mem(i * range_size + range_size / 2, mapper);
	for (size_t i = 0; i < num_nodes; i++)
		test_in_mem((i + 1) * range_size, mapper);
	for (size_t i = 0; i < num_nodes; i++)
		test_in_mem(i * range_size + 1, mapper);

	NUMA_mapper mapper1(1, range_size_log);
	test_in_mem(range_size / 2, mapper);
	test_in_mem(range_size + range_size / 2, mapper);
}

void test_load_save(size_t length)
{
	printf("test load and save\n");
	size_t num_nodes = numa_num_configured_nodes();
	NUMA_mapper mapper(num_nodes, range_size_log);

	NUMA_buffer::ptr buf = create_buf(length, mapper);
	std::unique_ptr<char[]> raw_buf(new char[length]);
	buf->copy_to(raw_buf.get(), length, 0);

	char *tmp_file = tempnam("/tmp/", "test");
	buf->dump(tmp_file);

	NUMA_buffer::ptr buf1 = NUMA_buffer::load(tmp_file, mapper);
	std::unique_ptr<char[]> raw_buf1(new char[length]);
	buf1->copy_to(raw_buf1.get(), length, 0);

	for (size_t i = 0; i < length; i++)
		assert(raw_buf[i] == raw_buf1[i]);

	int ret = unlink(tmp_file);
	assert(ret == 0);
}

void test_load_save()
{
	size_t num_nodes = numa_num_configured_nodes();
	for (size_t i = 0; i < num_nodes + 5; i++)
		test_load_save(i * range_size + range_size / 2);
}

int main()
{
	test_in_mem();
	test_load_save();
}
