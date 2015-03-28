#include "NUMA_vector.h"
#include "matrix_config.h"
#include "mem_worker_thread.h"

const size_t num_eles = 1024 * 1024 * 10;
size_t num_nodes = 1;
size_t nthreads = 8;

using namespace fm;

void test_mapping()
{
	NUMA_vector::ptr vec = NUMA_vector::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::vector<size_t> lens = vec->get_mapper().cal_local_lengths(
			vec->get_length());
	size_t tot_len = 0;
	for (size_t i = 0; i < lens.size(); i++)
		tot_len += lens[i];
	assert(tot_len == vec->get_length());

	for (size_t i = 0; i < vec->get_length(); i++) {
		std::pair<int, size_t> phy_loc = vec->get_mapper().map2physical(i);
		size_t loc = vec->get_mapper().map2logical(phy_loc.first, phy_loc.second);
		assert(loc == i);
	}
}

void test_copy()
{
	NUMA_vector::ptr vec = NUMA_vector::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::unique_ptr<long[]> raw_arr(new long[vec->get_length()]);
	for (size_t i = 0; i < vec->get_length(); i++)
		raw_arr[i] = random();
	vec->copy_from((char *) raw_arr.get(), vec->get_length() * sizeof(long));
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(raw_arr[i] == vec->get<long>(i));
}

void test_sort()
{
	NUMA_vector::ptr vec = NUMA_vector::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::unique_ptr<long[]> raw_arr(new long[vec->get_length()]);
	for (size_t i = 0; i < vec->get_length(); i++)
		raw_arr[i] = random();
	vec->copy_from((char *) raw_arr.get(), vec->get_length() * sizeof(long));
	vec->sort();
	std::sort(raw_arr.get(), raw_arr.get() + vec->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(raw_arr[i] == vec->get<long>(i));
}

int main(int argc, char *argv[])
{
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		nthreads = atoi(argv[2]);
	}

	matrix_conf.set_num_nodes(num_nodes);
	matrix_conf.set_num_threads(nthreads);
	detail::mem_thread_pool::init_global_mem_threads(num_nodes,
			nthreads / num_nodes);
	test_mapping();
	test_copy();
	test_sort();
}
