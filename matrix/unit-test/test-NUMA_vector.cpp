#include "NUMA_vector.h"
#include "matrix_config.h"
#include "mem_worker_thread.h"

const size_t num_eles = 1024 * 1024 * 10;
size_t num_nodes = 1;
size_t nthreads = 8;

using namespace fm;
using namespace detail;

class set_seq_vec: public type_set_vec_operate<long>
{
public:
	virtual void set(long *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_idx + i;
	}
};

class set_rand_vec: public type_set_vec_operate<long>
{
public:
	virtual void set(long *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random();
	}
};

void test_init()
{
	printf("test init NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	vec->set_data(set_seq_vec());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert((size_t) vec->get<long>(i) == i);
}

void test_mapping()
{
	printf("test mapping in NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
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
	printf("test copying NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::unique_ptr<long[]> raw_arr(new long[vec->get_length()]);
	for (size_t i = 0; i < vec->get_length(); i++)
		raw_arr[i] = random();
	vec->copy_from((char *) raw_arr.get(), vec->get_length() * sizeof(long));
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(raw_arr[i] == vec->get<long>(i));

	NUMA_vec_store::ptr vec1 = NUMA_vec_store::create(vec->get_length(), num_nodes,
			get_scalar_type<long>());
	bool ret = vec1->copy_from(*vec);
	assert(ret);
	assert(vec1->get_length() == vec->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(vec1->get<long>(i) == vec->get<long>(i));
}

void test_deep_copy()
{
	printf("test deep copying NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	vec->set_data(set_seq_vec());

	NUMA_vec_store::ptr vec1 = NUMA_vec_store::cast(vec->deep_copy());
	NUMA_vec_store::ptr vec2 = NUMA_vec_store::cast(vec->deep_copy());
	vec->set_data(set_rand_vec());
	assert(vec1->get_length() == vec2->get_length());
	for (size_t i = 0; i < vec1->get_length(); i++)
		assert(vec1->get<long>(i) == vec2->get<long>(i));
}

void test_sort()
{
	printf("test sorting NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
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
	test_init();
	test_mapping();
	test_copy();
	test_deep_copy();
	test_sort();
}
