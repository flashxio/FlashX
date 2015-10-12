#include "NUMA_vector.h"
#include "matrix_config.h"
#include "vector.h"
#include "mem_worker_thread.h"

const size_t num_eles = 1024 * 1024 * 1024;

using namespace fm;

int num_nodes = 1;
int nthreads = 8;

class rand_set: public type_set_vec_operate<long>
{
public:
	void set(long *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random();
	}
};

void test_sort()
{
	size_t len = (random() % num_eles) + num_eles;
	vector::ptr mem_vec = vector::create(len, get_scalar_type<long>(), true,
			rand_set());
	detail::NUMA_vec_store::ptr numa_vec = detail::NUMA_vec_store::create(len,
			num_nodes, get_scalar_type<long>());
	numa_vec->copy_from(
			dynamic_cast<const detail::mem_vec_store &>(mem_vec->get_data()).get_raw_arr(),
			mem_vec->get_length() * sizeof(long));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	numa_vec->sort();
	gettimeofday(&end, NULL);
	printf("sort %ld elements in NUMA vector takes %.3f\n", numa_vec->get_length(),
			time_diff(start, end));

	gettimeofday(&start, NULL);
	mem_vec->sort();
	gettimeofday(&end, NULL);
	printf("sort %ld elements in mem vector takes %.3f\n", numa_vec->get_length(),
			time_diff(start, end));
}

int main(int argc, char *argv[])
{
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		nthreads = atoi(argv[2]);
	}

	matrix_conf.set_num_nodes(num_nodes);
	matrix_conf.set_num_DM_threads(nthreads);
	detail::mem_thread_pool::init_global_mem_threads(num_nodes,
			nthreads / num_nodes);
	test_sort();
}
