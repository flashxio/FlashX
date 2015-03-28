#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"

using namespace fm;

int num_nodes = 1;
int num_threads = 8;

class set_row_operate: public type_set_operate<long>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(long *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
};

void test_init()
{
	NUMA_row_tall_dense_matrix::ptr mat = NUMA_row_tall_dense_matrix::create(
			10000000, 10, num_nodes, get_scalar_type<long>());
	mat->reset_data();
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			assert(mat->get<long>(i, j) == 0);
	}

	set_row_operate op(mat->get_num_cols());
	mat->set_data(op);
	size_t val = 0;
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++) {
			assert(mat->get<long>(i, j) == val);
			val++;
		}
	}
}

int main(int argc, char *argv[])
{
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		num_threads = atoi(argv[2]);
	}

	matrix_conf.set_num_nodes(num_nodes);
	matrix_conf.set_num_threads(num_threads);
	detail::mem_thread_pool::init_global_mem_threads(num_nodes,
			num_threads / num_nodes);
	test_init();
}
