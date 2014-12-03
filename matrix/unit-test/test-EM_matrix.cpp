#include <malloc.h>

#include "EM_dense_matrix.h"
#include "mem_dense_matrix.h"
#include "sparse_matrix.h"

using namespace fm;

class read_double_compute: public submatrix_compute
{
	size_t start_row;
public:
	read_double_compute(size_t start_row): submatrix_compute(start_row, 0) {
		this->start_row = start_row;
	}

	virtual void run(const mem_dense_matrix &subm) {
		size_t expected_v = start_row * subm.get_num_cols();
		for (size_t i = 0; i < subm.get_num_rows(); i++) {
			for (size_t j = 0; j < subm.get_num_cols(); j++) {
				double *v = (double *) subm.get(i, j);
				if (*v != expected_v)
					fprintf(stderr, "m[%ld,%ld]: %lf, expected: %ld\n",
							i, j, *v, expected_v);
				assert(*v == expected_v);
				expected_v++;
			}
		}
	}
};

void init_matrix(mem_dense_matrix &m, size_t start_val)
{
	for (size_t i = 0; i < m.get_num_rows(); i++) {
		for (size_t j = 0; j < m.get_num_cols(); j++) {
			double *v = (double *) m.get(i, j);
			*v = start_val;
			start_val++;
		}
	}
}

void test_fetch_set()
{
	EM_col_dense_matrix::ptr m = EM_col_dense_matrix::create(1024 * 1024, 50,
			sizeof(double));
	EM_dense_matrix_accessor::ptr accessor = m->create_accessor();
	size_t num_rows = 1024 * 8;
	for (size_t i = 0; i < m->get_num_rows(); i += num_rows) {
		mem_dense_matrix::ptr mem_m = mem_col_dense_matrix::create(num_rows,
				m->get_num_cols(), m->get_entry_size());
		init_matrix(*mem_m, i * m->get_num_cols());
		accessor->set_submatrix(i, 0, mem_m);
	}
	accessor->wait4all();

	for (size_t i = 0; i < m->get_num_rows(); i += num_rows) {
		accessor->fetch_submatrix(i, num_rows, 0, m->get_num_cols(),
				submatrix_compute::ptr(new read_double_compute(i)));
	}
	accessor->wait4all();
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	test_fetch_set();

	destroy_flash_matrix();
}
