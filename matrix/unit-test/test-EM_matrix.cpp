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
	EM_col_dense_matrix::ptr m = EM_col_dense_matrix::create(1024 * 1024, 5,
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
	printf("passed fetch/set test\n");
}

class verify_prod_compute: public submatrix_compute
{
	size_t start_row;
	const mem_dense_matrix &m;
public:
	verify_prod_compute(size_t start_row,
			const mem_dense_matrix &_m): submatrix_compute(start_row, 0), m(_m) {
		this->start_row = start_row;
	}

	virtual void run(const mem_dense_matrix &subm) {
		for (size_t i = 0; i < subm.get_num_rows(); i++) {
			for (size_t j = 0; j < subm.get_num_cols(); j++) {
				double *v = (double *) subm.get(i, j);
				double *expected_v = (double *) m.get(i + start_row, j);
				if (*v != *expected_v)
					printf("m[%ld,%ld]: %lf, expected: %lf\n",
							i, j, *v, *expected_v);
				assert(*v == *expected_v);
			}
		}
	}
};

void test_inner_prod()
{
	EM_col_dense_matrix::ptr em = EM_col_dense_matrix::create(1024 * 1024, 5,
			sizeof(double));
	mem_col_dense_matrix::ptr big_im = mem_col_dense_matrix::create(
			em->get_num_rows(), em->get_num_cols(), sizeof(double));
	mem_col_dense_matrix::ptr small_im = mem_col_dense_matrix::create(
			em->get_num_cols(), em->get_num_cols(), sizeof(double));

	// Init the big external-memory matrix
	size_t num_rows = 1024 * 8;
	EM_dense_matrix_accessor::ptr accessor = em->create_accessor();
	for (size_t i = 0; i < em->get_num_rows(); i += num_rows) {
		mem_dense_matrix::ptr mem_m = mem_col_dense_matrix::create(num_rows,
				em->get_num_cols(), em->get_entry_size());
		init_matrix(*mem_m, i * em->get_num_cols());
		accessor->set_submatrix(i, 0, mem_m);
	}
	accessor->wait4all();

	// Init the big in-memory matrix
	init_matrix(*big_im, 0);
	init_matrix(*small_im, 0);

	EM_dense_matrix::ptr em_res = multiply<double, double, double>(*em,
			*small_im);
	mem_dense_matrix::ptr im_res = multiply<double, double, double>(*big_im,
			*small_im);

	accessor = em_res->create_accessor();
	for (size_t i = 0; i < em_res->get_num_rows(); i += num_rows) {
		accessor->fetch_submatrix(i, num_rows, 0, em_res->get_num_cols(),
				submatrix_compute::ptr(new verify_prod_compute(i, *im_res)));
	}
	accessor->wait4all();
	// TODO we have to destroy the accessors first before we can destroy
	// the matrices and vectors.
	accessor = NULL;
	printf("passed inner product test\n");
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
	test_inner_prod();

	destroy_flash_matrix();
}
