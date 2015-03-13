#include <malloc.h>

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "matrix_config.h"
#include "EM_dense_matrix.h"
#include "mem_dense_matrix.h"
#include "sparse_matrix.h"

using namespace fm;

class set_col_operate: public type_set_operate<double>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(double *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

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

EM_dense_matrix::ptr test_EM_inner_prod(size_t nrow, size_t ncol)
{
	std::string matrix_name
		= boost::str(boost::format("test-%1%-%2%-D") % nrow % ncol);
	bool exist = EM_dense_matrix::exist(matrix_name);

	EM_col_dense_matrix::ptr em = EM_col_dense_matrix::create(nrow, ncol,
			get_scalar_type<double>(), matrix_name);
	mem_col_dense_matrix::ptr small_im = mem_col_dense_matrix::create(
			ncol, ncol, get_scalar_type<double>());

	// Init the big external-memory matrix
	if (!exist)
		em->set_data(set_col_operate(em->get_num_cols()));
	small_im->serial_set_data(set_col_operate(small_im->get_num_cols()));

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	struct timeval start, end;
	gettimeofday(&start, NULL);
	EM_dense_matrix::ptr em_res = multiply<double, double, double>(*em,
			*small_im);
	gettimeofday(&end, NULL);
	printf("multiply on EM matrix: %.3fs\n", time_diff(start, end));
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	return em_res;
}

mem_dense_matrix::ptr test_IM_inner_prod(size_t nrow, size_t ncol)
{
	mem_col_dense_matrix::ptr big_im = mem_col_dense_matrix::create(
			nrow, ncol, get_scalar_type<double>());
	mem_col_dense_matrix::ptr small_im = mem_col_dense_matrix::create(
			ncol, ncol, get_scalar_type<double>());
	// Init the big in-memory matrix
	big_im->set_data(set_col_operate(big_im->get_num_cols()));
	small_im->serial_set_data(set_col_operate(small_im->get_num_cols()));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	mem_dense_matrix::ptr im_res = par_multiply<double, double, double>(*big_im,
			*small_im);
	gettimeofday(&end, NULL);
	printf("multiply on IM matrix: %.3fs\n", time_diff(start, end));
	return im_res;
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

	test_EM_inner_prod(1024 * 1024 * 120, 5);
	test_IM_inner_prod(1024 * 1024 * 120, 5);

	destroy_flash_matrix();
}
