#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include "sparse_matrix.h"
#include "matrix/FG_sparse_matrix.h"

using namespace fm;

void int_handler(int sig_num)
{
#ifdef PROFILER
	printf("stop profiling\n");
	if (!fg::graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void test_SpMV(sparse_matrix::ptr mat)
{
	printf("test sparse matrix vector multiplication\n");
	struct timeval start, end;
	type_mem_vector<double>::ptr in_vec
		= type_mem_vector<double>::create(mat->get_num_cols());
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_vec->set(i, i);
	gettimeofday(&start, NULL);
	type_mem_vector<double>::ptr out = mat->multiply<double>(in_vec);
	gettimeofday(&end, NULL);
	double in_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_sum += in_vec->get(i);
	double out_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		out_sum += out->get(i);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in_sum, out_sum, time_diff(start, end));
}

class mat_init_operate: public type_set_operate<double>
{
	size_t num_rows;
	size_t num_cols;
public:
	mat_init_operate(size_t num_rows, size_t num_cols) {
		this->num_rows = num_rows;
		this->num_cols = num_cols;
	}

	virtual void set(double *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		double start_val = row_idx * num_cols + col_idx;
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_val++;
	}
};

void test_SpMM(sparse_matrix::ptr mat, size_t mat_width)
{
	printf("test sparse matrix dense matrix multiplication\n");
	struct timeval start, end;
	type_mem_dense_matrix<double>::ptr in
		= type_mem_dense_matrix<double>::create(mat->get_num_cols(),
				mat_width, matrix_layout_t::L_ROW);
	in->get_matrix()->set_data(mat_init_operate(in->get_num_rows(),
				in->get_num_cols()));
	gettimeofday(&start, NULL);
	dense_matrix::ptr out = mat->multiply<double>(in->get_matrix());
	gettimeofday(&end, NULL);
	printf("it takes %.3f seconds\n", time_diff(start, end));
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr,
				"test conf_file matrix_file index_file [width]\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string matrix_file = argv[2];
	std::string index_file = argv[3];
	size_t mat_width = 1;
	if (argc >= 5)
		mat_width = atoi(argv[4]);
	signal(SIGINT, int_handler);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	SpM_2d_index::ptr index = SpM_2d_index::load(index_file);
	SpM_2d_storage::ptr mat_store = SpM_2d_storage::load(matrix_file, index);

	sparse_matrix::ptr mat = sparse_matrix::create(index, mat_store);
	if (mat_width == 1)
		test_SpMV(mat);
	else
		test_SpMM(mat, mat_width);

	destroy_flash_matrix();
}
