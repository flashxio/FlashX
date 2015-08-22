#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include "matrix_config.h"
#include "io_interface.h"
#include "safs_file.h"
#include "sparse_matrix.h"
#include "NUMA_dense_matrix.h"
#include "matrix/FG_sparse_matrix.h"

using namespace fm;

void int_handler(int sig_num)
{
#ifdef PROFILER
	printf("stop profiling\n");
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void test_SpMV(sparse_matrix::ptr mat)
{
	printf("test sparse matrix vector multiplication\n");
	struct timeval start, end;
	detail::NUMA_vec_store::ptr in_vec = detail::NUMA_vec_store::create(
			mat->get_num_cols(), matrix_conf.get_num_nodes(),
			get_scalar_type<double>());
#pragma omp parallel for
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_vec->set<double>(i, i);
	printf("initialize the input vector\n");

	// Initialize the output vector and allocate pages for it.
	gettimeofday(&start, NULL);
	detail::NUMA_vec_store::ptr out = detail::NUMA_vec_store::create(
			mat->get_num_rows(), matrix_conf.get_num_nodes(),
			get_scalar_type<double>());
	out->reset_data();
	gettimeofday(&end, NULL);
	printf("initialize a vector of %ld entries takes %.3f seconds\n",
			out->get_length(), time_diff(start, end));

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	printf("start SpMV\n");
	gettimeofday(&start, NULL);
	mat->multiply<double>(*in_vec, *out);
	gettimeofday(&end, NULL);
	printf("SpMV completes\n");
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	double in_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_sum += in_vec->get<double>(i);
	double out_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		out_sum += out->get<double>(i);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in_sum, out_sum, time_diff(start, end));
}

class mat_init_operate: public type_set_operate<double>
{
public:
	mat_init_operate(size_t num_rows, size_t num_cols) {
	}

	virtual void set(double *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = row_idx * (i + col_idx);
	}
};

void test_SpMM(sparse_matrix::ptr mat, size_t mat_width)
{
	printf("test sparse matrix dense matrix multiplication\n");
	struct timeval start, end;
	detail::NUMA_row_tall_matrix_store::ptr in
		= detail::NUMA_row_tall_matrix_store::create(mat->get_num_cols(),
				mat_width, matrix_conf.get_num_nodes(),
				get_scalar_type<double>());
	in->set_data(mat_init_operate(in->get_num_rows(), in->get_num_cols()));
	printf("set input data\n");

	// Initialize the output matrix and allocate pages for it.
	detail::NUMA_row_tall_matrix_store::ptr out
		= detail::NUMA_row_tall_matrix_store::create(mat->get_num_rows(),
				mat_width, matrix_conf.get_num_nodes(),
				get_scalar_type<double>());
	out->reset_data();
	printf("reset output data\n");

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	printf("Start SpMM\n");
	gettimeofday(&start, NULL);
	mat->multiply<double>(*in, *out);
	gettimeofday(&end, NULL);
	printf("SpMM completes\n");
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("it takes %.3f seconds\n", time_diff(start, end));

	for (size_t k = 0; k < in->get_num_cols(); k++) {
		double in_sum = 0;
		for (size_t i = 0; i < in->get_num_rows(); i++)
			in_sum += *(double *) in->get(i, k);
		double out_sum = 0;
		for (size_t i = 0; i < mat->get_num_cols(); i++)
			out_sum += *(double *) out->get(i, k);
		printf("%ld: sum of input: %lf, sum of product: %lf\n",
				k, in_sum, out_sum);
	}
}

void print_usage()
{
	fprintf(stderr, "test conf_file matrix_file index_file [options]\n");
	fprintf(stderr,
			"-w matrix_width: the number of columns of the dense matrix\n");
	fprintf(stderr, "-o exec_order: hilbert or seq\n");
	fprintf(stderr, "-c cache_size: cpu cache size\n");
	fprintf(stderr, "-m: force to run in memory\n");
	fprintf(stderr, "-r number: the number of repeats\n");
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		print_usage();
		exit(1);
	}

	size_t mat_width = 0;
	std::string exec_order = "hilbert";
	size_t cpu_cache_size = 1024 * 1024;
	int opt;
	bool in_mem = false;
	size_t repeats = 1;
	while ((opt = getopt(argc, argv, "w:o:c:mr:")) != -1) {
		switch (opt) {
			case 'w':
				mat_width = atoi(optarg);
				break;
			case 'o':
				exec_order = optarg;
				break;
			case 'c':
				cpu_cache_size = atoi(optarg);
				break;
			case 'm':
				in_mem = true;
				break;
			case 'r':
				repeats = atoi(optarg);
				break;
			default:
				print_usage();
				abort();
		}
	}

	std::string conf_file = argv[argc - 3];
	std::string matrix_file = argv[argc - 2];
	std::string index_file = argv[argc - 1];
	signal(SIGINT, int_handler);

	if (exec_order == "seq")
		matrix_conf.set_hilbert_order(false);
	matrix_conf.set_cpu_cache_size(cpu_cache_size);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	SpM_2d_index::ptr index;
	safs::safs_file idx_f(safs::get_sys_RAID_conf(), index_file);
	if (idx_f.exist())
		index = SpM_2d_index::safs_load(index_file);
	else
		index = SpM_2d_index::load(index_file);
	printf("load the matrix index\n");

	sparse_matrix::ptr mat;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), matrix_file);
	if (mat_f.exist() && in_mem)
		mat = sparse_matrix::create(index,
				SpM_2d_storage::safs_load(matrix_file, index));
	else if (mat_f.exist())
		mat = sparse_matrix::create(index, safs::create_io_factory(
					matrix_file, safs::REMOTE_ACCESS));
	else
		mat = sparse_matrix::create(index,
				SpM_2d_storage::load(matrix_file, index));
	printf("load the matrix image\n");

	if (mat_width == 0) {
		for (size_t k = 0; k < repeats; k++)
			test_SpMV(mat);
		for (size_t i = 1; i <= 16; i *= 2)
			for (size_t k = 0; k < repeats; k++)
				test_SpMM(mat, i);
	}
	else {
		for (size_t k = 0; k < repeats; k++)
			test_SpMM(mat, mat_width);
	}

	destroy_flash_matrix();
}
