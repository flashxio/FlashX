#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include "in_mem_storage.h"
#include "matrix_config.h"
#include "io_interface.h"
#include "safs_file.h"
#include "sparse_matrix.h"
#include "NUMA_dense_matrix.h"
#include "EM_dense_matrix.h"
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

class mat_init_operate: public type_set_operate<double>
{
public:
	mat_init_operate(size_t num_rows, size_t num_cols) {
	}

	virtual void set(double *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			// Reduce the value of the elements to avoid float-point
			// rounding error.
			arr[i] = (row_idx % 100) * (i + col_idx + 1);
	}
};

void test_SpMM(sparse_matrix::ptr mat, size_t mat_width, int num_nodes,
		bool ext_mem_out)
{
	printf("test sparse matrix dense matrix multiplication\n");
	struct timeval start, end;
	detail::mem_matrix_store::ptr in
		= detail::mem_matrix_store::create(mat->get_num_cols(), mat_width,
				matrix_layout_t::L_ROW, get_scalar_type<double>(), num_nodes);
	if (num_nodes < 0) {
		// This forces all memory is allocated in a single NUMA node.
		for (size_t i = 0; i < in->get_num_rows(); i++)
			for (size_t j = 0; j < in->get_num_cols(); j++)
				in->set<double>(i, j, (i % 100) * (j + 1));
	}
	else
		in->set_data(mat_init_operate(in->get_num_rows(), in->get_num_cols()));
	printf("set input data\n");

	// Initialize the output matrix and allocate pages for it.
	detail::matrix_store::ptr out;
	if (ext_mem_out)
		out = detail::EM_matrix_store::create(mat->get_num_rows(), mat_width,
				matrix_layout_t::L_ROW, get_scalar_type<double>());
	else
		out = detail::mem_matrix_store::create(mat->get_num_rows(), mat_width,
				matrix_layout_t::L_ROW, get_scalar_type<double>(), num_nodes);
	if (num_nodes < 0 && out->is_in_mem()) {
		detail::mem_matrix_store::ptr mem_out
			= detail::mem_matrix_store::cast(out);
		// This forces all memory is allocated in a single NUMA node.
		for (size_t i = 0; i < in->get_num_rows(); i++)
			for (size_t j = 0; j < in->get_num_cols(); j++)
				mem_out->set<double>(i, j, 0);
	}
	else if (out->is_in_mem())
		out->reset_data();
	printf("reset output data\n");
	printf("in mat is on %d nodes and out mat is on %d nodes\n",
			in->get_num_nodes(), out->get_num_nodes());

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	printf("Start SpMM\n");
	gettimeofday(&start, NULL);
	mat->multiply<double, float>(in, out);
	gettimeofday(&end, NULL);
	printf("SpMM completes\n");
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("it takes %.3f seconds\n", time_diff(start, end));

	dense_matrix::ptr in_mat = dense_matrix::create(in);
	dense_matrix::ptr out_mat = dense_matrix::create(out);
	std::vector<double> in_col_sum = in_mat->col_sum()->conv2std<double>();
	std::vector<double> out_col_sum = out_mat->col_sum()->conv2std<double>();
	for (size_t k = 0; k < in->get_num_cols(); k++) {
		printf("%ld: sum of input: %lf, sum of product: %lf\n",
				k, in_col_sum[k], out_col_sum[k]);
	}
}

sparse_matrix::ptr load_2d_matrix(const std::string &matrix_file,
		const std::string &index_file, bool in_mem)
{
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

	return mat;
}

sparse_matrix::ptr load_fg_matrix(const std::string &matrix_file,
		const std::string &index_file, bool in_mem, config_map::ptr configs,
		const std::string &entry_type)
{
	fg::FG_graph::ptr fg = fg::FG_graph::create(matrix_file, index_file, configs);
	if (entry_type.empty())
		return sparse_matrix::create(fg, NULL);
	else if (entry_type == "I")
		return sparse_matrix::create(fg, &get_scalar_type<int>());
	else if (entry_type == "L")
		return sparse_matrix::create(fg, &get_scalar_type<long>());
	else if (entry_type == "F")
		return sparse_matrix::create(fg, &get_scalar_type<float>());
	else if (entry_type == "D")
		return sparse_matrix::create(fg, &get_scalar_type<double>());
	else {
		fprintf(stderr, "unknown entry type\n");
		return sparse_matrix::ptr();
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
	fprintf(stderr, "-g: the matrix is stored in FlashGraph format\n");
	fprintf(stderr, "-n number: the number of NUMA nodes\n");
	fprintf(stderr, "-t type: the type of non-zero entries\n");
	fprintf(stderr, "-e: output matrix in external memory\n");
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
	bool use_fg = false;
	int num_nodes = 0;
	bool ext_mem_out = false;
	std::string entry_type;
	while ((opt = getopt(argc, argv, "w:o:c:mr:gn:t:e")) != -1) {
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
			case 'g':
				use_fg = true;
				break;
			case 'n':
				num_nodes = atoi(optarg);
				break;
			case 't':
				entry_type = optarg;
				break;
			case 'e':
				ext_mem_out = true;
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

	if (ext_mem_out && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS must be enabled if output matrix is in external memory\n");
		return -1;
	}

	sparse_matrix::ptr mat;
	if (use_fg)
		mat = load_fg_matrix(matrix_file, index_file, in_mem, configs,
				entry_type);
	else
		mat = load_2d_matrix(matrix_file, index_file, in_mem);

	if (num_nodes == 0)
		num_nodes = matrix_conf.get_num_nodes();

	if (mat_width == 0) {
		for (size_t i = 1; i <= 16; i *= 2)
			for (size_t k = 0; k < repeats; k++)
				test_SpMM(mat, i, num_nodes, ext_mem_out);
	}
	else {
		for (size_t k = 0; k < repeats; k++)
			test_SpMM(mat, mat_width, num_nodes, ext_mem_out);
	}

	destroy_flash_matrix();
}
