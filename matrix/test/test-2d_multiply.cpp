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

bool spmm_debug = false;

void int_handler(int sig_num)
{
#ifdef PROFILER
	printf("stop profiling\n");
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

typedef double mat_ele_t;

class mat_init_operate: public type_set_operate<mat_ele_t>
{
public:
	mat_init_operate(size_t num_rows, size_t num_cols) {
	}

	virtual void set(mat_ele_t *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			// Reduce the value of the elements to avoid float-point
			// rounding error.
			arr[i] = (row_idx % 100) * (i + col_idx + 1);
	}
};

void test_SpMM(sparse_matrix::ptr mat, size_t mat_width, size_t indiv_mat_width,
		int num_nodes, bool ext_mem, size_t repeats)
{
	struct timeval start, end;
	std::vector<detail::matrix_store::ptr> ins(mat_width / indiv_mat_width);
	for (size_t i = 0; i < ins.size(); i++) {
		detail::matrix_store::ptr in;
		if (ext_mem)
			in = detail::EM_matrix_store::create(mat->get_num_cols(),
					indiv_mat_width, matrix_layout_t::L_ROW,
					get_scalar_type<mat_ele_t>());
		else
			in = detail::mem_matrix_store::create(mat->get_num_cols(),
					indiv_mat_width, matrix_layout_t::L_ROW,
					get_scalar_type<mat_ele_t>(), num_nodes);
		if (!ext_mem && num_nodes < 0) {
			detail::mem_matrix_store::ptr mem_in
				= detail::mem_matrix_store::cast(in);
			// This forces all memory is allocated in a single NUMA node.
			for (size_t i = 0; i < in->get_num_rows(); i++)
				for (size_t j = 0; j < in->get_num_cols(); j++)
					mem_in->set<mat_ele_t>(i, j, (i % 100) * (j + 1));
		}
		else
			in->set_data(mat_init_operate(in->get_num_rows(), in->get_num_cols()));
		ins[i] = in;
	}

	std::vector<detail::matrix_store::ptr> outs(ins.size());
	for (size_t i = 0; i < outs.size(); i++) {
		detail::matrix_store::ptr out;
		// Initialize the output matrix and allocate pages for it.
		if (ext_mem)
			out = detail::EM_matrix_store::create(mat->get_num_rows(),
					indiv_mat_width, matrix_layout_t::L_ROW,
					get_scalar_type<mat_ele_t>());
		else
			out = detail::mem_matrix_store::create(mat->get_num_rows(),
					indiv_mat_width, matrix_layout_t::L_ROW,
					get_scalar_type<mat_ele_t>(), num_nodes);
		if (num_nodes < 0 && out->is_in_mem()) {
			detail::mem_matrix_store::ptr mem_out
				= detail::mem_matrix_store::cast(out);
			// This forces all memory is allocated in a single NUMA node.
			for (size_t i = 0; i < out->get_num_rows(); i++)
				for (size_t j = 0; j < out->get_num_cols(); j++)
					mem_out->set<mat_ele_t>(i, j, 0);
		}
		else if (out->is_in_mem())
			out->reset_data();
		outs[i] = out;
	}
	printf("in mat is on %d nodes and out mat is on %d nodes\n",
			ins[0]->get_num_nodes(), outs[0]->get_num_nodes());

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	printf("Start SpMM\n");
	for (size_t k = 0; k < repeats; k++) {
		gettimeofday(&start, NULL);
		for (size_t i = 0; i < ins.size(); i++)
			mat->multiply<mat_ele_t, float>(ins[i], outs[i]);
		gettimeofday(&end, NULL);
		printf("it takes %.3f seconds\n", time_diff(start, end));
	}
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	if (spmm_debug) {
		for (size_t i = 0; i < ins.size(); i++) {
			dense_matrix::ptr in_mat = dense_matrix::create(ins[i]);
			dense_matrix::ptr out_mat = dense_matrix::create(outs[i]);
			dense_matrix::ptr sum = in_mat->col_sum();
			vector::ptr sum_vec = sum->get_col(0);
			std::vector<mat_ele_t> in_col_sum = sum_vec->conv2std<mat_ele_t>();
			sum = out_mat->col_sum();
			sum_vec = sum->get_col(0);
			std::vector<mat_ele_t> out_col_sum = sum_vec->conv2std<mat_ele_t>();
			for (size_t k = 0; k < in_mat->get_num_cols(); k++) {
				printf("%ld: sum of input: %lf, sum of product: %lf\n",
						k, in_col_sum[k], out_col_sum[k]);
			}
		}
	}
}

sparse_matrix::ptr load_2d_matrix(const std::string &matrix_file,
		const std::string &index_file, bool in_mem)
{
	SpM_2d_index::ptr index;
	if (safs::exist_safs_file(index_file))
		index = SpM_2d_index::safs_load(index_file);
	else
		index = SpM_2d_index::load(index_file);
	printf("load the matrix index\n");

	sparse_matrix::ptr mat;
	if (safs::exist_safs_file(matrix_file) && in_mem)
		mat = sparse_matrix::create(index,
				SpM_2d_storage::safs_load(matrix_file, index));
	else if (safs::exist_safs_file(matrix_file))
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
	fprintf(stderr, "-i width: the number of columns of the individual dense matrix\n");
	fprintf(stderr, "-d: enable debug\n");
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		print_usage();
		exit(1);
	}

	size_t mat_width = 1;
	size_t indiv_mat_width = 0;
	std::string exec_order = "hilbert";
	size_t cpu_cache_size = 1024 * 1024;
	int opt;
	bool in_mem = false;
	size_t repeats = 1;
	bool use_fg = false;
	int num_nodes = 0;
	bool ext_mem = false;
	std::string entry_type;
	while ((opt = getopt(argc, argv, "w:o:c:mr:gn:t:ei:d")) != -1) {
		switch (opt) {
			case 'w':
				mat_width = atoi(optarg);
				break;
			case 'i':
				indiv_mat_width = atoi(optarg);
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
				ext_mem = true;
				break;
			case 'd':
				spmm_debug = true;
				break;
			default:
				print_usage();
				abort();
		}
	}

	if (indiv_mat_width == 0)
		indiv_mat_width = mat_width;

	std::string conf_file = argv[argc - 3];
	std::string matrix_file = argv[argc - 2];
	std::string index_file = argv[argc - 1];
	signal(SIGINT, int_handler);

	if (exec_order == "seq")
		matrix_conf.set_hilbert_order(false);
	matrix_conf.set_cpu_cache_size(cpu_cache_size);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	if (ext_mem && !safs::is_safs_init()) {
		fprintf(stderr,
				"SAFS must be enabled if output matrix is in external memory\n");
		return -1;
	}

	sparse_matrix::ptr mat;
	try {
		if (use_fg)
			mat = load_fg_matrix(matrix_file, index_file, in_mem, configs,
					entry_type);
		else
			mat = load_2d_matrix(matrix_file, index_file, in_mem);
	} catch (std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(-1);
	}

	if (num_nodes == 0)
		num_nodes = matrix_conf.get_num_nodes();

	printf("SpMM on %s with matrix width of %ld\n", matrix_file.c_str(), mat_width);
	test_SpMM(mat, mat_width, indiv_mat_width, num_nodes, ext_mem, repeats);

	destroy_flash_matrix();
}
