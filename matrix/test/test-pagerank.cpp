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
#include "data_frame.h"

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

class mat_init_one: public type_set_operate<float>
{
public:
	virtual void set(float *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1;
	}
};

class mat_init_operate: public type_set_operate<float>
{
	size_t num_vertices;
public:
	mat_init_operate(size_t num_vertices) {
		this->num_vertices = num_vertices;
	}

	virtual void set(float *arr, size_t num_eles, off_t row_idx,
			            off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1.0 / num_vertices;
	}
};

detail::matrix_store::ptr get_out_degree(sparse_matrix::ptr mat, bool in_mem)
{
	int num_nodes = matrix_conf.get_num_nodes();
	detail::matrix_store::ptr out_deg = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, in_mem);
	detail::matrix_store::ptr one = detail::matrix_store::create(
			mat->get_num_cols(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, true);
	one->set_data(mat_init_one());
	mat->multiply<float, float>(one, out_deg);
	return out_deg;
}

dense_matrix::ptr test_PageRank(sparse_matrix::ptr mat, size_t max_niters,
		float damping_factor, size_t num_in_mem)
{
	int num_nodes = matrix_conf.get_num_nodes();
	assert(num_in_mem >= 1);
	// If we allow 3 or more vectors in memory, we will keep
	// the out-degree vector in memory.
	dense_matrix::ptr out_deg = dense_matrix::create(get_out_degree(mat,
				num_in_mem >= 3));
	dense_matrix::ptr one = dense_matrix::create_const<float>(1,
			out_deg->get_num_rows(), out_deg->get_num_cols(),
			matrix_layout_t::L_ROW);
	out_deg = out_deg->pmax(*one);

	detail::matrix_store::ptr tmp = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			get_scalar_type<float>(), num_nodes, num_in_mem >= 2);
	tmp->set_data(mat_init_operate(tmp->get_num_rows()));
	dense_matrix::ptr pr1 = dense_matrix::create(tmp);

	detail::matrix_store::ptr out = detail::matrix_store::create(
			mat->get_num_rows(), 1, matrix_layout_t::L_ROW,
			// If we allow 2 or more vectors in memory, we will keep
			// the out-degree vector in memory.
			get_scalar_type<float>(), num_nodes, num_in_mem >= 2);

	size_t converge = 0;
	size_t num_iters = 0;
	const size_t num_vertices = mat->get_num_rows();
	mat = mat->transpose();
	printf("after transpose\n");
	const float eps = 0.001 / num_vertices;
	struct timeval start, end;
	while (converge < num_vertices && num_iters < max_niters) {
		struct timeval start1, end1;
		gettimeofday(&start, NULL);
		start1 = start;
		dense_matrix::ptr in = pr1->div(*out_deg);
		in = in->cast_ele_type(get_scalar_type<float>());
		// This is guaranteed to be in memory.
		in->move_store(true, num_nodes);
		pr1 = in->multiply_ele(*out_deg);
		gettimeofday(&end1, NULL);
		printf("generate input takes %.3f seconds\n", time_diff(start1, end1));
		start1 = end1;
		mat->multiply<float, float>(in->get_raw_store(), out);
		gettimeofday(&end1, NULL);
		printf("SpMM takes %.3f seconds\n", time_diff(start1, end1));
		start1 = end1;
		// pr2 is a virtual matrix, all computation below has to be computed
		// twice. Once to test the converged vertices, once to compute pr1.
		// TODO we may instruct to materialize pr2 directly.
		dense_matrix::ptr pr2 = dense_matrix::create(out);
		pr2 = pr2->multiply_scalar<float>(damping_factor);
		pr2 = pr2->add_scalar<float>((1.0 - damping_factor) / num_vertices);
		dense_matrix::ptr diff = pr1->minus(*pr2)->abs();
		gettimeofday(&end1, NULL);
		printf("virtual computation takes %.3f seconds\n", time_diff(start1, end1));
		scalar_variable::ptr convg = diff->lt_scalar<float>(eps)->sum();
		assert(convg->get_type() == get_scalar_type<size_t>());
		converge = *(const size_t *) convg->get_raw();
		pr1 = pr2;
		num_iters++;
		gettimeofday(&end, NULL);
		printf("computing converge takes %.3f seconds\n", time_diff(start1, end));
		printf("iter %ld takes %.3f seconds, not converged: %ld\n",
				num_iters, time_diff(start, end), num_vertices - converge);
	}

	return pr1;
}

sparse_matrix::ptr load_2d_matrix(const std::string &matrix_file,
		const std::string &index_file, const std::string &t_matrix_file,
		const std::string &t_index_file, bool in_mem)
{
	SpM_2d_index::ptr index;
	{
		safs::safs_file idx_f(safs::get_sys_RAID_conf(), index_file);
		if (idx_f.exist())
			index = SpM_2d_index::safs_load(index_file);
		else
			index = SpM_2d_index::load(index_file);
	}

	SpM_2d_index::ptr t_index;
	{
		safs::safs_file idx_f(safs::get_sys_RAID_conf(), t_index_file);
		if (idx_f.exist())
			t_index = SpM_2d_index::safs_load(t_index_file);
		else
			t_index = SpM_2d_index::load(t_index_file);
	}
	printf("load the matrix index\n");

	assert(index && t_index);
	sparse_matrix::ptr mat;
	safs::safs_file mat_f(safs::get_sys_RAID_conf(), matrix_file);
	safs::safs_file t_mat_f(safs::get_sys_RAID_conf(), t_matrix_file);
	assert(mat_f.exist() == t_mat_f.exist());
	if (mat_f.exist() && in_mem)
		mat = sparse_matrix::create(index,
				SpM_2d_storage::safs_load(matrix_file, index), t_index,
				SpM_2d_storage::safs_load(t_matrix_file, t_index));
	else if (mat_f.exist())
		mat = sparse_matrix::create(index, safs::create_io_factory(
					matrix_file, safs::REMOTE_ACCESS), t_index,
				safs::create_io_factory(t_matrix_file, safs::REMOTE_ACCESS));
	else
		mat = sparse_matrix::create(index,
				SpM_2d_storage::load(matrix_file, index), t_index,
				SpM_2d_storage::load(t_matrix_file, t_index));
	printf("load the matrix image\n");

	return mat;
}

void print_usage()
{
	fprintf(stderr, "test-pagerank conf_file matrix_file index_file t_matrix_file t_index_file [options]\n");
	fprintf(stderr, "-n number: the number of vectors in memory\n");
	fprintf(stderr, "-i number: the max number of iterations\n");
	fprintf(stderr, "-m: force to run in memory\n");
}

int main(int argc, char *argv[])
{
	int opt;
	size_t num_in_mem = std::numeric_limits<int>::max();
	size_t max_niters = std::numeric_limits<int>::max();
	bool in_mem = false;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "n:i:m")) != -1) {
		num_opts++;
		switch (opt) {
			case 'n':
				num_in_mem = atoi(optarg);
				num_opts++;
				break;
			case 'i':
				max_niters = atoi(optarg);
				num_opts++;
				break;
			case 'm':
				in_mem = true;
				break;
			default:
				print_usage();
				abort();
		}
	}

	argc -= num_opts + 1;
	argv += num_opts + 1;
	if (argc < 5) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string matrix_file = argv[1];
	std::string index_file = argv[2];
	std::string t_matrix_file = argv[3];
	std::string t_index_file = argv[4];
	signal(SIGINT, int_handler);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	{
		sparse_matrix::ptr mat = load_2d_matrix(matrix_file, index_file,
				t_matrix_file, t_index_file, in_mem);
		struct timeval start, end;
		gettimeofday(&start, NULL);
		dense_matrix::ptr pr = test_PageRank(mat, max_niters, 0.85, num_in_mem);
		gettimeofday(&end, NULL);
		printf("PageRank takes %.3f seconds\n", time_diff(start, end));
		pr->move_store(true, -1);
		data_frame::ptr sorted = pr->get_col(0)->sort_with_index();
		std::vector<off_t> first10(10);
		size_t num_vertices = mat->get_num_rows();
		for (size_t i = 0; i < first10.size(); i++)
			first10[i] = num_vertices - i - 1;
		vector::ptr vids = vector::create(sorted->get_vec("idx"))->get(first10);
		if (vids == NULL)
			return -1;
		vector::ptr prs = vector::create(sorted->get_vec("val"))->get(first10);
		if (prs == NULL)
			return -1;

		std::vector<off_t> std_vids = vids->conv2std<off_t>();
		std::vector<float> std_prs = prs->conv2std<float>();
		for (size_t i = 0; i < std_vids.size(); i++)
			printf("%ld: %f\n", std_vids[i], std_prs[i]);
	}

	destroy_flash_matrix();
}
