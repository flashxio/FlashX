#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include "io_interface.h"
#include "safs_file.h"

#include "matrix_config.h"
#include "data_frame.h"

#include "matrix_algs.h"

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

sparse_matrix::ptr load_2d_matrix(const std::string &matrix_file,
		const std::string &index_file, const std::string &t_matrix_file,
		const std::string &t_index_file, bool in_mem)
{
	SpM_2d_index::ptr index;
	if (!index_file.empty()) {
		safs::safs_file idx_f(safs::get_sys_RAID_conf(), index_file);
		if (idx_f.exist())
			index = SpM_2d_index::safs_load(index_file);
		else
			index = SpM_2d_index::load(index_file);
	}

	SpM_2d_index::ptr t_index;
	if (!t_index_file.empty()) {
		safs::safs_file idx_f(safs::get_sys_RAID_conf(), t_index_file);
		if (idx_f.exist())
			t_index = SpM_2d_index::safs_load(t_index_file);
		else
			t_index = SpM_2d_index::load(t_index_file);
	}
	printf("load the matrix index\n");

	sparse_matrix::ptr mat;
	if (!matrix_file.empty() && !t_matrix_file.empty()) {
		assert(index && t_index);
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
	}
	else {
		assert(index);
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
	}
	printf("load the matrix image\n");

	return mat;
}

void test_pagerank(sparse_matrix::ptr mat, size_t max_niters, size_t num_in_mem)
{
	struct timeval start, end;
	gettimeofday(&start, NULL);
	dense_matrix::ptr pr = alg::PageRank(mat, max_niters, 0.85, num_in_mem);
	gettimeofday(&end, NULL);
	printf("PageRank takes %.3f seconds\n", time_diff(start, end));
	pr->move_store(true, -1);
#if 0
	data_frame::ptr sorted = pr->get_col(0)->sort_with_index();
	std::vector<off_t> first10(10);
	size_t num_vertices = mat->get_num_rows();
	for (size_t i = 0; i < first10.size(); i++)
		first10[i] = num_vertices - i - 1;
	vector::ptr vids = vector::create(sorted->get_vec("idx"))->get(first10);
	if (vids == NULL) {
		fprintf(stderr, "can't get index vector from the sort result\n");
		return;
	}
	vector::ptr prs = vector::create(sorted->get_vec("val"))->get(first10);
	if (prs == NULL) {
		fprintf(stderr, "can't get value vector from the sort result\n");
		return;
	}

	std::vector<off_t> std_vids = vids->conv2std<off_t>();
	std::vector<float> std_prs = prs->conv2std<float>();
	for (size_t i = 0; i < std_vids.size(); i++)
		printf("%ld: %f\n", std_vids[i], std_prs[i]);
#endif
}

void test_nmf(sparse_matrix::ptr mat, size_t k, size_t max_niters,
		size_t num_in_mem)
{
	auto res = fm::alg::NMF(mat, k, max_niters, num_in_mem);
}

void print_usage()
{
	fprintf(stderr, "test-algs [options] conf_file matrix_file index_file [t_matrix_file t_index_file]\n");
	fprintf(stderr, "-k rank: the rank of factorized matrices in NMF\n");
	fprintf(stderr, "-a name: test algorithm\n");
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
	int rank = 1;
	std::string alg;
	while ((opt = getopt(argc, argv, "k:n:i:ma:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'k':
				rank = atoi(optarg);
				num_opts++;
				break;
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
			case 'a':
				alg = optarg;
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	argc -= num_opts + 1;
	argv += num_opts + 1;
	if (argc < 3) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string matrix_file = argv[1];
	std::string index_file = argv[2];
	std::string t_matrix_file;
	std::string t_index_file;
	if (argc == 5) {
		t_matrix_file = argv[3];
		t_index_file = argv[4];
	}
	signal(SIGINT, int_handler);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	sparse_matrix::ptr mat = load_2d_matrix(matrix_file, index_file,
			t_matrix_file, t_index_file, in_mem);
	if (alg == "pagerank")
		test_pagerank(mat, max_niters, num_in_mem);
	else if (alg == "nmf")
		test_nmf(mat, rank, max_niters, num_in_mem);
	else {
		fprintf(stderr, "unknown algorithm %s\n", alg.c_str());
	}
	mat = NULL;

	destroy_flash_matrix();
}
