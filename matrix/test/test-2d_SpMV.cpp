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

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr,
				"test conf_file matrix_file index_file [t_matrix_file t_index_file]\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string matrix_file = argv[2];
	std::string index_file = argv[3];
	std::string t_matrix_file;
	std::string t_index_file;
	if (argc == 6) {
		t_matrix_file = argv[4];
		t_index_file = argv[5];
	}
	signal(SIGINT, int_handler);

	struct timeval start, end;
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	SpM_2d_index::ptr index = SpM_2d_index::load(index_file);
	SpM_2d_storage::ptr mat_store = SpM_2d_storage::load(matrix_file, index);

	sparse_matrix::ptr mat;
	if (t_matrix_file.empty())
		mat = sparse_matrix::create(index, mat_store);
	else {
		SpM_2d_index::ptr t_index = SpM_2d_index::load(t_index_file);
		mat = sparse_matrix::create(index, mat_store, t_index,
				SpM_2d_storage::load(t_matrix_file, t_index));
	}

	type_mem_vector<double>::ptr in_vec
		= type_mem_vector<double>::create(mat->get_num_cols());
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_vec->set(i, i);
	gettimeofday(&start, NULL);
	type_mem_vector<double>::ptr out = mat->multiply<double>(in_vec);
	double in_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		in_sum += in_vec->get(i);
	double out_sum = 0;
	for (size_t i = 0; i < mat->get_num_cols(); i++)
		out_sum += out->get(i);
	gettimeofday(&end, NULL);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in_sum, out_sum, time_diff(start, end));
	destroy_flash_matrix();
}
