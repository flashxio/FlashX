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
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "test conf_file graph_file index_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];
	signal(SIGINT, int_handler);

	struct timeval start, end;
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	fg::FG_graph::ptr fg = fg::FG_graph::create(graph_file, index_file, configs);
	fg::FG_adj_matrix::ptr fg_m = fg::FG_adj_matrix::create(fg);
	fg::FG_vector<double>::ptr in = fg::FG_vector<double>::create(fg_m->get_num_cols());
	in->init_rand(1000 * 1000);
	printf("sum of input: %lf\n", in->sum());

	fg::FG_vector<double>::ptr fg_out = fg::FG_vector<double>::create(
			fg_m->get_num_rows());
	gettimeofday(&start, NULL);
	fg_m->multiply<double>(*in, *fg_out);
	gettimeofday(&end, NULL);
	printf("sum of input: %lf, sum of FG product: %lf, it takes %.3f seconds\n",
			in->sum(), fg_out->sum(), time_diff(start, end));

	sparse_matrix::ptr m = sparse_matrix::create(fg);
	gettimeofday(&start, NULL);
	fg::FG_vector<double>::ptr out = m->multiply<double>(in);
	gettimeofday(&end, NULL);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in->sum(), out->sum(), time_diff(start, end));
	destroy_flash_matrix();
}
