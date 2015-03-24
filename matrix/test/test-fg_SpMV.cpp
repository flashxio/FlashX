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
	for (size_t i = 0; i < in->get_size(); i++)
		in->set(i, i);
	printf("sum of input: %lf\n", in->sum());

	gettimeofday(&start, NULL);
	fg::FG_vector<double>::ptr fg_out = fg::FG_vector<double>::create(
			fg_m->get_num_rows());
	gettimeofday(&end, NULL);
	printf("initialize FG_vector of %ld entries takes %.3f seconds\n",
			fg_out->get_size(), time_diff(start, end));

	gettimeofday(&start, NULL);
	fg_m->multiply<double>(*in, *fg_out);
	gettimeofday(&end, NULL);
	printf("sum of input: %lf, sum of FG product: %lf, it takes %.3f seconds\n",
			in->sum(), fg_out->sum(), time_diff(start, end));

	NUMA_vector::ptr in_vec = NUMA_vector::create(fg_m->get_num_cols(),
			get_scalar_type<double>());
	for (size_t i = 0; i < fg_m->get_num_cols(); i++)
		in_vec->set(i, i);
	sparse_matrix::ptr m = sparse_matrix::create(fg);
	gettimeofday(&start, NULL);
	NUMA_vector::ptr out = m->multiply<double>(in_vec);
	gettimeofday(&end, NULL);
	double sum = 0;
	for (size_t i = 0; i < fg_m->get_num_cols(); i++)
		sum += out->get<double>(i);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in->sum(), sum, time_diff(start, end));
	destroy_flash_matrix();
}
