#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif
#include "matrix_config.h"
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

#ifdef PROFILER
	if (!fg::graph_conf.get_prof_file().empty())
		ProfilerStart(fg::graph_conf.get_prof_file().c_str());
#endif
	gettimeofday(&start, NULL);
	fg_m->multiply<double>(*in, *fg_out);
	gettimeofday(&end, NULL);
#ifdef PROFILER
	if (!fg::graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	printf("sum of input: %lf, sum of FG product: %lf, it takes %.3f seconds\n",
			in->sum(), fg_out->sum(), time_diff(start, end));

	sparse_matrix::ptr m = sparse_matrix::create(fg);
	detail::NUMA_vec_store::ptr in_vec = detail::NUMA_vec_store::create(
			m->get_num_cols(), matrix_conf.get_num_nodes(),
			get_scalar_type<double>());
	for (size_t i = 0; i < m->get_num_cols(); i++)
		in_vec->set<double>(i, i);
	detail::NUMA_vec_store::ptr out = detail::NUMA_vec_store::create(
			m->get_num_rows(), matrix_conf.get_num_nodes(),
			get_scalar_type<double>());

#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStart(matrix_conf.get_prof_file().c_str());
#endif
	gettimeofday(&start, NULL);
	m->multiply<double>(*in_vec, *out);
	gettimeofday(&end, NULL);
#ifdef PROFILER
	if (!matrix_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	double sum = 0;
	for (size_t i = 0; i < fg_m->get_num_cols(); i++)
		sum += out->get<double>(i);
	printf("sum of input: %lf, sum of product: %lf, it takes %.3f seconds\n",
			in->sum(), sum, time_diff(start, end));
	destroy_flash_matrix();
}
