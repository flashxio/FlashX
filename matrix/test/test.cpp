#include "sparse_matrix.h"

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "test conf_file graph_file index_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	std::string graph_file = argv[2];
	std::string index_file = argv[3];

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);
	FG_graph::ptr fg = FG_graph::create(graph_file, index_file, configs);
	sparse_matrix::ptr m = sparse_matrix::create(fg);
	FG_vector<double>::ptr in = FG_vector<double>::create(m->get_num_cols());
	in->init_rand(1000 * 1000);
	FG_vector<double>::ptr out = m->multiply<double>(in);
	printf("sum of product: %lf\n", out->sum());
	destroy_flash_matrix();
}
