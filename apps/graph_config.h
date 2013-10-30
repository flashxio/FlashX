#ifndef __GRAPH_CONFIG_H__
#define __GRAPH_CONFIG_H__

#include "parameters.h"

class graph_config
{
	int num_threads;
	std::string prof_file;
public:
	graph_config() {
		num_threads = 4;
	}

	void print_help();
	void print();

	void init(const config_map &map);

	const std::string &get_prof_file() const {
		return prof_file;
	}

	int get_num_threads() const {
		return num_threads;
	}
};

void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
}

void graph_config::print()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: %d\n", num_threads);
}

void graph_config::init(const config_map &map)
{
	map.read_option_int("threads", num_threads);
	map.read_option("prof_file", prof_file);
}

#endif
