#ifndef __GRAPH_CONFIG_H__
#define __GRAPH_CONFIG_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "parameters.h"

class graph_config
{
	int num_threads;
	std::string prof_file;
	bool print_io_stat;
	std::string trace_file;
public:
	graph_config() {
		num_threads = 4;
		print_io_stat = false;
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

	bool get_print_io_stat() const {
		return print_io_stat;
	}

	const std::string &get_trace_file() const {
		return trace_file;
	}
};

inline void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
	printf("\tprint_io_stat: print the status of IO instances after the program completes\n");
	printf("\ttrace_file: log IO requests\n");
}

inline void graph_config::print()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: %d\n", num_threads);
	printf("\tprint_io_stat: %d\n", print_io_stat);
	printf("\ttrace_file: %s\n", trace_file.c_str());
}

inline void graph_config::init(const config_map &map)
{
	map.read_option_int("threads", num_threads);
	map.read_option("prof_file", prof_file);
	map.read_option_bool("print_io_stat", print_io_stat);
	map.read_option("trace_file", trace_file);
}

extern graph_config graph_conf;

#endif
