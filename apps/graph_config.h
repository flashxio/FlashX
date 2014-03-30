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
	int max_processing_vertices;
	bool enable_elevator;
public:
	graph_config() {
		num_threads = 4;
		print_io_stat = false;
		max_processing_vertices = 2000;
		enable_elevator = false;
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

	bool get_elevator_enabled() const {
		return enable_elevator;
	}

	const std::string &get_trace_file() const {
		return trace_file;
	}

	int get_max_processing_vertices() const {
		return max_processing_vertices;
	}
};

inline void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
	printf("\tprint_io_stat: print the status of IO instances after the program completes\n");
	printf("\ttrace_file: log IO requests\n");
	printf("\tmax_processing_vertices: the max number of vertices being processed\n");
	printf("\tenable_elevator: enable the elevator algorithm for scheduling vertices\n");
}

inline void graph_config::print()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: %d\n", num_threads);
	printf("\tprint_io_stat: %d\n", print_io_stat);
	printf("\ttrace_file: %s\n", trace_file.c_str());
	printf("\tmax_processing_vertices: %d\n", max_processing_vertices);
	printf("\tenable_elevator: %d\n", enable_elevator);
}

inline void graph_config::init(const config_map &map)
{
	map.read_option_int("threads", num_threads);
	map.read_option("prof_file", prof_file);
	map.read_option_bool("print_io_stat", print_io_stat);
	map.read_option("trace_file", trace_file);
	map.read_option_int("max_processing_vertices", max_processing_vertices);
	map.read_option_bool("enable_elevator", enable_elevator);
}

extern graph_config graph_conf;

#endif
