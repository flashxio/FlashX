#ifndef __GRAPH_CONFIG_H__
#define __GRAPH_CONFIG_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "parameters.h"

class graph_config
{
	int num_threads;
	std::string prof_file;
	std::string trace_file;
	int max_processing_vertices;
	bool enable_elevator;
	int part_range_size_log;
	bool _preload;
public:
	graph_config() {
		num_threads = 4;
		max_processing_vertices = 2000;
		enable_elevator = false;
		part_range_size_log = 10;
		_preload = false;
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

	bool get_elevator_enabled() const {
		return enable_elevator;
	}

	const std::string &get_trace_file() const {
		return trace_file;
	}

	int get_max_processing_vertices() const {
		return max_processing_vertices;
	}

	int get_part_range_size_log() const {
		return part_range_size_log;
	}

	bool preload() const {
		return _preload;
	}
};

inline void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
	printf("\ttrace_file: log IO requests\n");
	printf("\tmax_processing_vertices: the max number of vertices being processed\n");
	printf("\tenable_elevator: enable the elevator algorithm for scheduling vertices\n");
	printf("\tpart_range_size_log: the log2 of the range size in range partitioning\n");
	printf("\tpreload: preload the graph data to the page cache\n");
}

inline void graph_config::print()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: %d\n", num_threads);
	printf("\ttrace_file: %s\n", trace_file.c_str());
	printf("\tmax_processing_vertices: %d\n", max_processing_vertices);
	printf("\tenable_elevator: %d\n", enable_elevator);
	printf("\tpart_range_size_log: %d\n", part_range_size_log);
	printf("\tpreload: %d\n", _preload);
}

inline void graph_config::init(const config_map &map)
{
	map.read_option_int("threads", num_threads);
	map.read_option("prof_file", prof_file);
	map.read_option("trace_file", trace_file);
	map.read_option_int("max_processing_vertices", max_processing_vertices);
	map.read_option_bool("enable_elevator", enable_elevator);
	map.read_option_int("part_range_size_log", part_range_size_log);
	map.read_option_bool("preload", _preload);
}

extern graph_config graph_conf;

#endif
