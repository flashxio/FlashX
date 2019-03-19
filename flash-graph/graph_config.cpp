/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include <math.h>

#include "thread.h"
#include "common.h"

#include "parameters.h"

#include "graph_config.h"

namespace fg
{

void graph_config::print_help()
{
	printf("Configuration parameters in graph algorithm.\n");
	printf("\tthreads: the number of threads processing the graph\n");
	printf("\tprof_file: the output file containing CPU profiling\n");
	printf("\ttrace_file: log IO requests\n");
	printf("\tmax_processing_vertices: the max number of vertices being processed\n");
	printf("\tenable_elevator: enable the elevator algorithm for scheduling vertices\n");
	printf("\tpart_range_size_log: the log2 of the range size in range partitioning\n");
	printf("\tpreload: preload the graph data to the page cache\n");
	printf("\tindex_file_weight: the weight for the graph index file\n");
	printf("\tin_mem_graph: indicate whether to load the entire graph to memory in advance\n");
	printf("\tnum_vparts: the number of vertical partitions\n");
	printf("\tmin_vpart_degree: the min degree of a vertex to perform vertical partitioning\n");
	printf("\tserial_run: run the user code on a vertex in serial\n");
	printf("\tvertex_merge_gap: the gap size allowed when merging two vertex requests\n");
}

void graph_config::print()
{
	std::cout << "Configuration parameters in graph algorithm.\n";
	std::cout << "\tthreads: " << num_threads << std::endl;
	std::cout << "\tprof_file: " << prof_file << std::endl;
	std::cout << "\ttrace_file: " << trace_file << std::endl;
	std::cout << "\tmax_processing_vertices: " <<
        max_processing_vertices << std::endl;
	std::cout << "\tenable_elevator: " << enable_elevator << std::endl;
	std::cout << "\tpart_range_size_log: " << part_range_size_log << std::endl;
	std::cout << "\tpreload: " << _preload << std::endl;
	std::cout << "\tindex_file_weight: " << index_file_weight << std::endl;
	std::cout << "\tin_mem_graph: " << _in_mem_graph << std::endl;
	std::cout << "\tnum_vparts: " << num_vparts << std::endl;
	std::cout << "\tmin_vpart_degree: " << min_vpart_degree << std::endl;
	std::cout << "\tserial_run: " << serial_run << std::endl;
	std::cout << "\tvertex_merge_gap: " << vertex_merge_gap << std::endl;
}

void graph_config::init(config_map::ptr map)
{
	if (map->has_option("threads"))
		map->read_option_int("threads", num_threads);
#ifdef USE_HWLOC
	else {
		num_threads = cpus.get_num_cores();
		num_threads = 1 << (int) ceil(log2(num_threads));
	}
#endif
	std::cout << "FlashGraph runs on " << num_threads << " threads and " <<
        safs::params.get_num_nodes() <<" nodes\n";
	if (!power2(num_threads))
		throw conf_exception("The number of worker threads has to be 2^n");
	map->read_option("prof_file", prof_file);
	map->read_option("trace_file", trace_file);
	map->read_option_int("max_processing_vertices", max_processing_vertices);
	map->read_option_bool("enable_elevator", enable_elevator);
	map->read_option_int("part_range_size_log", part_range_size_log);
	map->read_option_bool("preload", _preload);
	map->read_option_int("index_file_weight", index_file_weight);
	map->read_option_bool("in_mem_graph", _in_mem_graph);
	map->read_option_int("num_vparts", num_vparts);
	map->read_option_int("min_vpart_degree", min_vpart_degree);
	map->read_option_bool("serial_run", serial_run);
	map->read_option_int("vertex_merge_gap", vertex_merge_gap);
}

}
