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

#include "log.h"
#include "config_map.h"
#include "graph_exception.h"

/**
 * The data structure contains the configurations for the graph engine.
 */
class graph_config
{
	int num_threads;
	std::string prof_file;
	std::string trace_file;
	int max_processing_vertices;
	bool enable_elevator;
	int part_range_size_log;
	bool _preload;
	int index_file_weight;
	bool _in_mem_index;
	bool _in_mem_graph;
	int num_vparts;
	int min_vpart_degree;
	bool serial_run;
public:
	/**
	 * \brief The default constructor that set all configurations to
	 * their default values.
	 */
	graph_config() {
		num_threads = 4;
		max_processing_vertices = 2000;
		enable_elevator = false;
		part_range_size_log = 10;
		_preload = false;
		index_file_weight = 10;
		_in_mem_index = false;
		_in_mem_graph = false;
		num_vparts = 1;
		min_vpart_degree = std::numeric_limits<int>::max();
		serial_run = false;
	}

	/**
	 * \brief Print the explanations of all configurations.
	 */
	void print_help();
	/**
	 * \brief Print the current values of all configurations.
	 */
	void print();

	/**
	 * \brief Set the configurations to the user-defined values.
	 */
	void init(config_map::ptr map);

	/**
	 * \brief Get the output file containing CPU profiling.
	 * \return the file name.
	 */
	const std::string &get_prof_file() const {
		return prof_file;
	}

	/**
	 * \brief Get the number of worker threads used by the graph engine.
	 * \return The number of worker threads.
	 */
	int get_num_threads() const {
		return num_threads;
	}

	/**
	 * \brief Determine whehter to use the elevator algorithm in accessing
	 * the adjacency lists of vertices from SSDs.
	 * \return true if the elevator algorithm is used.
	 */
	bool get_elevator_enabled() const {
		return enable_elevator;
	}

	/**
	 * \brief Get the I/O trace file that records all I/O requests generated
	 * by the graph engine.
	 * \return the I/O trace file anme.
	 */
	const std::string &get_trace_file() const {
		return trace_file;
	}

	/**
	 * \brief Get the maximal number of vertices being processed by
	 * a worker thread.
	 * \return the maximal number of vertices.
	 */
	int get_max_processing_vertices() const {
		return max_processing_vertices;
	}

	/**
	 * \brief Get the size of a partition range at log scale.
	 * \return the size of a partition range at log scale.
	 */
	int get_part_range_size_log() const {
		return part_range_size_log;
	}

	/**
	 * \brief Determine whether to preload the graph data to the page cache.
	 * \return true if the graph is preloaded; else false.
	 */
	bool preload() const {
		return _preload;
	}

	/**
	 * \brief Get the weight for graph index file.
	 * A SAFS file has weight that is used in the page cache. The pages of
	 * the files with higher weight should are more likely to stay in the
	 * page cache.
	 * \return the weight for graph index file.
	 */
	int get_index_file_weight() const {
		return index_file_weight;
	}

	/**
	 * \brief Determine whether to use in-mem vertex index.
	 * \return true if the graph engine uses in-mem vertex index.
	 */
	bool use_in_mem_index() const {
		return _in_mem_index;
	}

	/**
	 * \brief Determine whether to use in-mem graph data.
	 * \return true if the graph engine loads the entire graph data in memory
	 * in advance.
	 */
	bool use_in_mem_graph() const {
		return _in_mem_graph;
	}

	/**
	 * \brief Determine whether to run the user code on a vertex in serial.
	 * \return true if the graph engine runs the user code on a vertex in serial.
	 */
	bool use_serial_run() const {
		return serial_run;
	}

	/**
	 * \brief Get the number of vertical partitions.
	 * \return The number of vertical partitions.
	 */
	int get_num_vparts() const {
		return num_vparts;
	}

	/**
	 * \brief Get the min degree of a vertex to perform vertical partitioning.
	 * \return The min degree of a vertex to perform vertical partitioning.
	 */
	int get_min_vpart_degree() const {
		return min_vpart_degree;
	}
};

inline void graph_config::print_help()
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
	printf("\tin_mem_index: indicate whether to use in-mem vertex index\n");
	printf("\tin_mem_graph: indicate whether to load the entire graph to memory in advance\n");
	printf("\tnum_vparts: the number of vertical partitions\n");
	printf("\tmin_vpart_degree: the min degree of a vertex to perform vertical partitioning\n");
	printf("\tserial_run: run the user code on a vertex in serial\n");
}

inline void graph_config::print()
{
	BOOST_LOG_TRIVIAL(info) << "Configuration parameters in graph algorithm.";
	BOOST_LOG_TRIVIAL(info) << "\tthreads: " << num_threads;
	BOOST_LOG_TRIVIAL(info) << "\tprof_file: " << prof_file;
	BOOST_LOG_TRIVIAL(info) << "\ttrace_file: " << trace_file;
	BOOST_LOG_TRIVIAL(info) << "\tmax_processing_vertices: " << max_processing_vertices;
	BOOST_LOG_TRIVIAL(info) << "\tenable_elevator: " << enable_elevator;
	BOOST_LOG_TRIVIAL(info) << "\tpart_range_size_log: " << part_range_size_log;
	BOOST_LOG_TRIVIAL(info) << "\tpreload: " << _preload;
	BOOST_LOG_TRIVIAL(info) << "\tindex_file_weight: " << index_file_weight;
	BOOST_LOG_TRIVIAL(info) << "\tin_mem_index: " << _in_mem_index;
	BOOST_LOG_TRIVIAL(info) << "\tin_mem_graph: " << _in_mem_graph;
	BOOST_LOG_TRIVIAL(info) << "\tnum_vparts: " << num_vparts;
	BOOST_LOG_TRIVIAL(info) << "\tmin_vpart_degree: " << min_vpart_degree;
	BOOST_LOG_TRIVIAL(info) << "\tserial_run: " << serial_run;
}

inline void graph_config::init(config_map::ptr map)
{
	map->read_option_int("threads", num_threads);
	if (!power2(num_threads))
		throw conf_exception("The number of worker threads has to be 2^n");
	map->read_option("prof_file", prof_file);
	map->read_option("trace_file", trace_file);
	map->read_option_int("max_processing_vertices", max_processing_vertices);
	map->read_option_bool("enable_elevator", enable_elevator);
	map->read_option_int("part_range_size_log", part_range_size_log);
	map->read_option_bool("preload", _preload);
	map->read_option_int("index_file_weight", index_file_weight);
	map->read_option_bool("in_mem_index", _in_mem_index);
	map->read_option_bool("in_mem_graph", _in_mem_graph);
	map->read_option_int("num_vparts", num_vparts);
	map->read_option_int("min_vpart_degree", min_vpart_degree);
	map->read_option_bool("serial_run", serial_run);
}

extern graph_config graph_conf;

#endif
