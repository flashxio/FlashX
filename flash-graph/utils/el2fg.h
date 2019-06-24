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

#ifndef __FG_EL2FG_H_
#define __FG_EL2FG_H_

#include <libgen.h>

#include "common.h"
#include "config_map.h"
#include "native_file.h"
#include "io_interface.h"
#include "utils.h"

using namespace safs;

namespace fg { namespace utils {
struct str2int_pair {
	std::string str;
	int number;
};

static struct str2int_pair edge_type_map[] = {
	{"count", utils::EDGE_COUNT},
	{"timestamp", utils::EDGE_TIMESTAMP},
};

static int type_map_size = sizeof(edge_type_map) / sizeof(edge_type_map[0]);

static inline int conv_edge_type_str2int(const std::string &type_str)
{
	for (int i = 0; i < type_map_size; i++) {
		if (edge_type_map[i].str == type_str) {
			return edge_type_map[i].number;
		}
	}
	return utils::DEFAULT_TYPE;
}

static bool check_graph = false;

 void el2fg(std::vector<std::string> edge_list_files,
         std::string adjacency_list_file,
        std::string index_file, bool directed=true,
        int num_threads=1,
        char* type_str=NULL, bool merge_graph=false,
        bool write_graph=true, bool on_disk=false,
        size_t sort_buf_size=0, size_t write_buf_size=0,
        std::string conf_file="", std::string work_dir=".") {

	utils::set_num_threads(num_threads);
	if (sort_buf_size > 0)
		utils::set_sort_buf_size(sort_buf_size);
	if (write_buf_size > 0)
		utils::set_write_buf_size(write_buf_size);

	int edge_attr_type = utils::DEFAULT_TYPE;
	if (type_str) {
		edge_attr_type = conv_edge_type_str2int(type_str);
	}

	utils::large_io_creator::ptr creator;
	if (conf_file.empty())
		creator = utils::large_io_creator::create(false, work_dir);
	else {
		config_map::ptr configs = config_map::create(conf_file);
		configs->add_options("writable=1");
		safs::init_io_system(configs);
		creator = utils::large_io_creator::create(true, ".");
	}
	if (merge_graph) {
		utils::edge_graph::ptr edge_g = utils::parse_edge_lists(edge_list_files,
				edge_attr_type, directed, !on_disk);
		if (edge_g) {
			utils::disk_serial_graph::ptr g
				= std::static_pointer_cast<
                utils::disk_serial_graph, utils::serial_graph>(
						utils::construct_graph(edge_g, creator));
			// Write the constructed individual graph to a file.
			if (write_graph) {
                if (file_exist(adjacency_list_file))
                    throw std::runtime_error(std::string("Cannot write ") +
                                std::string(adjacency_list_file) +
                                std::string(" -- already exists"));

                if (file_exist(index_file))
                    throw std::runtime_error(std::string("Cannot write ") +
                                std::string(index_file) +
                                std::string(" -- already exists"));

				g->dump(index_file, adjacency_list_file, true);
			}
			printf("There are %ld vertices, %ld non-empty vertices"
                    " and %ld edges\n",
					g->get_num_vertices(), g->get_num_non_empty_vertices(),
					g->get_num_edges());
			if (fg::utils::check_graph) {
				struct timeval start, end;
				gettimeofday(&start, NULL);
				g->check_ext_graph(*edge_g, index_file,
						creator->create_reader(adjacency_list_file));
				gettimeofday(&end, NULL);
				printf("verifying a graph takes %.2f seconds\n",
						time_diff(start, end));
			}
		}
	}
	else {
		std::vector<std::string> graph_files;
		std::vector<std::string> index_files;
		if (edge_list_files.size() > 1) {
			for (size_t i = 0; i < edge_list_files.size(); i++) {
				graph_files.push_back(adjacency_list_file + "-" + itoa(i));
				index_files.push_back(index_file + "-" + itoa(i));
			}
		}
		else {
			graph_files.push_back(adjacency_list_file);
			index_files.push_back(index_file);
		}

		for (size_t i = 0; i < edge_list_files.size(); i++) {
			// construct individual graphs.
			std::vector<std::string> files(1);
			files[0] = edge_list_files[i];

			utils::edge_graph::ptr edge_g = utils::parse_edge_lists(files,
					edge_attr_type, directed, !on_disk);
			if (edge_g == NULL)
				continue;
			utils::disk_serial_graph::ptr g
				= std::static_pointer_cast<
                utils::disk_serial_graph, utils::serial_graph>(
						utils::construct_graph(edge_g, creator));
			// Write the constructed individual graph to a file.
			if (write_graph) {
                if (file_exist(graph_files[i]))
                    throw std::runtime_error(std::string("Cannot write ") +
                                std::string(graph_files[i]) +
                                std::string(" -- already exists"));

                if (file_exist(index_files[i]))
                    throw std::runtime_error(std::string("Cannot write ") +
                                std::string(index_files[i]) +
                                std::string(" -- already exists"));
				g->dump(index_files[i], graph_files[i], true);
			}
			printf("There are %ld vertices, %ld non-empty "
                    "vertices and %ld edges\n",
					g->get_num_vertices(), g->get_num_non_empty_vertices(),
					g->get_num_edges());
			if (fg::utils::check_graph) {
				struct timeval start, end;
				gettimeofday(&start, NULL);
				g->check_ext_graph(*edge_g, index_files[i],
						creator->create_reader(graph_files[i]));
				gettimeofday(&end, NULL);
				printf("verifying a graph takes %.2f seconds\n",
						time_diff(start, end));
			}
		}
	}

	if (is_safs_init())
		destroy_io_system();
 }
} } // End namespace fg::utils
#endif
