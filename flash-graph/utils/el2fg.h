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

 void el2fg_(std::vector<std::string> edge_list_files,
         std::string adjacency_list_file,
        std::string index_file, bool directed=true,
        int num_threads=1,
        char* type_str=NULL, bool merge_graph=false,
        bool write_graph=true, bool on_disk=false,
        size_t sort_buf_size=0, size_t write_buf_size=0,
        std::string conf_file=0, std::string work_dir=".");

void el2fg(std::vector<std::string> extfns,
        std::string safs_adj_fn,
        std::string safs_index_fn, bool directed=true,
        int nthread=4) {
    el2fg_(extfns, safs_adj_fn, safs_index_fn, directed, nthread);
}
} } // End namespace fg::utils
#endif
