/*
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

/*
 * This is a test file to print a small graph.
 */

#include <stdio.h>

#include <string>

#include "vertex.h"
#include "native_file.h"
#include "vertex_index.h"
#include "FGlib.h"
#include "fg_utils.h"

using namespace fg;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "print_graph adj_list_file index_file [attr_type]\n");
		return -1;
	}

	const std::string adj_file_name = argv[1];
	const std::string index_file_name = argv[2];
	std::string type;
	if (argc == 4)
		type = argv[3];

	FG_graph::ptr graph = FG_graph::create(adj_file_name, index_file_name,
			config_map::ptr());
	print_graph_el(graph, "\t", type, stdout);
}
