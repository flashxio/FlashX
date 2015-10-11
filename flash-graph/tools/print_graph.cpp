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

using namespace fg;

class empty_attr
{
};

template<class AttrType>
void print_vertex(const ext_mem_undirected_vertex &v)
{
	for (size_t i = 0; i < v.get_num_edges(); i++) {
		if (v.has_edge_data())
			std::cout << v.get_id() << "\t" << v.get_neighbor(i) << "\t"
				<< v.get_edge_data<AttrType>(i) << std::endl;
		else
			std::cout << v.get_id() << "\t" << v.get_neighbor(i) << std::endl;
	}
}

template<>
void print_vertex<empty_attr>(const ext_mem_undirected_vertex &v)
{
	for (size_t i = 0; i < v.get_num_edges(); i++)
		std::cout << v.get_id() << "\t" << v.get_neighbor(i) << std::endl;
}

template<class AttrType>
void print_directed_graph(const char *adj_list, vertex_index::ptr index)
{
	in_mem_cdirected_vertex_index::ptr qindex
		= in_mem_cdirected_vertex_index::create(*index);
	for (size_t i = 0; i < index->get_num_vertices(); i++) {
		off_t off = qindex->get_vertex(i).get_in_off();
		int size = qindex->get_in_size(i);
		ext_mem_undirected_vertex *v
			= (ext_mem_undirected_vertex *) (adj_list + off);
		assert(i == v->get_id());
		assert(size == (int) v->get_size());
		if (v->get_num_edges() > 0)
			print_vertex<AttrType>(*v);

		off = qindex->get_vertex(i).get_out_off();
		size = qindex->get_out_size(i);
		v = (ext_mem_undirected_vertex *) (adj_list + off);
		assert(i == v->get_id());
		assert(size == (int) v->get_size());
		if (v->get_num_edges() > 0)
			print_vertex<AttrType>(*v);
	}
}

template<class AttrType>
void print_undirected_graph(const char *adj_list, vertex_index::ptr index)
{
	in_mem_cundirected_vertex_index::ptr qindex
		= in_mem_cundirected_vertex_index::create(*index);
	for (size_t i = 0; i < index->get_num_vertices(); i++) {
		off_t off = qindex->get_vertex(i).get_off();
		int size = qindex->get_size(i);
		ext_mem_undirected_vertex *v
			= (ext_mem_undirected_vertex *) (adj_list + off);
		assert(i == v->get_id());
		assert(size == (int) v->get_size());
		if (v->get_num_edges() > 0)
			print_vertex<AttrType>(*v);
	}
}

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

	safs::native_file adj_file(adj_file_name);
	ssize_t adj_file_size = adj_file.get_size();
	char *adj_list = new char[adj_file_size];
	FILE *f = fopen(adj_file_name.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		return -1;
	}
	size_t ret = fread(adj_list, adj_file_size, 1, f);
	assert(ret == 1);
	fclose(f);

	vertex_index::ptr index = vertex_index::load(index_file_name);
	graph_header *header = (graph_header *) adj_list;
	header->verify();
	if (header->is_directed_graph()) {
		if (type.empty())
			print_directed_graph<empty_attr>(adj_list, index);
		else if (type == "I") {
			assert(sizeof(int) == index->get_graph_header().get_edge_data_size());
			print_directed_graph<int>(adj_list, index);
		}
		else if (type == "L") {
			assert(sizeof(long) == index->get_graph_header().get_edge_data_size());
			print_directed_graph<long>(adj_list, index);
		}
		else if (type == "F") {
			assert(sizeof(float) == index->get_graph_header().get_edge_data_size());
			print_directed_graph<float>(adj_list, index);
		}
		else if (type == "D") {
			assert(sizeof(double) == index->get_graph_header().get_edge_data_size());
			print_directed_graph<double>(adj_list, index);
		}
		else
			fprintf(stderr, "unsupported edge attribute type\n");
	}
	else {
		if (type.empty())
			print_undirected_graph<empty_attr>(adj_list, index);
		else if (type == "I") {
			assert(sizeof(int) == index->get_graph_header().get_edge_data_size());
			print_undirected_graph<int>(adj_list, index);
		}
		else if (type == "L") {
			assert(sizeof(long) == index->get_graph_header().get_edge_data_size());
			print_undirected_graph<long>(adj_list, index);
		}
		else if (type == "F") {
			assert(sizeof(float) == index->get_graph_header().get_edge_data_size());
			print_undirected_graph<float>(adj_list, index);
		}
		else if (type == "D") {
			assert(sizeof(double) == index->get_graph_header().get_edge_data_size());
			print_undirected_graph<double>(adj_list, index);
		}
		else
			fprintf(stderr, "unsupported edge attribute type\n");
	}
}
