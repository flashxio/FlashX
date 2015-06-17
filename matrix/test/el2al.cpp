/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include <string.h>
#include <assert.h>

#include <vector>
#include <string>

#include "log.h"
#include "FG_basic_types.h"
#include "vertex.h"
#include "graph_file_header.h"
#include "vertex_index.h"
#include "in_mem_storage.h"

#include "generic_type.h"
#include "data_io.h"
#include "data_frame.h"
#include "vector.h"
#include "vector_vector.h"
#include "fm_utils.h"
#include "matrix_config.h"
#include "sparse_matrix.h"

using namespace fm;

/*
 * This class parses a line into an edge (source, destination).
 */
class edge_parser: public line_parser
{
public:
	size_t parse(const std::vector<std::string> &lines, data_frame &df) const;
	size_t get_num_cols() const {
		return 2;
	}

	std::string get_col_name(off_t idx) const {
		if (idx == 0)
			return "source";
		else
			return "dest";
	}

	const scalar_type &get_col_type(off_t idx) const {
		return get_scalar_type<fg::vertex_id_t>();
	}
};

size_t edge_parser::parse(const std::vector<std::string> &lines,
		data_frame &df) const
{
	detail::smp_vec_store::ptr froms = detail::smp_vec_store::create(lines.size(),
			get_scalar_type<fg::vertex_id_t>());
	detail::smp_vec_store::ptr tos = detail::smp_vec_store::create(lines.size(),
			get_scalar_type<fg::vertex_id_t>());
	for (size_t i = 0; i < lines.size(); i++) {
		const char *line = lines[i].c_str();
		int len = strlen(line);
		const char *first = line;
		// Skip space
		for (; isspace(*first); first++);
		if (*first == '#')
			continue;

		// Make sure we get a number.
		if (!isdigit(*first)) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("the first entry isn't a number: ") + first;
			continue;
		}
		fg::vertex_id_t from = atol(first);
		assert(from >= 0 && from < fg::MAX_VERTEX_ID);

		const char *second = first;
		// Go to the end of the first number.
		for (; isdigit(*second); second++);
		if (second - line == len) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("there isn't second entry: ") + line;
			continue;
		}
		// Skip space between two numbers.
		for (; isspace(*second); second++);
		// Make sure we get a number.
		if (!isdigit(*second)) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("the second entry isn't a number: ") + second;
			continue;
		}
		fg::vertex_id_t to = atol(second);
		assert(to >= 0 && to < fg::MAX_VERTEX_ID);

		froms->set(i, from);
		tos->set(i, to);
	}

	df.get_vec(0)->append(*froms);
	df.get_vec(1)->append(*tos);
	return froms->get_length();
}

void print_usage()
{
	fprintf(stderr, "convert an edge list to adjacency lists\n");
	fprintf(stderr, "el2al [options] conf_file edge_file graph_name\n");
	fprintf(stderr, "-u: undirected graph\n");
}

int main(int argc, char *argv[])
{
	bool directed = true;
	int opt;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "u")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				directed = false;
				break;
			default:
				print_usage();
				exit(1);
		}
	}

	argv += 1 + num_opts;
	argc -= 1 + num_opts;
	if (argc < 2) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string file_name = argv[1];
	std::string graph_name = argv[2];
	std::string adj_file = graph_name + ".adj";
	std::string index_file = graph_name + ".index";
	std::vector<std::string> files;
	files.push_back(file_name);

	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	edge_parser parser;
	data_frame::ptr df = read_lines(files, parser);
	printf("There are %ld edges\n", df->get_num_entries());

	fg::vertex_id_t max_vid = 0;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		vector::ptr vec = vector::create(df->get_vec(i));
		max_vid = std::max(max_vid, vec->max<fg::vertex_id_t>());
	}
	printf("max id: %d\n", max_vid);
	detail::vec_store::ptr seq_vec = detail::create_vec_store<fg::vertex_id_t>(
			0, max_vid, 1);
	detail::vec_store::ptr rep_vec = detail::create_vec_store<fg::vertex_id_t>(
			max_vid + 1, fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());

	fg::FG_graph::ptr graph;
	if (directed) {
		// I artificially add an invalid out-edge for each vertex, so it's
		// guaranteed that each vertex exists in the adjacency lists.
		data_frame::ptr new_df = data_frame::create();
		new_df->add_vec(df->get_vec_name(0), seq_vec);
		new_df->add_vec(df->get_vec_name(1), rep_vec);
		df->append(new_df);

		// I artificially add an invalid in-edge for each vertex.
		new_df = data_frame::create();
		new_df->add_vec(df->get_vec_name(1), seq_vec);
		new_df->add_vec(df->get_vec_name(0), rep_vec);
		df->append(new_df);

		graph = create_fg_graph(graph_name, df, true);
	}
	else {
		// I artificially add an invalid out-edge for each vertex, so it's
		// guaranteed that each vertex exists in the adjacency lists.
		data_frame::ptr new_df = data_frame::create();
		new_df->add_vec(df->get_vec_name(0), seq_vec);
		new_df->add_vec(df->get_vec_name(1), rep_vec);
		df->append(new_df);

		detail::vec_store::ptr vec0 = df->get_vec(0)->deep_copy();
		detail::vec_store::ptr vec1 = df->get_vec(1)->deep_copy();
		vec0->append(df->get_vec_ref(1));
		vec1->append(df->get_vec_ref(0));
		new_df = data_frame::create();
		new_df->add_vec(df->get_vec_name(0), vec0);
		new_df->add_vec(df->get_vec_name(1), vec1);
		df = new_df;

		graph = create_fg_graph(graph_name, df, false);
	}

	if (graph->get_index_data())
		graph->get_index_data()->dump(index_file);
	if (graph->get_graph_data())
		graph->get_graph_data()->dump(adj_file);

	destroy_flash_matrix();

	return 0;
}
