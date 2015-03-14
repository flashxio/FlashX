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

#include "generic_type.h"
#include "data_io.h"
#include "data_frame.h"
#include "mem_data_frame.h"
#include "mem_vector.h"
#include "vector_vector.h"
#include "fm_utils.h"

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
	type_mem_vector<fg::vertex_id_t>::ptr froms
		= type_mem_vector<fg::vertex_id_t>::create(lines.size());
	type_mem_vector<fg::vertex_id_t>::ptr tos
		= type_mem_vector<fg::vertex_id_t>::create(lines.size());
	for (size_t i = 0; i < lines.size(); i++) {
		const char *line = lines[i].c_str();
		int len = strlen(line);
		const char *first = line;
		// Skip space
		for (; isspace(*first); first++);
		// Make sure we get a number.
		if (!isdigit(*first)) {
			BOOST_LOG_TRIVIAL(error)
				<< std::string("the first entry isn't a number: ") + first;
			continue;
		}
		long from = atol(first);
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
		long to = atol(second);
		assert(to >= 0 && to < fg::MAX_VERTEX_ID);

		froms->set(i, from);
		tos->set(i, to);
	}

	df.get_vec(0)->append(*froms);
	df.get_vec(1)->append(*tos);
	return froms->get_length();
}

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "el2al edge_file graph_file\n");
		return -1;
	}

	std::string file_name = argv[1];
	std::string adj_file = std::string(argv[2]) + ".adj";
	std::string index_file = std::string(argv[2]) + ".index";
	std::vector<std::string> files;
	files.push_back(file_name);

	edge_parser parser;
	data_frame::ptr df = read_lines(files, parser);
	size_t num_edges = df->get_num_entries();
	printf("There are %ld edges\n", num_edges);
	fg::vertex_id_t max_vid = 0;
	for (size_t i = 0; i < df->get_num_vecs(); i++) {
		type_mem_vector<fg::vertex_id_t>::ptr vec
			= type_mem_vector<fg::vertex_id_t>::cast(df->get_vec(i));
		max_vid = std::max(max_vid, vec->max());
	}
	printf("max id: %d\n", max_vid);

	vector::ptr seq_vec = create_vector<fg::vertex_id_t>(0, max_vid, 1);
	vector::ptr rep_vec = create_vector<fg::vertex_id_t>(max_vid + 1,
			fg::INVALID_VERTEX_ID);
	assert(seq_vec->get_length() == rep_vec->get_length());
	// I artificially add an invalid out-edge for each vertex, so it's
	// guaranteed that each vertex exists in the adjacency lists.
	mem_data_frame::ptr new_df = mem_data_frame::create();
	new_df->add_vec(parser.get_col_name(0), seq_vec);
	new_df->add_vec(parser.get_col_name(1), rep_vec);
	df->append(new_df);

	// I artificially add an invalid in-edge for each vertex.
	new_df = mem_data_frame::create();
	new_df->add_vec(parser.get_col_name(1), seq_vec);
	new_df->add_vec(parser.get_col_name(0), rep_vec);
	df->append(new_df);

	// For in-edge adjacency lists, all edges share the same destination vertex
	// should be stored together.
	mem_data_frame::ptr tmp = mem_data_frame::create();
	tmp->add_vec("dest", df->get_vec("dest"));
	tmp->add_vec("source", df->get_vec("source"));
	df = tmp;
	vector_vector::ptr in_adjs = fm::create_1d_matrix(df);
	printf("There are %ld in-edge adjacency lists and they use %ld bytes in total\n",
			in_adjs->get_num_vecs(), in_adjs->get_tot_num_entries());

	// For out-edge adjacency lists, all edges share the same source vertex
	// should be stored together.
	tmp = mem_data_frame::create();
	tmp->add_vec("source", df->get_vec("source"));
	tmp->add_vec("dest", df->get_vec("dest"));
	df = tmp;
	vector_vector::ptr out_adjs = create_1d_matrix(df);
	printf("There are %ld out-edge adjacency lists and they use %ld bytes in total\n",
			out_adjs->get_num_vecs(), out_adjs->get_tot_num_entries());

	assert(out_adjs->get_num_vecs() == in_adjs->get_num_vecs());
	size_t num_vertices = out_adjs->get_num_vecs();
	type_mem_vector<fg::vsize_t>::ptr num_in_edges
		= type_mem_vector<fg::vsize_t>::create(num_vertices);
	type_mem_vector<fg::vsize_t>::ptr num_out_edges
		= type_mem_vector<fg::vsize_t>::create(num_vertices);
	for (size_t i = 0; i < num_vertices; i++) {
		num_in_edges->set(i, in_adjs->get_length(i));
		num_out_edges->set(i, out_adjs->get_length(i));
	}
	fg::graph_header header(fg::graph_type::DIRECTED, num_vertices, num_edges, 0);

	// Construct the vertex index.
	// The vectors that contains the numbers of edges have the length of #V + 1
	// because we add -1 to the edge lists artificially and the last entries
	// are the number of vertices.
	fg::cdirected_vertex_index::ptr vindex
		= fg::cdirected_vertex_index::construct(num_vertices,
				(const fg::vsize_t *) num_in_edges->get_raw_arr(),
				(const fg::vsize_t *) num_out_edges->get_raw_arr(),
				header);
	vindex->dump(index_file);

	// Construct the file for the adjacency list file.
	FILE *f_graph = fopen(adj_file.c_str(), "w");
	if (f_graph == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format("open %1%: %2%")
			% adj_file % strerror(errno);
		return -1;
	}
	fwrite(&header, sizeof(header), 1, f_graph);
	bool ret = mem_vector::cast(in_adjs->cat())->export2(f_graph);
	assert(ret);
	ret = mem_vector::cast(out_adjs->cat())->export2(f_graph);
	assert(ret);
	fclose(f_graph);

	return 0;
}
