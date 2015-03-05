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

#include "generic_type.h"
#include "data_io.h"
#include "data_frame.h"
#include "mem_data_frame.h"
#include "mem_vector.h"
#include "vector_vector.h"

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

	df.get_vec(0)->append(froms);
	df.get_vec(1)->append(tos);
	return froms->get_length();
}

/*
 * This applies to a vector of values corresponding to the same key,
 * and generates an adjacency list.
 */
class adj_apply_operate: public gr_apply_operate<data_frame>
{
	fg::edge_type etype;
public:
	adj_apply_operate(fg::edge_type etype) {
		this->etype = etype;
	}

	void run(const void *key, const data_frame &val, mem_vector &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<fg::vertex_id_t>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<char>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void adj_apply_operate::run(const void *key, const data_frame &val,
		mem_vector &out) const
{
	fg::vertex_id_t vid = *(const fg::vertex_id_t *) key;
	if (vid == fg::INVALID_VERTEX_ID) {
		out.resize(0);
		return;
	}

	// Right now, we just assume there aren't attributes.
	size_t edge_data_size = 0;
	assert(val.get_num_vecs() == 2);

	assert(out.is_type<char>());
	type_mem_vector<fg::vertex_id_t>::ptr vec;
	if (etype == fg::edge_type::OUT_EDGE)
		vec = type_mem_vector<fg::vertex_id_t>::cast(val.get_vec("dest"));
	else
		vec = type_mem_vector<fg::vertex_id_t>::cast(val.get_vec("source"));
	assert(vec);

	vec->sort();
	// I added an invalid edge for each vertex.
	// The invalid edge is the maximal integer. After sorting, it's at the end
	// of the vector. So it's easy to remove it.
	fg::vsize_t num_edges = vec->get_length() - 1;
	size_t size = fg::ext_mem_undirected_vertex::num_edges2vsize(num_edges,
			edge_data_size);
	out.resize(size);

	// Even if we generate a directed, we still can use undirected vertex to
	// store one type of edges of a vertex.
	fg::in_mem_undirected_vertex<> v(vid, edge_data_size > 0);
	for (size_t i = 0; i < num_edges; i++)
		v.add_edge(fg::edge<>(vid, vec->get(i)));
	fg::ext_mem_undirected_vertex::serialize(v, out.get_raw_arr(), size, etype);
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "el2al edge_file\n");
		return -1;
	}

	std::string file_name = argv[1];
	std::vector<std::string> files;
	files.push_back(file_name);

	edge_parser parser;
	data_frame::ptr df = read_lines(files, parser);
	printf("There are %ld edges\n", df->get_num_entries());
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

	adj_apply_operate out_adj_op(fg::edge_type::OUT_EDGE);
	vector_vector::ptr out_adjs = df->groupby("source", out_adj_op);
	printf("There are %ld adjacency lists and they use %ld bytes in total\n",
			out_adjs->get_num_vecs(), out_adjs->get_tot_num_entries());

	adj_apply_operate in_adj_op(fg::edge_type::IN_EDGE);
	vector_vector::ptr in_adjs = df->groupby("dest", in_adj_op);
	printf("There are %ld adjacency lists and they use %ld bytes in total\n",
			in_adjs->get_num_vecs(), in_adjs->get_tot_num_entries());
}
