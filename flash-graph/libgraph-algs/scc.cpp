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
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <atomic>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "graph_engine.h"
#include "graph_config.h"
#include "FG_vector.h"
#include "FGlib.h"

namespace {

typedef std::unordered_map<vertex_id_t, vsize_t> comp_map_t;

class trim1_message: public vertex_message
{
	edge_type type;
public:
	trim1_message(edge_type type): vertex_message(
			sizeof(trim1_message), true) {
		this->type = type;
	}

	edge_type get_type() const {
		return type;
	}
};

class trim2_message: public vertex_message
{
	vertex_id_t comp_id;
public:
	trim2_message(vertex_id_t comp_id): vertex_message(
			sizeof(trim1_message), false) {
		this->comp_id = comp_id;
	}

	vertex_id_t get_comp_id() const {
		return comp_id;
	}
};

class fwbw_message: public vertex_message
{
	uint64_t color;
	vertex_id_t pivot;
	bool forward;
public:
	fwbw_message(vertex_id_t color, vertex_id_t pivot,
			bool forward): vertex_message(sizeof(fwbw_message), true) {
		this->color = color;
		this->pivot = pivot;
		this->forward = forward;
	}

	vertex_id_t get_pivot() const {
		return pivot;
	}

	uint64_t get_color() const {
		return color;
	}

	bool is_forward() const {
		return forward;
	}
};

class wcc_id_t
{
	vertex_id_t id;
public:
	wcc_id_t() {
		id = 0;
	}

	wcc_id_t(vertex_id_t id) {
		this->id = id;
	}

	bool operator<(const wcc_id_t &id) const {
		return this->id < id.id;
	}

	vertex_id_t get_id() const {
		return id;
	}
};

class wcc_comp_message: public vertex_message
{
	wcc_id_t id;
	uint64_t color;
public:
	wcc_comp_message(const wcc_id_t &id, uint64_t color): vertex_message(
			sizeof(wcc_comp_message), true) {
		this->id = id;
		this->color = color;
	}

	uint64_t get_color() const {
		return color;
	}

	const wcc_id_t &get_wcc_id() const {
		return id;
	}
};

enum scc_stage_t {
	// Initialize #edges.
	INIT,
	// Trim vertices with only in-edges or out-edges
	TRIM1,
	// Trim vertices in a SCC of size 2.
	TRIM2,
	// Additional trimming before each WCC.
	TRIM3,
	FWBW,
	// After the FWBW phase, we need to partition the remaining vertices.
	PARTITION,
	// Use label propagation with only in-edges.
	IN_WCC,
	// Use label propagation with only out-edges.
	OUT_WCC,
} scc_stage;

template<class T>
class bit_flags
{
	T v;
public:
	bit_flags() {
		v = 0;
	}

	void set_flag(int flag) {
		v = v | (0x1UL << flag);
	}

	void clear_flag(int flag) {
		v = v & (~(0x1UL << flag));
	}

	bool test_flag(int flag) const {
		return v & (0x1UL << flag);
	}
};

class fwbw_state
{
	enum {
		FW_COLOR,
		BW_COLOR,
		FW_BFS,
		BW_BFS,
		ASSIGNED,
		FW_VISITED,
		BW_VISITED,
		WCC_UPDATED,
	};

	static const int COLOR_OFF = 60;

	vertex_id_t base_color;
	vertex_id_t pivot;
	bit_flags<short> flags;
public:
	fwbw_state() {
		base_color = 0;
		pivot = INVALID_VERTEX_ID;
	}

	uint64_t get_color() const {
		return ((uint64_t) base_color)
			| ((((uint64_t) flags.test_flag(FW_COLOR)) << (1 + COLOR_OFF))
					| (((uint64_t) flags.test_flag(BW_COLOR)) << COLOR_OFF));
	}

	void set_pivot(vertex_id_t pivot) {
		this->pivot = pivot;
	}

	vertex_id_t get_pivot() const {
		return pivot;
	}

	vertex_id_t get_comp_id() const {
		assert(is_assigned());
		return pivot;
	}

	// Test if the vertex is assigned to a component.
	bool is_assigned() const {
		return flags.test_flag(ASSIGNED);
	}

	void assign_new_fw_color() {
		base_color = pivot;
		flags.clear_flag(BW_COLOR);
		flags.set_flag(FW_COLOR);
	}

	void assign_new_bw_color() {
		base_color = pivot;
		flags.clear_flag(FW_COLOR);
		flags.set_flag(BW_COLOR);
	}

	void assign_new_color(vertex_id_t new_color) {
		base_color = new_color;
		flags.clear_flag(FW_COLOR);
		flags.clear_flag(BW_COLOR);
	}

	void clear_flags() {
		flags.clear_flag(FW_BFS);
		flags.clear_flag(BW_BFS);
		flags.clear_flag(FW_VISITED);
		flags.clear_flag(BW_VISITED);
	}

	bool has_fw_visited() const {
		return flags.test_flag(FW_VISITED);
	}

	bool has_bw_visited() const {
		return flags.test_flag(BW_VISITED);
	}

	void set_fw_visited() {
		flags.set_flag(FW_VISITED);
	}

	void set_bw_visited() {
		flags.set_flag(BW_VISITED);
	}

	void set_fw() {
		return flags.set_flag(FW_BFS);
	}

	void set_bw() {
		return flags.set_flag(BW_BFS);
	}

	bool is_fw() const {
		return flags.test_flag(FW_BFS);
	}

	bool is_bw() const {
		return flags.test_flag(BW_BFS);
	}

	bool is_wcc_updated() const {
		return flags.test_flag(WCC_UPDATED);
	}

	void set_wcc_updated() {
		flags.set_flag(WCC_UPDATED);
	}

	void clear_wcc_updated() {
		flags.clear_flag(WCC_UPDATED);
	}
};

struct trim1_state
{
	// for trimming
	vsize_t num_in_edges;
	vsize_t num_out_edges;
};

struct wcc_state
{
	fwbw_state fwbw;
	wcc_id_t wcc_min;
};

class scc_vertex: public compute_directed_vertex
{
	vertex_id_t id;
	vsize_t comp_id;
	union scc_state {
		trim1_state trim1;
		fwbw_state fwbw;
		wcc_state wcc;

		scc_state() {
			memset(this, 0, sizeof(*this));
		}
	} state;
	vsize_t num_in_edges;
	vsize_t num_out_edges;

public:
	scc_vertex(vertex_id_t id): compute_directed_vertex(id) {
		this->id = id;
		comp_id = INVALID_VERTEX_ID;
		num_in_edges = 0;
		num_out_edges = 0;
	}

	vertex_id_t get_id() const {
		return id;
	}

	vsize_t get_degree() const {
		return num_in_edges + num_out_edges;
	}

	vsize_t get_num_in_edges() const {
		return num_in_edges;
	}

	vsize_t get_num_out_edges() const {
		return num_out_edges;
	}

	bool is_assigned() const {
		return comp_id != INVALID_VERTEX_ID;
	}

	vertex_id_t get_comp_id() const {
		return comp_id;
	}

	uint64_t get_color() const {
		return state.fwbw.get_color();
	}

	void init_trim1() {
		state.trim1.num_out_edges = get_num_out_edges();
		state.trim1.num_in_edges = get_num_in_edges();
	}

	void init_wcc() {
		state.wcc.wcc_min = wcc_id_t(get_id());
		state.fwbw.set_wcc_updated();
	}

	void reset_for_fwbw() {
		state.fwbw = fwbw_state();
	}

	void init_fwbw() {
		state.fwbw.set_fw();
		state.fwbw.set_bw();
		state.fwbw.set_pivot(get_id());
	}

	void post_wcc_init() {
		assert(!state.fwbw.has_fw_visited());
		assert(!state.fwbw.has_bw_visited());
		state.fwbw.assign_new_color(state.wcc.wcc_min.get_id());
	}

	void run(vertex_program &prog) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::INIT:
				run_stage_init(prog);
				break;
			case scc_stage_t::TRIM1:
				run_stage_trim1(prog);
				break;
			case scc_stage_t::TRIM2:
				run_stage_trim2(prog);
				break;
			case scc_stage_t::TRIM3:
				run_stage_trim3(prog);
				break;
			case scc_stage_t::FWBW:
				run_stage_FWBW(prog);
				break;
			case scc_stage_t::PARTITION:
				run_stage_part(prog);
				break;
			case scc_stage_t::IN_WCC:
			case scc_stage_t::OUT_WCC:
				run_stage_wcc(prog);
				break;
			default:
				ABORT_MSG("wrong SCC stage");
		}
	}

	void run_stage_init(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertex_headers(&id, 1);
	}

	void run_stage_trim1(vertex_program &prog);
	void run_stage_trim2(vertex_program &prog);
	void run_stage_trim3(vertex_program &prog);
	void run_stage_FWBW(vertex_program &prog);
	void run_stage_part(vertex_program &prog);
	void run_stage_wcc(vertex_program &prog);

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::TRIM1:
				run_stage_trim1(prog, vertex);
				break;
			case scc_stage_t::TRIM2:
				run_stage_trim2(prog, vertex);
				break;
			case scc_stage_t::TRIM3:
				run_stage_trim3(prog, vertex);
				break;
			case scc_stage_t::FWBW:
				run_stage_FWBW(prog, vertex);
				break;
			case scc_stage_t::PARTITION:
				run_stage_part(prog, vertex);
				break;
			case scc_stage_t::IN_WCC:
			case scc_stage_t::OUT_WCC:
				run_stage_wcc(prog, vertex);
				break;
			default:
				ABORT_MSG("wrong SCC stage");
		}
	}

	void run_stage_trim1(vertex_program &prog, const page_vertex &vertex);
	void run_stage_trim2(vertex_program &prog, const page_vertex &vertex);
	void run_stage_trim3(vertex_program &prog, const page_vertex &vertex);
	void run_stage_FWBW(vertex_program &prog, const page_vertex &vertex);
	void run_stage_part(vertex_program &prog, const page_vertex &vertex);
	void run_stage_wcc(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::TRIM1:
				run_on_message_stage_trim1(prog, msg);
				break;
			case scc_stage_t::TRIM2:
				run_on_message_stage_trim2(prog, msg);
				break;
			case scc_stage_t::TRIM3:
				run_on_message_stage_trim3(prog, msg);
				break;
			case scc_stage_t::FWBW:
				run_on_message_stage_FWBW(prog, msg);
				break;
			case scc_stage_t::PARTITION:
				run_on_message_stage_part(prog, msg);
				break;
			case scc_stage_t::IN_WCC:
			case scc_stage_t::OUT_WCC:
				run_on_message_stage_wcc(prog, msg);
				break;
			default:
				ABORT_MSG("wrong SCC stage");
		}
	}

	void run_on_message_stage_trim1(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_trim2(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_trim3(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_FWBW(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_part(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_wcc(vertex_program &prog, const vertex_message &msg);

	void run_on_vertex_header(vertex_program &prog, const vertex_header &header) {
		const directed_vertex_header &dheader = (const directed_vertex_header &) header;
		this->num_in_edges = dheader.get_num_in_edges();
		this->num_out_edges = dheader.get_num_out_edges();
	}

	vertex_id_t get_result() const {
		if (get_degree() > 0)
			return get_comp_id();
		else
			return INVALID_VERTEX_ID;
	}
};

class part_vertex_program: public vertex_program_impl<scc_vertex>
{
	std::vector<vertex_id_t> remain_vertices;
	comp_map_t comp_sizes;
	size_t num_assigned;
public:
	typedef std::shared_ptr<part_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<part_vertex_program, vertex_program>(prog);
	}

	part_vertex_program() {
		num_assigned = 0;
	}

	void add_remain_vertex(vertex_id_t id) {
		remain_vertices.push_back(id);
	}

	void assign_vertex(vertex_id_t comp_id) {
		num_assigned++;
		comp_map_t::iterator it = comp_sizes.find(comp_id);
		if (it == comp_sizes.end())
			comp_sizes.insert(comp_map_t::value_type(comp_id, 1));
		else
			it->second++;
	}

	size_t get_num_assigned() const {
		return num_assigned;
	}

	const std::vector<vertex_id_t> &get_remain_vertices() const {
		return remain_vertices;
	}

	void merge_comp_size(comp_map_t &merged_comp_sizes) const {
		for (comp_map_t::const_iterator it = comp_sizes.begin();
				it != comp_sizes.end(); it++) {
			vertex_id_t comp_id = it->first;
			size_t size = it->second;
			comp_map_t::iterator it1 = merged_comp_sizes.find(comp_id);
			if (it1 != merged_comp_sizes.end())
				it1->second += size;
			else
				merged_comp_sizes.insert(comp_map_t::value_type(comp_id, size));
		}
	}
};

class part_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new part_vertex_program());
	}
};

class trim_vertex_program: public vertex_program_impl<scc_vertex>
{
	size_t num_trims;
public:
	typedef std::shared_ptr<trim_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<trim_vertex_program, vertex_program>(prog);
	}

	trim_vertex_program() {
		num_trims = 0;
	}

	void trim_vertex(int num) {
		num_trims += num;
	}

	size_t get_num_trimmed() const {
		return num_trims;
	}
};

class trim_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new trim_vertex_program());
	}
};

void scc_vertex::run_stage_trim1(vertex_program &prog)
{
	if (state.trim1.num_in_edges == 0 || state.trim1.num_out_edges == 0) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);

		// This vertex has to be a SCC itself.
		comp_id = id;
		((trim_vertex_program &) prog).trim_vertex(1);
	}
}

void scc_vertex::run_stage_trim1(vertex_program &prog, const page_vertex &vertex)
{
	// The vertices on the other side of the edges should reduce their degree
	// by 1. They have the opposite direction of the edges.
	edge_type type = edge_type::NONE;
	if (vertex.get_num_edges(edge_type::IN_EDGE) > 0) {
		assert(vertex.get_num_edges(edge_type::OUT_EDGE) == 0);
		type = edge_type::OUT_EDGE;
	}
	else if (vertex.get_num_edges(edge_type::OUT_EDGE) > 0) {
		assert(vertex.get_num_edges(edge_type::IN_EDGE) == 0);
		type = edge_type::IN_EDGE;
	}
	if (type != edge_type::NONE) {
		trim1_message msg(type);
		int num_edges = vertex.get_num_edges(BOTH_EDGES);
		edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0,
				num_edges);
		prog.multicast_msg(it, msg);
	}
}

void scc_vertex::run_on_message_stage_trim1(vertex_program &prog,
		const vertex_message &msg1)
{
	const trim1_message &msg = (const trim1_message &) msg1;
	switch(msg.get_type()) {
		case edge_type::IN_EDGE:
			assert(state.trim1.num_in_edges > 0);
			state.trim1.num_in_edges--;
			break;
		case edge_type::OUT_EDGE:
			assert(state.trim1.num_out_edges > 0);
			state.trim1.num_out_edges--;
			break;
		default:
			ABORT_MSG("wrong message type");
	}
}

void scc_vertex::run_stage_trim2(vertex_program &prog)
{
	vertex_id_t id = get_id();
	if (get_num_in_edges() == 1) {
		// TODO requesting partial vertices causes errors.
		request_vertices(&id, 1);
	}
	else if (get_num_out_edges() == 1) {
		// TODO requesting partial vertices causes errors.
		request_vertices(&id, 1);
	}
}

static inline bool contain_edge(const page_vertex &vertex, edge_type type,
		vertex_id_t id)
{
	return std::binary_search(vertex.get_neigh_begin(type),
			vertex.get_neigh_end(type), id);
}

void scc_vertex::run_stage_trim2(vertex_program &prog, const page_vertex &vertex)
{
	assert(vertex.get_id() == get_id());
	// Ideally, we should use the remaining in-edges or out-edges,
	// but we don't know which edges have been removed, so we just
	// use the original number of edges.
	if (get_num_in_edges() == 1) {
		page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(edge_type::IN_EDGE);
		vertex_id_t neighbor = *it;
		// If the only in-edge is to itself, it's a SCC itself.
		if (neighbor == get_id()) {
			comp_id = get_id();
			((trim_vertex_program &) prog).trim_vertex(1);
		}
		else {
			scc_vertex &neigh_v = (scc_vertex &) prog.get_graph().get_vertex(neighbor);
			// If the vertex's out-edge list contains the neighbor,
			// it means the neighbor's only in-edge connect to this vertex.
			if (get_id() < neighbor
					&& neigh_v.get_num_in_edges() == 1
					&& contain_edge(vertex, edge_type::OUT_EDGE, neighbor)) {
				comp_id = get_id();
				trim2_message msg(get_id());
				prog.send_msg(neighbor, msg);
				((trim_vertex_program &) prog).trim_vertex(2);
			}
		}
	}
	else if (get_num_out_edges() == 1) {
		page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(edge_type::OUT_EDGE);
		vertex_id_t neighbor = *it;
		// If the only in-edge is to itself, it's a SCC itself.
		if (neighbor == get_id()) {
			comp_id = get_id();
			((trim_vertex_program &) prog).trim_vertex(1);
		}
		else {
			scc_vertex &neigh_v = (scc_vertex &) prog.get_graph().get_vertex(neighbor);
			// The same as above.
			if (get_id() < neighbor
					&& neigh_v.get_num_out_edges() == 1
					&& contain_edge(vertex, edge_type::IN_EDGE, neighbor)) {
				comp_id = get_id();
				trim2_message msg(get_id());
				prog.send_msg(neighbor, msg);
				((trim_vertex_program &) prog).trim_vertex(2);
			}
		}
	}
	else
		assert(0);
}

void scc_vertex::run_on_message_stage_trim2(vertex_program &prog,
		const vertex_message &msg1)
{
	const trim2_message &msg = (const trim2_message &) msg1;
	comp_id = msg.get_comp_id();
}

void scc_vertex::run_stage_trim3(vertex_program &prog)
{
	vertex_id_t id = get_id();
	request_vertices(&id, 1);
}

std::atomic_long trim3_vertices;

void scc_vertex::run_stage_trim3(vertex_program &prog, const page_vertex &vertex)
{
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(IN_EDGE);
	stack_array<vertex_id_t, 1024> in_neighs(vertex.get_num_edges(IN_EDGE));
	int num_in_neighs = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(IN_EDGE); it != end_it; ++it) {
		vertex_id_t id = *it;
		scc_vertex &neigh = (scc_vertex &) prog.get_graph().get_vertex(id);
		// We should ignore the neighbors that has been assigned to a component.
		// or has a different color.
		if (neigh.is_assigned()
				|| neigh.state.fwbw.get_color() != state.fwbw.get_color())
			continue;

		in_neighs[num_in_neighs++] = id;
	}

	end_it = vertex.get_neigh_end(OUT_EDGE);
	stack_array<vertex_id_t, 1024> out_neighs(vertex.get_num_edges(OUT_EDGE));
	int num_out_neighs = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(OUT_EDGE); it != end_it; ++it) {
		vertex_id_t id = *it;
		scc_vertex &neigh = (scc_vertex &) prog.get_graph().get_vertex(id);
		// We should ignore the neighbors that has been assigned to a component.
		// or has a different color.
		if (neigh.is_assigned()
				|| neigh.state.fwbw.get_color() != state.fwbw.get_color())
			continue;

		out_neighs[num_out_neighs++] = id;
	}

	if (num_in_neighs == 0 || num_out_neighs == 0) {
		trim3_vertices++;
		// This vertex has been isolated, it can assign to a SCC now.
		comp_id = get_id();
		if (num_in_neighs > 0)
			prog.activate_vertices(in_neighs.data(), num_in_neighs);
		if (num_out_neighs > 0)
			prog.activate_vertices(out_neighs.data(), num_out_neighs);
	}
}

void scc_vertex::run_on_message_stage_trim3(vertex_program &prog,
		const vertex_message &msg)
{
}

void scc_vertex::run_stage_FWBW(vertex_program &prog)
{
	// If the vertex has been visited in both directions,
	// we don't need to do anything.
	if (state.fwbw.has_fw_visited() && state.fwbw.has_bw_visited())
		return;
	// If the vertex has been visisted in forward direction, and it doesn't
	// need to visit other in the backwoard direction, then we don't need to
	// do anything.
	if (state.fwbw.has_fw_visited() && !state.fwbw.is_bw())
		return;
	// The same for the other direction.
	if (state.fwbw.has_bw_visited() && !state.fwbw.is_fw())
		return;

	// It's possible that the vertex is activated by another vertex of
	// a different color. If that is the case, the vertex may not have
	// the forward BFS flag nor the backward BFS flag. Do nothing.
	if (!state.fwbw.is_bw() && !state.fwbw.is_fw())
		return;

	vertex_id_t id = get_id();
	request_vertices(&id, 1);
}

void scc_vertex::run_stage_FWBW(vertex_program &prog, const page_vertex &vertex)
{
	bool do_some = false;

	if (state.fwbw.is_bw()) {
		do_some = true;
		state.fwbw.set_bw_visited();
		fwbw_message msg(state.fwbw.get_color(), state.fwbw.get_pivot(), false);
		int num_edges = vertex.get_num_edges(IN_EDGE);
		edge_seq_iterator it = vertex.get_neigh_seq_it(IN_EDGE, 0,
				num_edges);
		prog.multicast_msg(it, msg);
	}

	if (state.fwbw.is_fw()) {
		do_some = true;
		state.fwbw.set_fw_visited();
		fwbw_message msg(state.fwbw.get_color(), state.fwbw.get_pivot(), true);
		int num_edges = vertex.get_num_edges(OUT_EDGE);
		edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0,
				num_edges);
		prog.multicast_msg(it, msg);
	}
	assert(do_some);
}

void scc_vertex::run_on_message_stage_FWBW(vertex_program &prog,
		const vertex_message &msg1)
{
	uint64_t color = state.fwbw.get_color();
	const fwbw_message &msg = (const fwbw_message &) msg1;
	// If the current vertex has a different color, it means it's in
	// a different partition. The vertex can just ignore the message.
	if (msg.get_color() != color)
		return;

	state.fwbw.set_pivot(msg.get_pivot());
	if (msg.is_forward())
		state.fwbw.set_fw();
	else
		state.fwbw.set_bw();
}

void scc_vertex::run_stage_part(vertex_program &prog)
{
	if (state.fwbw.is_fw() && state.fwbw.is_bw()) {
		comp_id = state.fwbw.get_pivot();
		((part_vertex_program &) prog).assign_vertex(comp_id);
	}
	else if (state.fwbw.is_fw())
		state.fwbw.assign_new_fw_color();
	else if (state.fwbw.is_bw())
		state.fwbw.assign_new_bw_color();
	state.fwbw.clear_flags();

	if (!is_assigned())
		((part_vertex_program &) prog).add_remain_vertex(get_id());
}

void scc_vertex::run_stage_part(vertex_program &prog, const page_vertex &vertex)
{
}

void scc_vertex::run_on_message_stage_part(vertex_program &prog,
		const vertex_message &msg)
{
}

void scc_vertex::run_stage_wcc(vertex_program &prog)
{
	if (state.fwbw.is_wcc_updated()) {
		state.fwbw.clear_wcc_updated();
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}
}

void scc_vertex::run_stage_wcc(vertex_program &prog, const page_vertex &vertex)
{
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	edge_type type;
	if (scc_stage == IN_WCC)
		type = IN_EDGE;
	else
		type = OUT_EDGE;
	wcc_comp_message msg(state.wcc.wcc_min, state.fwbw.get_color());
	int num_edges = vertex.get_num_edges(type);
	edge_seq_iterator it = vertex.get_neigh_seq_it(type, 0, num_edges);
	prog.multicast_msg(it, msg);
}

void scc_vertex::run_on_message_stage_wcc(vertex_program &prog,
		const vertex_message &msg1)
{
	wcc_comp_message &msg = (wcc_comp_message &) msg1;
	// If the current vertex has a different color, it means it's in
	// a different partition. The vertex can just ignore the message.
	if (msg.get_color() != state.fwbw.get_color())
		return;

	if (msg.get_wcc_id() < state.wcc.wcc_min) {
		state.wcc.wcc_min = msg.get_wcc_id();
		state.fwbw.set_wcc_updated();
	}
}

#if 0
class sec_fwbw_filter: public vertex_filter
{
public:
	virtual bool keep(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		bool activate = !scc_v.is_assigned()
			&& scc_v.wcc_min.get_id() == scc_v.get_id();
		// If the vertex hasn't been assigned to a component,
		// let's use the result of wcc (which is stored in pivot) as the color
		if (!scc_v.is_assigned())
			scc_v.fwbw_state.assign_new_color(scc_v.wcc_min.get_id());
		if (activate)
			scc_v.init_fwbw();
		return activate;
	}
};
#endif

class trim1_initializer: public vertex_initializer
{
public:
	void init(compute_vertex &v) {
		scc_vertex &sv = (scc_vertex &) v;
		sv.init_trim1();
	}
};

/**
 * This initializes the start vertices for forward-backward BFS.
 */
class fwbw_initializer: public vertex_initializer
{
public:
	void init(compute_vertex &v) {
		scc_vertex &sv = (scc_vertex &) v;
		assert(!sv.is_assigned());
		sv.init_fwbw();
	}
};

/**
 * This prepares all vertices in the graph for forward-backward BFS.
 */
class fwbw_reset: public vertex_initializer
{
public:
	void init(compute_vertex &v) {
		scc_vertex &sv = (scc_vertex &) v;
		sv.reset_for_fwbw();
	}
};

class in_wcc_initializer: public vertex_initializer
{
public:
	virtual void init(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		if (scc_v.is_assigned())
			return;
		scc_v.init_wcc();
	}
};

class out_wcc_initializer: public vertex_initializer
{
public:
	virtual void init(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		if (scc_v.is_assigned())
			return;
		// OUT_WCC runs after IN_WCC, so we need to do post-WCC
		// initialization.
		scc_v.post_wcc_init();
		assert(scc_v.get_color() < INVALID_VERTEX_ID);
		scc_v.init_wcc();
	}
};

class max_degree_query: public vertex_query
{
	vsize_t max_degree;
	vertex_id_t max_id;
public:
	max_degree_query() {
		max_degree = 0;
		max_id = INVALID_VERTEX_ID;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		vsize_t degree = scc_v.get_degree();
		if (degree > max_degree && !scc_v.is_assigned()) {
			max_degree = degree;
			max_id = scc_v.get_id();
		}
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		max_degree_query *mdq = (max_degree_query *) q.get();
		if (max_degree < mdq->max_degree) {
			max_degree = mdq->max_degree;
			max_id = mdq->max_id;
		}
	}

	virtual ptr clone() {
		return vertex_query::ptr(new max_degree_query());
	}

	vertex_id_t get_max_id() const {
		return max_id;
	}
};

class post_wcc_query: public vertex_query
{
	// The largest-degree vertices in each color
	typedef std::unordered_map<uint64_t, vertex_id_t> color_map_t;
	color_map_t max_ids;
public:
	virtual void run(graph_engine &graph, compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		// Ignore the assigned vertex
		if (scc_v.is_assigned())
			return;
		scc_v.post_wcc_init();
		assert(scc_v.get_color() < INVALID_VERTEX_ID);

		color_map_t::iterator it = max_ids.find(scc_v.get_color());
		// The color doesn't exist;
		if (it == max_ids.end()) {
			max_ids.insert(color_map_t::value_type(scc_v.get_color(),
						scc_v.get_id()));
		}
		else {
			vertex_id_t curr_max_id = it->second;
			if (scc_v.get_degree()
					> ((scc_vertex &) graph.get_vertex(curr_max_id)).get_degree())
				it->second = scc_v.get_id();
		}
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		post_wcc_query *pwq = (post_wcc_query *) q.get();
		for (color_map_t::const_iterator it = pwq->max_ids.begin();
				it != pwq->max_ids.end(); it++) {
			uint64_t color = it->first;
			vertex_id_t id = it->second;
			scc_vertex &scc_v = (scc_vertex &) graph.get_vertex(id);
			assert(!scc_v.is_assigned());
			color_map_t::iterator it1 = this->max_ids.find(color);
			// The same color exists.
			if (it1 != this->max_ids.end()) {
				// If the vertex of the same color in the other query is larger
				if (scc_v.get_degree()
						> ((scc_vertex &) graph.get_vertex(it1->second)).get_degree())
					it1->second = id;
			}
			else
				this->max_ids.insert(color_map_t::value_type(color, it->second));
		}
	}

	virtual ptr clone() {
		return vertex_query::ptr(new post_wcc_query());
	}

	size_t get_max_ids(std::vector<vertex_id_t> &ids) const {
		for (color_map_t::const_iterator it = max_ids.begin();
				it != max_ids.end(); it++)
			ids.push_back(it->second);
		return max_ids.size();
	}
};

class comp_size_compare
{
public:
	bool operator()(const std::pair<vsize_t, int> &p1,
			const std::pair<vsize_t, int> &p2) {
		return p1.first < p2.first;
	}
};

}

#include "save_result.h"
FG_vector<vertex_id_t>::ptr compute_scc(FG_graph::ptr fg)
{
	graph_index::ptr index = NUMA_graph_index<scc_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	BOOST_LOG_TRIVIAL(info) << "SCC starts";
	BOOST_LOG_TRIVIAL(info) << "prof_file: " << graph_conf.get_prof_file();
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	std::vector<vertex_id_t> active_vertices;
	vertex_id_t max_v = 0;
	size_t num_comp1 = 0;

	struct timeval start, end, scc_start;
	scc_stage = scc_stage_t::INIT;
	gettimeofday(&start, NULL);
	scc_start = start;
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("init takes %1% seconds.") % time_diff(start, end);

	scc_stage = scc_stage_t::TRIM1;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(new trim1_initializer()),
			vertex_program_creater::ptr(new trim_vertex_program_creater()));
	graph->wait4complete();
	std::vector<vertex_program::ptr> trim_vprogs;
	graph->get_vertex_programs(trim_vprogs);
	BOOST_FOREACH(vertex_program::ptr vprog, trim_vprogs) {
		trim_vertex_program::ptr trim_vprog = trim_vertex_program::cast2(vprog);
		num_comp1 += trim_vprog->get_num_trimmed();
	}
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("trim1 takes %1% seconds. It trims %2% vertices")
		% time_diff(start, end) % num_comp1;

#if 0
	scc_stage = scc_stage_t::TRIM2;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(),
			vertex_program_creater::ptr(new trim_vertex_program_creater()));
	graph->wait4complete();
	graph->get_vertex_programs(trim_vprogs);
	BOOST_FOREACH(vertex_program::ptr vprog, trim_vprogs) {
		trim_vertex_program::ptr trim_vprog = trim_vertex_program::cast2(vprog);
		num_comp2 += trim_vprog->get_num_trimmed();
	}
	gettimeofday(&end, NULL);
	printf("trim2 takes %f seconds. It trims %ld vertices\n",
			time_diff(start, end), num_comp2);
#endif

	vertex_query::ptr mdq(new max_degree_query());
	graph->query_on_all(mdq);
	max_v = ((max_degree_query *) mdq.get())->get_max_id();
	scc_stage = scc_stage_t::FWBW;
	gettimeofday(&start, NULL);
	graph->init_all_vertices(vertex_initializer::ptr(new fwbw_reset()));
	scc_vertex &v = (scc_vertex &) index->get_vertex(max_v);
	v.init_fwbw();
	graph->start(&max_v, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("FWBW takes %1% seconds") % time_diff(start, end);

	scc_stage = scc_stage_t::PARTITION;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initializer::ptr(),
			vertex_program_creater::ptr(new part_vertex_program_creater()));
	graph->wait4complete();

	std::vector<vertex_program::ptr> part_vprogs;
	graph->get_vertex_programs(part_vprogs);
	size_t largest_comp_size = 0;
	BOOST_FOREACH(vertex_program::ptr vprog, part_vprogs) {
		part_vertex_program::ptr part_vprog = part_vertex_program::cast2(vprog);
		largest_comp_size += part_vprog->get_num_assigned();
		active_vertices.insert(active_vertices.begin(),
				part_vprog->get_remain_vertices().begin(),
				part_vprog->get_remain_vertices().end());
	}
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("partition takes %1% seconds. Assign %2% vertices to components.")
		% time_diff(start, end) % largest_comp_size;
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("after partition, finding %1% active vertices takes %2% seconds.")
		% active_vertices.size() % time_diff(start, end);

	do {
		scc_stage = scc_stage_t::TRIM3;
		trim3_vertices = 0;
		graph->start(active_vertices.data(), active_vertices.size());
		graph->wait4complete();
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("trim3 %1% vertices") % trim3_vertices.load();
		num_comp1 += trim3_vertices.load();

		scc_stage = scc_stage_t::IN_WCC;
		graph->start(active_vertices.data(), active_vertices.size(),
				std::shared_ptr<vertex_initializer>(new in_wcc_initializer()));
		graph->wait4complete();

		scc_stage = scc_stage_t::OUT_WCC;
		graph->start(active_vertices.data(), active_vertices.size(),
				std::shared_ptr<vertex_initializer>(new out_wcc_initializer()));
		graph->wait4complete();
		vertex_query::ptr mdq1(new post_wcc_query());
		graph->query_on_all(mdq1);

		std::vector<vertex_id_t> fwbw_starts;
		((post_wcc_query *) mdq1.get())->get_max_ids(fwbw_starts);
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("FWBW starts on %1% vertices") % fwbw_starts.size();
		scc_stage = scc_stage_t::FWBW;
		graph->start(fwbw_starts.data(), fwbw_starts.size(),
				vertex_initializer::ptr(new fwbw_initializer()));
		graph->wait4complete();

		scc_stage = scc_stage_t::PARTITION;
		graph->start(active_vertices.data(), active_vertices.size(),
				vertex_initializer::ptr(),
				vertex_program_creater::ptr(new part_vertex_program_creater()));
		graph->wait4complete();

		std::vector<vertex_program::ptr> part_vprogs;
		graph->get_vertex_programs(part_vprogs);
		active_vertices.clear();
		size_t fwbw_vertices = 0;
		BOOST_FOREACH(vertex_program::ptr vprog, part_vprogs) {
			part_vertex_program::ptr part_vprog = part_vertex_program::cast2(vprog);
			// Count the number of vertices assigned to a component by FWBW.
			fwbw_vertices += part_vprog->get_num_assigned();
			// We figure out here the vertices that haven't been assigned to
			// a component.
			active_vertices.insert(active_vertices.begin(),
					part_vprog->get_remain_vertices().begin(),
					part_vprog->get_remain_vertices().end());
		}
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("partitioning assigns %1% vertices to components.")
			% fwbw_vertices;
		BOOST_LOG_TRIVIAL(info)
			<< boost::format("There are %1% vertices left unassigned")
			% active_vertices.size();
	} while (!active_vertices.empty());
	BOOST_LOG_TRIVIAL(info)
			<< boost::format("scc takes %1% seconds") % time_diff(scc_start, end);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif

	FG_vector<vertex_id_t>::ptr vec = FG_vector<vertex_id_t>::create(graph);
	graph->query_on_all(vertex_query::ptr(
				new save_query<vertex_id_t, scc_vertex>(vec)));
	return vec;
}
