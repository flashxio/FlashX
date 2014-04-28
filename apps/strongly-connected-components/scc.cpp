/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <signal.h>
#include <google/profiler.h>

#include <atomic>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

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
	vsize_t deg;
	vertex_id_t id;
public:
	wcc_id_t() {
		deg = 0;
		id = 0;
	}

	wcc_id_t(vsize_t deg, vertex_id_t id) {
		this->deg = deg;
		this->id = id;
	}

	bool operator>(const wcc_id_t &id) const {
		if (this->deg == id.deg)
			return this->id > id.id;
		else
			return this->deg > id.deg;
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
	// Trim vertices with only in-edges or out-edges
	TRIM1,
	// Trim vertices in a SCC of size 2.
	TRIM2,
	// Additional trimming before each WCC.
	TRIM3,
	FWBW,
	// After the FWBW phase, we need to partition the remaining vertices.
	PARTITION,
	WCC,
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
	wcc_id_t wcc_max;
};

class scc_vertex: public compute_directed_vertex
{
	vsize_t comp_id;
	union scc_state {
		trim1_state trim1;
		fwbw_state fwbw;
		wcc_state wcc;

		scc_state() {
			memset(this, 0, sizeof(*this));
		}
	} state;

public:
	scc_vertex() {
		comp_id = INVALID_VERTEX_ID;
	}

	scc_vertex(vertex_id_t id, const vertex_index *index1): compute_directed_vertex(
			id, index1) {
		comp_id = INVALID_VERTEX_ID;
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
		state.wcc.wcc_max = wcc_id_t(get_num_out_edges() + get_num_in_edges(), get_id());
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
		state.fwbw.assign_new_color(state.wcc.wcc_max.get_id());
	}

	void run(vertex_program &prog) {
		if (is_assigned())
			return;

		switch(scc_stage) {
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
			case scc_stage_t::WCC:
				run_stage_wcc(prog);
				break;
			default:
				assert(0);
		}
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
			case scc_stage_t::WCC:
				run_stage_wcc(prog, vertex);
				break;
			default:
				assert(0);
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
			case scc_stage_t::WCC:
				run_on_message_stage_wcc(prog, msg);
				break;
			default:
				assert(0);
		}
	}

	void run_on_message_stage_trim1(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_trim2(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_trim3(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_FWBW(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_part(vertex_program &prog, const vertex_message &msg);
	void run_on_message_stage_wcc(vertex_program &prog, const vertex_message &msg);
};

class part_vertex_program: public vertex_program_impl<scc_vertex>
{
	std::vector<vertex_id_t> remain_vertices;
public:
	typedef std::shared_ptr<part_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<part_vertex_program, vertex_program>(prog);
	}

	void add_remain_vertex(vertex_id_t id) {
		remain_vertices.push_back(id);
	}

	const std::vector<vertex_id_t> &get_remain_vertices() const {
		return remain_vertices;
	}
};

class part_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new part_vertex_program());
	}
};

std::atomic_ulong trim1_vertices;

void scc_vertex::run_stage_trim1(vertex_program &prog)
{
	if (state.trim1.num_in_edges == 0 || state.trim1.num_out_edges == 0) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);

		// This vertex has to be a SCC itself.
		comp_id = id;
		trim1_vertices++;
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
			assert(0);
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

std::atomic_ulong trim2_vertices;

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
			trim2_vertices++;
		}
		else {
			scc_vertex &neigh_v = (scc_vertex &) prog.get_graph().get_vertex(neighbor);
			// If the vertex's out-edge list contains the neighbor,
			// it means the neighbor's only in-edge connect to this vertex.
			if (get_id() < neighbor
					&& neigh_v.get_num_in_edges() == 1
					&& vertex.contain_edge(edge_type::OUT_EDGE, neighbor)) {
				comp_id = get_id();
				trim2_message msg(get_id());
				prog.send_msg(neighbor, msg);
				trim2_vertices += 2;
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
			trim2_vertices++;
		}
		else {
			scc_vertex &neigh_v = (scc_vertex &) prog.get_graph().get_vertex(neighbor);
			// The same as above.
			if (get_id() < neighbor
					&& neigh_v.get_num_out_edges() == 1
					&& vertex.contain_edge(edge_type::IN_EDGE, neighbor)) {
				comp_id = get_id();
				trim2_message msg(get_id());
				prog.send_msg(neighbor, msg);
				trim2_vertices += 2;
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

std::atomic_ulong fwbw_vertices;

void scc_vertex::run_stage_part(vertex_program &prog)
{
	if (state.fwbw.is_fw() && state.fwbw.is_bw()) {
		comp_id = state.fwbw.get_pivot();
		fwbw_vertices++;
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
	wcc_comp_message msg(state.wcc.wcc_max, state.fwbw.get_color());
	int num_edges = vertex.get_num_edges(BOTH_EDGES);
	edge_seq_iterator it = vertex.get_neigh_seq_it(BOTH_EDGES, 0,
			num_edges);
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

	if (msg.get_wcc_id() > state.wcc.wcc_max) {
		state.wcc.wcc_max = msg.get_wcc_id();
		state.fwbw.set_wcc_updated();
	}
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

#if 0
class sec_fwbw_filter: public vertex_filter
{
public:
	virtual bool keep(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		bool activate = !scc_v.is_assigned()
			&& scc_v.wcc_max.get_id() == scc_v.get_id();
		// If the vertex hasn't been assigned to a component,
		// let's use the result of wcc (which is stored in pivot) as the color
		if (!scc_v.is_assigned())
			scc_v.fwbw_state.assign_new_color(scc_v.wcc_max.get_id());
		if (activate)
			scc_v.init_fwbw();
		return activate;
	}
};
#endif

class trim1_initiator: public vertex_initiator
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
class fwbw_initiator: public vertex_initiator
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
class fwbw_reset: public vertex_initiator
{
public:
	void init(compute_vertex &v) {
		scc_vertex &sv = (scc_vertex &) v;
		sv.reset_for_fwbw();
	}
};

class wcc_initiator: public vertex_initiator
{
public:
	virtual void init(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		if (scc_v.is_assigned())
			return;
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
		if (graph.get_vertex_edges(v.get_id()) > max_degree
				&& !scc_v.is_assigned()) {
			max_degree = graph.get_vertex_edges(v.get_id());
			max_id = v.get_id();
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
	std::unordered_set<vertex_id_t> start_vertices;
public:
	virtual void run(graph_engine &graph, compute_vertex &v) {
		scc_vertex &sv = (scc_vertex &) v;
		// Ignore the assigned vertex
		if (sv.is_assigned())
			return;
		sv.post_wcc_init();
		assert(sv.get_color() < INVALID_VERTEX_ID);
		start_vertices.insert(sv.get_color());
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		post_wcc_query *pwq = (post_wcc_query *) q.get();
		start_vertices.insert(pwq->start_vertices.begin(), pwq->start_vertices.end());
	}

	virtual ptr clone() {
		return vertex_query::ptr(new post_wcc_query());
	}

	size_t get_max_ids(std::vector<vertex_id_t> &ids) const {
		ids.insert(ids.begin(), start_vertices.begin(), start_vertices.end());
		return start_vertices.size();
	}
};

void print_usage()
{
	fprintf(stderr,
			"scc [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-s size: the output min component size\n");
	fprintf(stderr, "-o file: output the component size to the file\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	struct timeval start, end, scc_start;
	int opt;
	std::string confs;
	std::string output_file;
	size_t min_comp_size = 0;
	int num_opts = 0;

	while ((opt = getopt(argc, argv, "c:s:o:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 's':
				min_comp_size = atoi(optarg);
				num_opts++;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 3) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];

	config_map configs(conf_file);
	configs.add_options(confs);
	graph_conf.init(configs);
	graph_conf.print();

	signal(SIGINT, int_handler);
	init_io_system(configs);

	graph_index *index = NUMA_graph_index<scc_vertex>::create(index_file,
			graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(graph_conf.get_num_threads(),
			params.get_num_nodes(), graph_file, index);
	printf("SCC starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	scc_stage = scc_stage_t::TRIM1;
	gettimeofday(&start, NULL);
	scc_start = start;
	graph->start_all(vertex_initiator::ptr(new trim1_initiator()));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("trim1 takes %f seconds. It trims %ld vertices\n",
			time_diff(start, end), trim1_vertices.load());

	scc_stage = scc_stage_t::TRIM2;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("trim2 takes %f seconds. It trims %ld vertices\n",
			time_diff(start, end), trim2_vertices.load());

	vertex_query::ptr mdq(new max_degree_query());
	graph->query_on_all(mdq);
	vertex_id_t max_v = ((max_degree_query *) mdq.get())->get_max_id();
	scc_stage = scc_stage_t::FWBW;
	gettimeofday(&start, NULL);
	graph->init_all_vertices(vertex_initiator::ptr(new fwbw_reset()));
	scc_vertex &v = (scc_vertex &) index->get_vertex(max_v);
	v.init_fwbw();
	graph->start(&max_v, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("FWBW takes %f seconds\n", time_diff(start, end));

	scc_stage = scc_stage_t::PARTITION;
	fwbw_vertices = 0;
	gettimeofday(&start, NULL);
	graph->start_all(vertex_initiator::ptr(),
			vertex_program_creater::ptr(new part_vertex_program_creater()));
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("partition takes %f seconds. Assign %ld vertices to components.\n",
			time_diff(start, end), fwbw_vertices.load());

	gettimeofday(&start, NULL);
	std::vector<vertex_program::ptr> part_vprogs;
	graph->get_vertex_programs(part_vprogs);
	std::vector<vertex_id_t> active_vertices;
	BOOST_FOREACH(vertex_program::ptr vprog, part_vprogs) {
		part_vertex_program::ptr part_vprog = part_vertex_program::cast2(vprog);
		active_vertices.insert(active_vertices.begin(),
				part_vprog->get_remain_vertices().begin(),
				part_vprog->get_remain_vertices().end());
	}
	gettimeofday(&end, NULL);
	printf("after partition, finding %ld active vertices takes %f seconds.\n",
			active_vertices.size(), time_diff(start, end));

	do {
		scc_stage = scc_stage_t::TRIM3;
		trim3_vertices = 0;
		gettimeofday(&start, NULL);
		graph->start(active_vertices.data(), active_vertices.size());
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("trim3 takes %f seconds, and trime %ld vertices\n",
				time_diff(start, end),
				trim3_vertices.load(std::memory_order_relaxed));

		scc_stage = scc_stage_t::WCC;
		gettimeofday(&start, NULL);
		graph->start(active_vertices.data(), active_vertices.size(),
				std::shared_ptr<vertex_initiator>(new wcc_initiator()));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("WCC takes %f seconds\n", time_diff(start, end));

		gettimeofday(&start, NULL);
		vertex_query::ptr mdq1(new post_wcc_query());
		graph->query_on_all(mdq1);
		std::vector<vertex_id_t> fwbw_starts;
		((post_wcc_query *) mdq1.get())->get_max_ids(fwbw_starts);
		printf("FWBW starts on %ld vertices\n", fwbw_starts.size());
		scc_stage = scc_stage_t::FWBW;
		graph->start(fwbw_starts.data(), fwbw_starts.size(),
				vertex_initiator::ptr(new fwbw_initiator()));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("FWBW takes %f seconds\n", time_diff(start, end));

		scc_stage = scc_stage_t::PARTITION;
		fwbw_vertices = 0;
		gettimeofday(&start, NULL);
		graph->start(active_vertices.data(), active_vertices.size(),
				vertex_initiator::ptr(),
				vertex_program_creater::ptr(new part_vertex_program_creater()));
		graph->wait4complete();

		std::vector<vertex_program::ptr> part_vprogs;
		graph->get_vertex_programs(part_vprogs);
		active_vertices.clear();
		BOOST_FOREACH(vertex_program::ptr vprog, part_vprogs) {
			part_vertex_program::ptr part_vprog = part_vertex_program::cast2(vprog);
			active_vertices.insert(active_vertices.begin(),
					part_vprog->get_remain_vertices().begin(),
					part_vprog->get_remain_vertices().end());
		}
		gettimeofday(&end, NULL);
		printf("partition takes %f seconds. Assign %ld vertices to components.\n",
				time_diff(start, end), fwbw_vertices.load());
		printf("There are %ld vertices left unassigned\n", active_vertices.size());
	} while (!active_vertices.empty());

	NUMA_graph_index<scc_vertex>::const_iterator it
		= ((NUMA_graph_index<scc_vertex> *) index)->begin();
	NUMA_graph_index<scc_vertex>::const_iterator end_it
		= ((NUMA_graph_index<scc_vertex> *) index)->end();
	it = ((NUMA_graph_index<scc_vertex> *) index)->begin();
	size_t max_comp_size = 0;
	size_t num_assigned = 0;
	for (; it != end_it; ++it) {
		const scc_vertex &v = (const scc_vertex &) *it;
		if (v.is_assigned())
			num_assigned++;
		if (v.is_assigned() && v.get_comp_id() == max_v)
			max_comp_size++;
	}
	printf("%ld vertices are assigned to components. max SCC has %ld vertices\n",
			num_assigned, max_comp_size);

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();

	// Compute the summary of the result.
	typedef std::unordered_map<vertex_id_t, size_t> comp_map_t;
	comp_map_t comp_counts;
	it = index->begin();
	for (; it != end_it; ++it) {
		const scc_vertex &v = (const scc_vertex &) *it;
//		if (v.get_num_in_edges() + v.get_num_out_edges() == 0)
//			continue;

		// If a vertex hasn't been assigned to a component, ignore it.
		if (!v.is_assigned())
			continue;

		comp_map_t::iterator map_it = comp_counts.find(v.get_comp_id());
		if (map_it == comp_counts.end()) {
			comp_counts.insert(std::pair<vertex_id_t, size_t>(
						v.get_comp_id(), 1));
		}
		else {
			map_it->second++;
		}
	}
	printf("There are %ld components\n", comp_counts.size());

	// Output the summary of the result.
	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		assert(f);
		BOOST_FOREACH(comp_map_t::value_type &p, comp_counts) {
			if (p.second >= min_comp_size)
				fprintf(f, "component %u: %ld\n", p.first, p.second);
		}
		fclose(f);
	}

	printf("SCC takes %f seconds\n", time_diff(scc_start, end));
}
