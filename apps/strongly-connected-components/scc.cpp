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

class wcc_comp_message: public vertex_message
{
	vertex_id_t id;
	uint64_t color;
public:
	wcc_comp_message(vertex_id_t id, uint64_t color): vertex_message(
			sizeof(wcc_comp_message), true) {
		this->id = id;
		this->color = color;
	}

	uint64_t get_color() const {
		return color;
	}

	vertex_id_t get_id() const {
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

class fwbw_state_t
{
	enum {
		FW_COLOR,
		BW_COLOR,
		FW_BFS,
		BW_BFS,
		ASSIGNED,
		FW_VISITED,
		BW_VISITED,
	};

	static const int COLOR_OFF = 60;

	vertex_id_t base_color;
	vertex_id_t pivot;
	bit_flags<short> flags;
public:
	fwbw_state_t() {
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

	void assign_comp(vertex_id_t comp) {
		pivot = comp;
		flags.set_flag(ASSIGNED);
	}

	void assign_comp() {
		flags.set_flag(ASSIGNED);
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
};

class scc_vertex: public compute_directed_vertex
{
	// for trimming
	vsize_t num_in_edges;
	vsize_t num_out_edges;

	vertex_id_t wcc_neigh_min;
public:
	fwbw_state_t fwbw_state;

	scc_vertex() {
		num_in_edges = 0;
		num_out_edges = 0;
		wcc_neigh_min = UINT_MAX;
	}

	scc_vertex(vertex_id_t id, const vertex_index *index): compute_directed_vertex(
			id, index) {
		num_in_edges = get_num_in_edges();
		num_out_edges = get_num_out_edges();
		wcc_neigh_min = id;
	}

	bool is_assigned() const {
		return fwbw_state.is_assigned();
	}

	vertex_id_t get_comp_id() const {
		return fwbw_state.get_comp_id();
	}

	void init_wcc() {
		wcc_neigh_min = get_id();
		fwbw_state.set_pivot(UINT_MAX);
	}

	void init_fwbw() {
		fwbw_state.set_fw();
		fwbw_state.set_bw();
		fwbw_state.set_pivot(get_id());
	}

	void run(graph_engine &graph) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::TRIM1:
				run_stage_trim1(graph);
				break;
			case scc_stage_t::TRIM2:
				run_stage_trim2(graph);
				break;
			case scc_stage_t::TRIM3:
				run_stage_trim3(graph);
				break;
			case scc_stage_t::FWBW:
				run_stage_FWBW(graph);
				break;
			case scc_stage_t::PARTITION:
				run_stage_part(graph);
				break;
			case scc_stage_t::WCC:
				run_stage_wcc(graph);
				break;
			default:
				assert(0);
		}
	}

	void run_stage_trim1(graph_engine &graph);
	void run_stage_trim2(graph_engine &graph);
	void run_stage_trim3(graph_engine &graph);
	void run_stage_FWBW(graph_engine &graph);
	void run_stage_part(graph_engine &graph);
	void run_stage_wcc(graph_engine &graph);

	void run(graph_engine &graph, const page_vertex &vertex) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::TRIM1:
				run_stage_trim1(graph, vertex);
				break;
			case scc_stage_t::TRIM2:
				run_stage_trim2(graph, vertex);
				break;
			case scc_stage_t::TRIM3:
				run_stage_trim3(graph, vertex);
				break;
			case scc_stage_t::FWBW:
				run_stage_FWBW(graph, vertex);
				break;
			case scc_stage_t::PARTITION:
				run_stage_part(graph, vertex);
				break;
			case scc_stage_t::WCC:
				run_stage_wcc(graph, vertex);
				break;
			default:
				assert(0);
		}
	}

	void run_stage_trim1(graph_engine &graph, const page_vertex &vertex);
	void run_stage_trim2(graph_engine &graph, const page_vertex &vertex);
	void run_stage_trim3(graph_engine &graph, const page_vertex &vertex);
	void run_stage_FWBW(graph_engine &graph, const page_vertex &vertex);
	void run_stage_part(graph_engine &graph, const page_vertex &vertex);
	void run_stage_wcc(graph_engine &graph, const page_vertex &vertex);

	virtual void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
		if (is_assigned())
			return;

		switch(scc_stage) {
			case scc_stage_t::TRIM1:
				run_on_messages_stage_trim1(graph, msgs, num);
				break;
			case scc_stage_t::TRIM2:
				run_on_messages_stage_trim2(graph, msgs, num);
				break;
			case scc_stage_t::TRIM3:
				run_on_messages_stage_trim3(graph, msgs, num);
				break;
			case scc_stage_t::FWBW:
				run_on_messages_stage_FWBW(graph, msgs, num);
				break;
			case scc_stage_t::PARTITION:
				run_on_messages_stage_part(graph, msgs, num);
				break;
			case scc_stage_t::WCC:
				run_on_messages_stage_wcc(graph, msgs, num);
				break;
			default:
				assert(0);
		}
	}

	void run_on_messages_stage_trim1(graph_engine &graph,
			const vertex_message *msgs[], int num);
	void run_on_messages_stage_trim2(graph_engine &graph,
			const vertex_message *msgs[], int num);
	void run_on_messages_stage_trim3(graph_engine &graph,
			const vertex_message *msgs[], int num);
	void run_on_messages_stage_FWBW(graph_engine &graph,
			const vertex_message *msgs[], int num);
	void run_on_messages_stage_part(graph_engine &graph,
			const vertex_message *msgs[], int num);
	void run_on_messages_stage_wcc(graph_engine &graph,
			const vertex_message *msgs[], int num);
};

void scc_vertex::run_stage_trim1(graph_engine &graph)
{
	if (num_in_edges == 0 || num_out_edges == 0) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);

		// This vertex has to be a SCC itself.
		fwbw_state.assign_comp(id);
	}
}

void scc_vertex::run_stage_trim1(graph_engine &graph, const page_vertex &vertex)
{
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(BOTH_EDGES);
	stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(BOTH_EDGES));
	int num_dests = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(BOTH_EDGES); it != end_it; ++it) {
		vertex_id_t id = *it;
		dest_buf[num_dests++] = id;
	}

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
		graph.multicast_msg(dest_buf.data(), num_dests, msg);
	}
}

void scc_vertex::run_on_messages_stage_trim1(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
	for (int i = 0; i < num; i++) {
		const trim1_message *msg = (const trim1_message *) msgs[i];
		switch(msg->get_type()) {
			case edge_type::IN_EDGE:
				assert(num_in_edges > 0);
				num_in_edges--;
				break;
			case edge_type::OUT_EDGE:
				assert(num_out_edges > 0);
				num_out_edges--;
				break;
			default:
				assert(0);
		}
	}
}

void scc_vertex::run_stage_trim2(graph_engine &graph)
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

void scc_vertex::run_stage_trim2(graph_engine &graph, const page_vertex &vertex)
{
	if (vertex.get_id() == get_id()) {
		// Ideally, we should use the remaining in-edges or out-edges,
		// but we don't know which edges have been removed, so we just
		// use the original number of edges.
		if (get_num_in_edges() == 1) {
			page_byte_array::const_iterator<vertex_id_t> it
				= vertex.get_neigh_begin(edge_type::IN_EDGE);
			vertex_id_t neighbor = *it;
			// If the only in-edge is to itself, it's a SCC itself.
			if (neighbor == get_id()) {
				fwbw_state.assign_comp(get_id());
				trim2_vertices++;
			}
			else {
				scc_vertex &neigh_v = (scc_vertex &) graph.get_vertex(neighbor);
				// If the vertex's out-edge list contains the neighbor,
				// it means the neighbor's only in-edge connect to this vertex.
				if (fwbw_state.get_color() == neigh_v.fwbw_state.get_color()
						&& get_id() < neighbor
						&& neigh_v.get_num_in_edges() == 1
						&& vertex.contain_edge(edge_type::OUT_EDGE, neighbor)) {
					fwbw_state.assign_comp(get_id());
					trim2_message msg(get_id());
					graph.send_msg(neighbor, msg);
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
				fwbw_state.assign_comp(get_id());
				trim2_vertices++;
			}
			else {
				scc_vertex &neigh_v = (scc_vertex &) graph.get_vertex(neighbor);
				// The same as above.
				if (fwbw_state.get_color() == neigh_v.fwbw_state.get_color()
						&& get_id() < neighbor
						&& neigh_v.get_num_out_edges() == 1
						&& vertex.contain_edge(edge_type::IN_EDGE, neighbor)) {
					fwbw_state.assign_comp(get_id());
					trim2_message msg(get_id());
					graph.send_msg(neighbor, msg);
					trim2_vertices += 2;
				}
			}
		}
		else
			assert(0);
	}
}

void scc_vertex::run_on_messages_stage_trim2(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
	for (int i = 0; i < num; i++) {
		const trim2_message *msg = (const trim2_message *) msgs[i];
		fwbw_state.assign_comp(msg->get_comp_id());
	}
}

void scc_vertex::run_stage_trim3(graph_engine &graph)
{
	vertex_id_t id = get_id();
	request_vertices(&id, 1);
}

std::atomic_long trim3_vertices;

void scc_vertex::run_stage_trim3(graph_engine &graph, const page_vertex &vertex)
{
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(IN_EDGE);
	stack_array<vertex_id_t, 1024> in_neighs(vertex.get_num_edges(IN_EDGE));
	int num_in_neighs = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(IN_EDGE); it != end_it; ++it) {
		vertex_id_t id = *it;
		scc_vertex &neigh = (scc_vertex &) graph.get_vertex(id);
		// We should ignore the neighbors that has been assigned to a component.
		// or has a different color.
		if (neigh.is_assigned()
				|| neigh.fwbw_state.get_color() != fwbw_state.get_color())
			continue;

		in_neighs[num_in_neighs++] = id;
	}

	end_it = vertex.get_neigh_end(OUT_EDGE);
	stack_array<vertex_id_t, 1024> out_neighs(vertex.get_num_edges(OUT_EDGE));
	int num_out_neighs = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(OUT_EDGE); it != end_it; ++it) {
		vertex_id_t id = *it;
		scc_vertex &neigh = (scc_vertex &) graph.get_vertex(id);
		// We should ignore the neighbors that has been assigned to a component.
		// or has a different color.
		if (neigh.is_assigned()
				|| neigh.fwbw_state.get_color() != fwbw_state.get_color())
			continue;

		out_neighs[num_out_neighs++] = id;
	}

	if (num_in_neighs == 0 || num_out_neighs == 0) {
		trim3_vertices++;
		// This vertex has been isolated, it can assign to a SCC now.
		fwbw_state.assign_comp(get_id());
		if (num_in_neighs > 0)
			graph.activate_vertices(in_neighs.data(), num_in_neighs);
		if (num_out_neighs > 0)
			graph.activate_vertices(out_neighs.data(), num_out_neighs);
	}
}

void scc_vertex::run_on_messages_stage_trim3(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
}

void scc_vertex::run_stage_FWBW(graph_engine &graph)
{
	// If the vertex has been visited in both directions,
	// we don't need to do anything.
	if (fwbw_state.has_fw_visited() && fwbw_state.has_bw_visited())
		return;
	// If the vertex has been visisted in forward direction, and it doesn't
	// need to visit other in the backwoard direction, then we don't need to
	// do anything.
	if (fwbw_state.has_fw_visited() && !fwbw_state.is_bw())
		return;
	// The same for the other direction.
	if (fwbw_state.has_bw_visited() && !fwbw_state.is_fw())
		return;

	// It's possible that the vertex is activated by another vertex of
	// a different color. If that is the case, the vertex may not have
	// the forward BFS flag nor the backward BFS flag. Do nothing.
	if (!fwbw_state.is_bw() && !fwbw_state.is_fw())
		return;

	vertex_id_t id = get_id();
	request_vertices(&id, 1);
}

void scc_vertex::run_stage_FWBW(graph_engine &graph, const page_vertex &vertex)
{
	bool do_some = false;

	if (fwbw_state.is_bw()) {
		do_some = true;
		fwbw_state.set_bw_visited();
		page_byte_array::const_iterator<vertex_id_t> end_it
			= vertex.get_neigh_end(IN_EDGE);
		stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(IN_EDGE));
		int num_dests = 0;
		for (page_byte_array::const_iterator<vertex_id_t> it
				= vertex.get_neigh_begin(IN_EDGE); it != end_it; ++it) {
			vertex_id_t id = *it;
			dest_buf[num_dests++] = id;
		}
		fwbw_message msg(fwbw_state.get_color(), fwbw_state.get_pivot(), false);
		graph.multicast_msg(dest_buf.data(), num_dests, msg);
	}

	if (fwbw_state.is_fw()) {
		do_some = true;
		fwbw_state.set_fw_visited();
		page_byte_array::const_iterator<vertex_id_t> end_it
			= vertex.get_neigh_end(OUT_EDGE);
		stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(OUT_EDGE));
		int num_dests = 0;
		for (page_byte_array::const_iterator<vertex_id_t> it
				= vertex.get_neigh_begin(OUT_EDGE); it != end_it; ++it) {
			vertex_id_t id = *it;
			dest_buf[num_dests++] = id;
		}
		fwbw_message msg(fwbw_state.get_color(), fwbw_state.get_pivot(), true);
		graph.multicast_msg(dest_buf.data(), num_dests, msg);
	}
	assert(do_some);
}

void scc_vertex::run_on_messages_stage_FWBW(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
	uint64_t color = fwbw_state.get_color();
	for (int i = 0; i < num; i++) {
		const fwbw_message *msg = (const fwbw_message *) msgs[i];
		// If the current vertex has a different color, it means it's in
		// a different partition. The vertex can just ignore the message.
		if (msg->get_color() != color)
			continue;

		fwbw_state.set_pivot(msg->get_pivot());
		if (msg->is_forward())
			fwbw_state.set_fw();
		else
			fwbw_state.set_bw();
	}
}

void scc_vertex::run_stage_part(graph_engine &graph)
{
	if (fwbw_state.is_fw() && fwbw_state.is_bw())
		fwbw_state.assign_comp();
	else if (fwbw_state.is_fw())
		fwbw_state.assign_new_fw_color();
	else if (fwbw_state.is_bw())
		fwbw_state.assign_new_bw_color();
	fwbw_state.clear_flags();
}

void scc_vertex::run_stage_part(graph_engine &graph, const page_vertex &vertex)
{
}

void scc_vertex::run_on_messages_stage_part(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
}

void scc_vertex::run_stage_wcc(graph_engine &graph)
{
	if (wcc_neigh_min < fwbw_state.get_pivot()) {
		fwbw_state.set_pivot(wcc_neigh_min);
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}
}

void scc_vertex::run_stage_wcc(graph_engine &graph, const page_vertex &vertex)
{
	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	page_byte_array::const_iterator<vertex_id_t> end_it
		= vertex.get_neigh_end(BOTH_EDGES);
	stack_array<vertex_id_t, 1024> dest_buf(vertex.get_num_edges(BOTH_EDGES));
	int num_dests = 0;
	for (page_byte_array::const_iterator<vertex_id_t> it
			= vertex.get_neigh_begin(BOTH_EDGES); it != end_it; ++it) {
		vertex_id_t id = *it;
		dest_buf[num_dests++] = id;
	}
	wcc_comp_message msg(fwbw_state.get_pivot(), fwbw_state.get_color());
	graph.multicast_msg(dest_buf.data(), num_dests, msg);
}

void scc_vertex::run_on_messages_stage_wcc(graph_engine &graph,
		const vertex_message *msgs[], int num)
{
	for (int i = 0; i < num; i++) {
		wcc_comp_message *msg = (wcc_comp_message *) msgs[i];
		// If the current vertex has a different color, it means it's in
		// a different partition. The vertex can just ignore the message.
		if (msg->get_color() != fwbw_state.get_color())
			continue;
		wcc_neigh_min = min(wcc_neigh_min, msg->get_id());
	}
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

class scc_filter: public vertex_filter
{
public:
	virtual bool keep(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		return !scc_v.is_assigned();
	}
};

class wcc_filter: public vertex_filter
{
	std::atomic_ulong count;
public:
	wcc_filter() {
		count = 0;
	}

	virtual bool keep(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		bool activate = !scc_v.is_assigned();
		if (activate) {
			scc_v.init_wcc();
			count.fetch_add(1, std::memory_order_relaxed);
		}
		return activate;
	}

	unsigned long get_count() const {
		return count;
	}
};

class sec_fwbw_filter: public vertex_filter
{
public:
	virtual bool keep(compute_vertex &v) {
		scc_vertex &scc_v = (scc_vertex &) v;
		bool activate = !scc_v.is_assigned()
			&& scc_v.fwbw_state.get_pivot() == scc_v.get_id();
		// If the vertex hasn't been assigned to a component,
		// let's use the result of wcc (which is stored in pivot) as the color
		if (!scc_v.is_assigned())
			scc_v.fwbw_state.assign_new_color(scc_v.fwbw_state.get_pivot());
		if (activate)
			scc_v.init_fwbw();
		return activate;
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
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

	NUMA_graph_index<scc_vertex>::const_iterator it
		= ((NUMA_graph_index<scc_vertex> *) index)->begin();
	NUMA_graph_index<scc_vertex>::const_iterator end_it
		= ((NUMA_graph_index<scc_vertex> *) index)->end();
	size_t num_assigned = 0;
	vertex_id_t max_v = -1;
	vertex_id_t max_degree = 0;
	for (; it != end_it; ++it) {
		const scc_vertex &v = (const scc_vertex &) *it;
		if (v.is_assigned())
			num_assigned++;
		if (v.get_num_in_edges() + v.get_num_out_edges() > max_degree) {
			max_degree = v.get_num_in_edges() + v.get_num_out_edges();
			max_v = v.get_id();
		}
	}
	printf("trim1 takes %f seconds. It trims %ld vertices\n",
			time_diff(start, end), num_assigned);

	scc_stage = scc_stage_t::TRIM2;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("trim2 takes %f seconds. It trims %ld vertices\n",
			time_diff(start, end),
			trim2_vertices.load(std::memory_order_relaxed));

	scc_stage = scc_stage_t::FWBW;
	scc_vertex &v = (scc_vertex &) index->get_vertex(max_v);
	printf("the remaining max vertex: %u (%d, %d)\n", max_v,
			v.get_num_in_edges(), v.get_num_out_edges());
	v.init_fwbw();
	gettimeofday(&start, NULL);
	graph->start(&max_v, 1);
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("FWBW takes %f seconds\n", time_diff(start, end));

	scc_stage = scc_stage_t::PARTITION;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("partition takes %f seconds\n", time_diff(start, end));

	while (true) {
		scc_stage = scc_stage_t::TRIM3;
		trim3_vertices = 0;
		gettimeofday(&start, NULL);
		graph->start(std::shared_ptr<vertex_filter>(new scc_filter()));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("trim3 takes %f seconds, and trime %ld vertices\n",
				time_diff(start, end),
				trim3_vertices.load(std::memory_order_relaxed));

		std::shared_ptr<vertex_filter> wfilter
			= std::shared_ptr<vertex_filter>(new wcc_filter());
		scc_stage = scc_stage_t::WCC;
		gettimeofday(&start, NULL);
		graph->start(wfilter);
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("WCC takes %f seconds\n", time_diff(start, end));
		size_t count = std::static_pointer_cast<wcc_filter, vertex_filter>(
				wfilter)->get_count();
		if (count == 0)
			break;
		printf("There are %ld vertices at the stage of WCC\n", count);

		scc_stage = scc_stage_t::FWBW;
		gettimeofday(&start, NULL);
		graph->start(std::shared_ptr<vertex_filter>(new sec_fwbw_filter()));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("FWBW takes %f seconds\n", time_diff(start, end));

		scc_stage = scc_stage_t::PARTITION;
		gettimeofday(&start, NULL);
		graph->start(std::shared_ptr<vertex_filter>(new scc_filter()));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("partition takes %f seconds\n", time_diff(start, end));
	}

	it = ((NUMA_graph_index<scc_vertex> *) index)->begin();
	size_t max_comp_size = 0;
	num_assigned = 0;
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
	FILE *f = stdout;
	if (!output_file.empty()) {
		f = fopen(output_file.c_str(), "w");
		assert(f);
	}
	BOOST_FOREACH(comp_map_t::value_type &p, comp_counts) {
		if (p.second >= min_comp_size)
			fprintf(f, "component %u: %ld\n", p.first, p.second);
	}
	if (!output_file.empty()) {
		fclose(f);
	}

	printf("SCC takes %f seconds\n", time_diff(scc_start, end));
}
