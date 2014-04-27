/**
 * Copyright 2014 Da Zheng
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

#include <set>
#include <vector>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"

#include "graphlab/cuckoo_set_pow2.hpp"
#include "scan_graph.h"

const double BIN_SEARCH_RATIO = 100;

atomic_number<long> num_working_vertices;
atomic_number<long> num_completed_vertices;

size_t neighbor_list::count_edges_hash(const page_vertex *v,
		page_byte_array::const_iterator<vertex_id_t> other_it,
		page_byte_array::const_iterator<vertex_id_t> other_end,
		std::vector<vertex_id_t> *common_neighs)
{
	size_t num_local_edges = 0;

	while (other_it != other_end) {
		vertex_id_t neigh_neighbor = *other_it;
		if (neigh_neighbor != v->get_id()
				&& neigh_neighbor != this->get_id()) {
			if (this->contains(neigh_neighbor)) {
#if 0
				num_local_edges += (*other_data_it).get_count();
#endif
				num_local_edges++;
				if (common_neighs)
					common_neighs->push_back(neigh_neighbor);
			}
		}
		++other_it;
	}
	return num_local_edges;
}

#if 0
int neighbor_list::count_edges_bin_search_this(const page_vertex *v,
		std::vector<attributed_neighbor>::const_iterator this_it,
		std::vector<attributed_neighbor>::const_iterator this_end,
		page_byte_array::const_iterator<vertex_id_t> other_it,
		page_byte_array::const_iterator<vertex_id_t> other_end,
		std::vector<vertex_id_t> &common_neighs)
{
	int num_local_edges = 0;
	int num_v_edges = other_end - other_it;
	int size_log2 = log2(neighbors->size());
	num_rand_jumps += size_log2 * num_v_edges;
	scan_bytes += num_v_edges * sizeof(vertex_id_t);
	while (other_it != other_end) {
		vertex_id_t neigh_neighbor = *other_it;
		if (neigh_neighbor != v->get_id()
				&& neigh_neighbor != this->get_id()) {
			std::vector<attributed_neighbor>::const_iterator first
				= std::lower_bound(this_it, this_end,
						attributed_neighbor(neigh_neighbor), comp_edge());
			if (first != this_end && neigh_neighbor == first->get_id()) {
#if 0
				num_local_edges += (*other_data_it).get_count();
#endif
				num_local_edges++;
				common_neighs.push_back(first->get_id());
			}
		}
		++other_it;
#if 0
		++other_data_it;
#endif
	}
	return num_local_edges;
}
#endif

size_t neighbor_list::count_edges_bin_search_other(const page_vertex *v,
		neighbor_list::id_iterator this_it,
		neighbor_list::id_iterator this_end,
		page_byte_array::const_iterator<vertex_id_t> other_it,
		page_byte_array::const_iterator<vertex_id_t> other_end,
		std::vector<vertex_id_t> *common_neighs)
{
	size_t num_local_edges = 0;

	for (; this_it != this_end; this_it++) {
		vertex_id_t this_neighbor = *this_it;
		// We need to skip loops.
		if (this_neighbor == v->get_id()
				|| this_neighbor == this->get_id()) {
			continue;
		}

		page_byte_array::const_iterator<vertex_id_t> first
			= std::lower_bound(other_it, other_end, this_neighbor);
		// found it.
		if (first != other_end && !(this_neighbor < *first)) {
			int num_dups = 0;
			do {
#if 0
				page_byte_array::const_iterator<edge_count> data_it
					= other_data_it;
				data_it += first - other_it;
				// Edges in the v's neighbor lists may duplicated.
				// The duplicated neighbors need to be counted
				// multiple times.
				num_local_edges += (*data_it).get_count();
				++data_it;
#endif
				num_dups++;
				num_local_edges++;
				++first;
			} while (first != other_end && this_neighbor == *first);
			if (common_neighs)
				common_neighs->push_back(this_neighbor);
		}
	}
	return num_local_edges;
}

size_t neighbor_list::count_edges_scan(const page_vertex *v,
		neighbor_list::id_iterator this_it,
		neighbor_list::id_iterator this_end,
		page_byte_array::seq_const_iterator<vertex_id_t> other_it,
		std::vector<vertex_id_t> *common_neighs)
{
	size_t num_local_edges = 0;
	while (other_it.has_next() && this_it != this_end) {
		vertex_id_t this_neighbor = *this_it;
		vertex_id_t neigh_neighbor = other_it.curr();
		if (neigh_neighbor == v->get_id()
				|| neigh_neighbor == this->get_id()) {
			other_it.next();
#if 0
			++other_data_it;
#endif
			continue;
		}
		if (this_neighbor == neigh_neighbor) {
			if (common_neighs)
				common_neighs->push_back(*this_it);
			do {
				// Edges in the v's neighbor lists may duplicated.
				// The duplicated neighbors need to be counted
				// multiple times.
#if 0
				num_local_edges += (*other_data_it).get_count();
				++other_data_it;
#endif
				num_local_edges++;
				other_it.next();
			} while (other_it.has_next() && this_neighbor == other_it.curr());
			++this_it;
		}
		else if (this_neighbor < neigh_neighbor) {
			++this_it;
		}
		else {
			other_it.next();
#if 0
			++other_data_it;
#endif
		}
	}
	return num_local_edges;
}

size_t neighbor_list::count_edges(const page_vertex *v, edge_type type,
		std::vector<vertex_id_t> *common_neighs)
{
	size_t num_v_edges = v->get_num_edges(type);
	if (num_v_edges == 0)
		return 0;

#ifdef PV_STAT
	min_comps += min(num_v_edges, this->size());
#endif
	page_byte_array::const_iterator<vertex_id_t> other_it
		= v->get_neigh_begin(type);
#if 0
	page_byte_array::const_iterator<edge_count> other_data_it
		= v->get_edge_data_begin<edge_count>(type);
#endif
	page_byte_array::const_iterator<vertex_id_t> other_end
		= std::lower_bound(other_it, v->get_neigh_end(type),
				v->get_id());
	num_v_edges = other_end - other_it;
	if (num_v_edges == 0)
		return 0;

	neighbor_list::id_iterator this_it = this->get_id_begin();
	neighbor_list::id_iterator this_end = this->get_id_end();
	this_end = std::lower_bound(this_it, this_end,
			v->get_id());

	if (num_v_edges / this->size() > BIN_SEARCH_RATIO) {
#ifdef PV_STAT
		int size_log2 = log2(num_v_edges);
		scan_bytes += this->size() * sizeof(vertex_id_t);
		rand_jumps += size_log2 * this->size();
#endif
		return count_edges_bin_search_other(v, this_it, this_end,
				other_it, other_end, common_neighs);
	}
	else if (this->size() / num_v_edges > 16) {
#ifdef PV_STAT
		scan_bytes += num_v_edges * sizeof(vertex_id_t);
		rand_jumps += num_v_edges;
#endif
		return count_edges_hash(v, other_it, other_end, common_neighs);
	}
	else {
#ifdef PV_STAT
		scan_bytes += num_v_edges * sizeof(vertex_id_t);
		scan_bytes += this->size() * sizeof(vertex_id_t);
#endif
		return count_edges_scan(v, this_it, this_end,
				v->get_neigh_seq_it(type, 0, num_v_edges), common_neighs);
	}
}

size_t neighbor_list::count_edges(const page_vertex *v)
{
	assert(!this->empty());
	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0)
		return 0;

	size_t ret = count_edges(v, edge_type::IN_EDGE, NULL)
		+ count_edges(v, edge_type::OUT_EDGE, NULL);
	return ret;
}

void scan_vertex::run_on_itself(vertex_program &prog, const page_vertex &vertex)
{
	assert(!local_value.has_real_local());
	assert(!local_value.has_runtime_data());

	size_t num_local_edges = vertex.get_num_edges(edge_type::BOTH_EDGES);
	assert(num_local_edges == get_num_edges());
#ifdef PV_STAT
	num_all_edges = num_local_edges;
#endif
	if (num_local_edges == 0)
		return;

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);

	runtime_data_t *local_data = create_runtime(prog.get_graph(), *this, vertex);
	local_value.set_runtime_data(local_data);
#ifdef PV_STAT
	gettimeofday(&vertex_start, NULL);
	fprintf(stderr, "compute v%u (with %d edges, compute on %ld edges, potential %ld inter-edges) on thread %d at %.f seconds\n",
			get_id(), num_all_edges, get_runtime_data()->neighbors->size(),
			get_est_local_scan(graph, &vertex), thread::get_curr_thread()->get_id(),
			time_diff(graph_start, vertex_start));
#endif

	page_byte_array::const_iterator<vertex_id_t> it = vertex.get_neigh_begin(
			edge_type::BOTH_EDGES);
#if 0
	page_byte_array::const_iterator<edge_count> data_it
		= vertex->get_edge_data_begin<edge_count>(
				edge_type::BOTH_EDGES);
#endif
	page_byte_array::const_iterator<vertex_id_t> end = vertex.get_neigh_end(
			edge_type::BOTH_EDGES);
#if 0
	for (; it != end; ++it, ++data_it) {
#endif
	size_t tmp = 0;
	for (; it != end; ++it) {
		vertex_id_t id = *it;
		// Ignore loops
		if (id != vertex.get_id()) {
#if 0
			num_edges.inc((*data_it).get_count());
#endif
			tmp++;
		}
	}
	local_data->local_scan += tmp;

	if (local_data->neighbors->empty()) {
		destroy_runtime(*this, local_data);
		local_value.set_real_local(0);
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	std::vector<vertex_id_t> neighbors;
	local_data->neighbors->get_neighbors(neighbors);
	request_vertices(neighbors.data(), neighbors.size());
}

void scan_vertex::run_on_neighbor(vertex_program &prog, const page_vertex &vertex)
{
	assert(local_value.has_runtime_data());
	runtime_data_t *local_data = local_value.get_runtime_data();
	local_data->num_joined++;
#ifdef PV_STAT
	struct timeval start, end;
	gettimeofday(&start, NULL);
#endif
	size_t ret = local_data->neighbors->count_edges(&vertex);
	if (ret > 0)
		local_data->local_scan += ret;
#ifdef PV_STAT
	gettimeofday(&end, NULL);
	time_us += time_diff_us(start, end);
#endif

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (local_data->num_joined == local_data->neighbors->size()) {
		local_value.set_real_local(local_data->local_scan);

		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

#ifdef PV_STAT
		struct timeval curr;
		gettimeofday(&curr, NULL);
		fprintf(stderr,
				"v%u: # edges: %d, scan: %d, scan bytes: %ld, # rand jumps: %ld, # comps: %ld, time: %ldms\n",
				get_id(), num_all_edges, local_data->local_scan, scan_bytes, rand_jumps, min_comps, time_us / 1000);
#endif

		::finding_triangles_end(prog, *this, local_data);

		destroy_runtime(*this, local_data);
	}
}

class skip_self {
	vertex_id_t id;
	graph_engine &graph;
public:
	skip_self(graph_engine &_graph, vertex_id_t id): graph(_graph) {
		this->id = id;
	}

	bool operator()(attributed_neighbor &e) {
		return operator()(e.get_id());
	}

	bool operator()(vertex_id_t id) {
		return this->id == id;
	}
};

class merge_edge
{
public:
	attributed_neighbor operator()(const attributed_neighbor &e1,
			const attributed_neighbor &e2) {
		assert(e1.get_id() == e2.get_id());
		return attributed_neighbor(e1.get_id(),
				e1.get_num_dups() + e2.get_num_dups());
	}
};

runtime_data_t *default_create_runtime(graph_engine &graph, scan_vertex &scan_v,
		const page_vertex &pg_v)
{
	merge_edge merge;
	std::vector<attributed_neighbor> neighbors(
			pg_v.get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			pg_v.get_neigh_begin(edge_type::IN_EDGE),
			pg_v.get_neigh_end(edge_type::IN_EDGE),
			pg_v.get_neigh_begin(edge_type::OUT_EDGE),
			pg_v.get_neigh_end(edge_type::OUT_EDGE),
			skip_self(graph, pg_v.get_id()), merge,
			neighbors.begin());
	neighbors.resize(num_neighbors);
	return new runtime_data_t(std::unique_ptr<neighbor_list>(
				new neighbor_list(pg_v, neighbors)));
}

void default_destroy_runtime(scan_vertex &graph, runtime_data_t *data)
{
	delete data;
}

void (*finding_triangles_end)(vertex_program &, scan_vertex &, runtime_data_t *);
runtime_data_t *(*create_runtime)(graph_engine &, scan_vertex &,
		const page_vertex &) = default_create_runtime;
void (*destroy_runtime)(scan_vertex &,
		runtime_data_t *) = default_destroy_runtime;
