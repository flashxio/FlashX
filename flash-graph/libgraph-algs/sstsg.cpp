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

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <set>
#include <vector>

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"
#include "FG_vector.h"
#include "ts_graph.h"

using namespace fg;

namespace {

time_t timestamp;
time_t time_interval = 1;
int num_time_intervals = 1;

class scan_vertex: public compute_vertex
{
	// The number of vertices that have joined with the vertex.
	int num_joined;
	// Local scan in its neighborhood in different timestamps.
	std::vector<size_t> *local_scans;
	// All neighbors (in both in-edges and out-edges)
	// in the specified timestamp.
	std::vector<vertex_id_t> *neighbors;

	// The final result.
	double result;
public:
	scan_vertex(vertex_id_t id): compute_vertex(id) {
		num_joined = 0;
		local_scans = NULL;
		neighbors = NULL;
	}

	double get_result() const {
		return result;
	}

	size_t count_edges(vertex_program &prog, const page_directed_vertex &v,
			const std::vector<vertex_id_t> *neighbors, time_t timestamp,
			time_t time_interval);
	size_t count_edges(vertex_program &prog, const page_directed_vertex &v,
			const std::vector<vertex_id_t> *neighbors, time_t timestamp,
			time_t time_interval, edge_type type);

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == prog.get_vertex_id(*this))
			run_on_itself(prog, (const page_directed_vertex &) vertex);
		else
			run_on_neighbor(prog, (const page_directed_vertex &) vertex);
	}

	void run_on_itself(vertex_program &prog, const page_directed_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_directed_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

size_t scan_vertex::count_edges(vertex_program &prog, const page_directed_vertex &v,
		const std::vector<vertex_id_t> *neighbors, time_t timestamp,
		time_t time_interval, edge_type type)
{
	size_t num_local_edges = 0;
	edge_seq_iterator it = get_ts_iterator(v, type, timestamp, time_interval);
	// If there are no edges in the time interval.
	if (it.get_num_tot_entries() == 0)
		return 0;

	std::vector<vertex_id_t>::const_iterator this_it = neighbors->begin();
	std::vector<vertex_id_t>::const_iterator this_end
		= std::lower_bound(this_it, neighbors->end(), v.get_id());
	if (this_it == this_end)
		return 0;

	PAGE_FOREACH(vertex_id_t, neigh_neighbor, it) {
		if (neigh_neighbor != v.get_id()
				&& neigh_neighbor != prog.get_vertex_id(*this)) {
			if (std::binary_search(this_it, this_end, neigh_neighbor))
				num_local_edges++;
		}
	} PAGE_FOREACH_END
	return num_local_edges;
}

size_t scan_vertex::count_edges(vertex_program &prog, const page_directed_vertex &v,
		const std::vector<vertex_id_t> *neighbors, time_t timestamp,
		time_t time_interval)
{
	assert(!neighbors->empty());
	return count_edges(prog, v, neighbors, timestamp, time_interval,
			edge_type::IN_EDGE)
		+ count_edges(prog, v, neighbors, timestamp, time_interval,
				edge_type::OUT_EDGE);
}

template<class InputIterator1, class InputIterator2, class Skipper,
	class OutputIterator>
int unique_merge(InputIterator1 it1, InputIterator1 last1,
		InputIterator2 it2, InputIterator2 last2, Skipper skip,
		OutputIterator result)
{
	OutputIterator result_begin = result;
	while (it1 != last1 && it2 != last2) {
		if (*it1 > *it2) {
			typename std::iterator_traits<InputIterator2>::value_type v = *it2;
			if (!skip(v))
				*(result++) = v;
			while (it2 != last2 && *it2 == v)
				++it2;
		}
		else if (*it1 < *it2) {
			typename std::iterator_traits<InputIterator1>::value_type v = *it1;
			if (!skip(v))
				*(result++) = v;
			while (it1 != last1 && *it1 == v)
				++it1;
		}
		else {
			typename std::iterator_traits<InputIterator1>::value_type v = *it1;
			if (!skip(v))
				*(result++) = v;
			while (it1 != last1 && *it1 == v)
				++it1;
			while (it2 != last2 && *it2 == v)
				++it2;
		}
	}

	while (it1 != last1) {
		typename std::iterator_traits<InputIterator1>::value_type v = *it1;
		if (!skip(v))
			*(result++) = v;
		while (it1 != last1 && *it1 == v)
			++it1;
	}

	while (it2 != last2) {
		typename std::iterator_traits<InputIterator2>::value_type v = *it2;
		if (!skip(v))
			*(result++) = v;
		while (it2 != last2 && *it2 == v)
			++it2;
	}
	return result - result_begin;
}

size_t get_neighbors(const page_directed_vertex &v, edge_type type, time_t time_start,
		time_t time_interval, std::vector<vertex_id_t> &neighbors)
{
	edge_seq_iterator it = get_ts_iterator(v, type, time_start, time_interval);
	size_t ret = it.get_num_tot_entries();
	PAGE_FOREACH(vertex_id_t, id, it) {
		neighbors.push_back(id);
	} PAGE_FOREACH_END
	return ret;
}

void scan_vertex::run_on_itself(vertex_program &prog,
		const page_directed_vertex &vertex)
{
	assert(neighbors == NULL);
	assert(num_joined == 0);

	std::vector<vertex_id_t> in_neighbors;
	std::vector<vertex_id_t> out_neighbors;
	get_neighbors(vertex, edge_type::IN_EDGE, timestamp, time_interval,
			in_neighbors);
	get_neighbors(vertex, edge_type::OUT_EDGE, timestamp, time_interval,
			out_neighbors);
	std::sort(in_neighbors.begin(), in_neighbors.end());
	std::sort(out_neighbors.begin(), out_neighbors.end());
	if (in_neighbors.size() + out_neighbors.size() == 0)
		return;

	local_scans = new std::vector<size_t>(num_time_intervals);
	neighbors = new std::vector<vertex_id_t>(
			in_neighbors.size() + out_neighbors.size());

	class skip_self {
		vertex_id_t id;
	public:
		skip_self(vertex_id_t id) {
			this->id = id;
		}

		bool operator()(vertex_id_t id) {
			return this->id == id;
		}
	};

	vertex_id_t curr_id = prog.get_vertex_id(*this);
	int num_neighbors = unique_merge(
			in_neighbors.begin(), in_neighbors.end(),
			out_neighbors.begin(), out_neighbors.end(),
			skip_self(vertex.get_id()),
			neighbors->begin());
	neighbors->resize(num_neighbors);

	size_t lscan = 0;
	BOOST_FOREACH(vertex_id_t id, in_neighbors) {
		if (id != curr_id)
			lscan++;
	}
	BOOST_FOREACH(vertex_id_t id, out_neighbors) {
		if (id != curr_id)
			lscan++;
	}
	local_scans->at(0) = lscan;

	if (neighbors->size() == 0) {
		delete local_scans;
		delete neighbors;
		local_scans = NULL;
		neighbors = NULL;
		return;
	}

	// Count the degree of the vertex in other time intervals.
	for (int ts_idx = 1; ts_idx < num_time_intervals
			&& timestamp >= ts_idx * time_interval; ts_idx++) {
		time_t timestamp2 = timestamp - ts_idx * time_interval;

		// For in-edges.
		edge_seq_iterator it = get_ts_iterator(vertex, edge_type::IN_EDGE,
				timestamp2, time_interval);
		PAGE_FOREACH(vertex_id_t, id, it) {
			// Ignore loop
			if (id != vertex.get_id()
					// The neighbor needs to exist in timestamp 1.
					&& std::binary_search(neighbors->begin(),
						neighbors->end(), id)) {
				local_scans->at(ts_idx)++;
			}
		} PAGE_FOREACH_END

		// For out-edges.
		it = get_ts_iterator(vertex, edge_type::OUT_EDGE, timestamp2,
					time_interval);
		PAGE_FOREACH(vertex_id_t, id, it) {
			// Ignore loop
			if (id != vertex.get_id()
					// The neighbor needs to exist in timestamp 1.
					&& std::binary_search(neighbors->begin(),
						neighbors->end(), id)) {
				local_scans->at(ts_idx)++;
			}
		} PAGE_FOREACH_END
	}

	request_vertices(neighbors->data(), neighbors->size());
}

void scan_vertex::run_on_neighbor(vertex_program &prog,
		const page_directed_vertex &vertex)
{
	num_joined++;
	assert(neighbors);
	for (int j = 0; j < num_time_intervals
			&& timestamp >= j * time_interval; j++) {
		size_t ret = count_edges(prog, vertex, neighbors,
				timestamp - j * time_interval, time_interval);
		if (ret > 0)
			local_scans->at(j) += ret;
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == (int) neighbors->size()) {
		double avg = 0;
		if (num_time_intervals - 1 > 0) {
			double sum = 0;
			for (int i = 1; i < num_time_intervals; i++)
				sum += local_scans->at(i);
			avg = sum / (num_time_intervals - 1);
		}

		double deviation;
		if (num_time_intervals - 1 <= 1)
			deviation = 1;
		else {
			double sum = 0;
			for (int i = 1; i < num_time_intervals; i++) {
				sum += (local_scans->at(i)
						- avg) * (local_scans->at(i) - avg);
			}
			sum = sum / (num_time_intervals - 2);
			deviation = sqrt(sum);
			if (deviation < 1)
				deviation = 1;
		}
		result = (local_scans->at(0) - avg) / deviation;

		delete local_scans;
		delete neighbors;
		local_scans = NULL;
		neighbors = NULL;
	}
}

}

#include "save_result.h"

namespace fg
{

FG_vector<float>::ptr compute_sstsg(FG_graph::ptr fg, time_t start_time,
		time_t interval, int num_intervals)
{
	timestamp = start_time;
	time_interval = interval;
	num_time_intervals = num_intervals;

	graph_index::ptr index = NUMA_graph_index<scan_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
	assert(graph->get_graph_header().get_graph_type() == graph_type::DIRECTED);
	assert(graph->get_graph_header().has_edge_data());
	BOOST_LOG_TRIVIAL(info)
		<< boost::format("scan statistics starts, start: %1%, interval: %2%, #interval: %3%")
		% timestamp % time_interval % num_time_intervals;
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);

#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	BOOST_LOG_TRIVIAL(info)
			<< boost::format("It takes %1% seconds") % time_diff(start, end);

	FG_vector<float>::ptr vec = FG_vector<float>::create(graph);
	graph->query_on_all(vertex_query::ptr(new save_query<float, scan_vertex>(vec)));
	return vec;
}

}
