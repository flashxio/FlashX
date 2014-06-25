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
#include <google/profiler.h>
#endif

#include <set>
#include <vector>

#include "thread.h"
#include "io_interface.h"

#include "graph_engine.h"
#include "graph_config.h"

const double BIN_SEARCH_RATIO = 10;
const int HOUR_SECS = 3600;
const int DAY_SECS = HOUR_SECS * 24;
const int MONTH_SECS = DAY_SECS * 30;

atomic_number<long> num_working_vertices;
atomic_number<long> num_completed_vertices;

time_t timestamp;
time_t time_interval = 3600 * 24;
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

	size_t count_edges(const page_directed_vertex &v,
			const std::vector<vertex_id_t> *neighbors, time_t timestamp,
			time_t time_interval);
	size_t count_edges(const page_directed_vertex &v,
			const std::vector<vertex_id_t> *neighbors, time_t timestamp,
			time_t time_interval, edge_type type);

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex) {
		if (vertex.get_id() == get_id())
			run_on_itself(prog, (const page_directed_vertex &) vertex);
		else
			run_on_neighbor(prog, (const page_directed_vertex &) vertex);
	}

	void run_on_itself(vertex_program &prog, const page_directed_vertex &vertex);
	void run_on_neighbor(vertex_program &prog, const page_directed_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

page_byte_array::seq_const_iterator<vertex_id_t> get_ts_iterator(
		const page_directed_vertex &v, edge_type type, time_t time_start,
		time_t time_interval)
{
	page_byte_array::const_iterator<ts_edge_data> begin_it
		= v.get_data_begin<ts_edge_data>(type);
	page_byte_array::const_iterator<ts_edge_data> end_it
		= v.get_data_end<ts_edge_data>(type);
	page_byte_array::const_iterator<ts_edge_data> ts_it = std::lower_bound(
			begin_it, end_it, ts_edge_data(time_start));
	// All timestamps are smaller than time_start
	if (ts_it == end_it)
		return v.get_neigh_seq_it(type, 0, 0);
	size_t start = ts_it - begin_it;

	page_byte_array::const_iterator<ts_edge_data> ts_end_it = std::lower_bound(
			begin_it, end_it, time_start + time_interval);
	size_t end;
	if (ts_end_it == end_it
			|| (*ts_end_it).get_timestamp() >= time_start + time_interval)
		end = end_it - begin_it;
	else
		// The timestamp pointed by ts_end_t may be covered by the range.
		end = ts_end_it - begin_it + 1;

	return v.get_neigh_seq_it(type, start, end);
}

size_t scan_vertex::count_edges(const page_directed_vertex &v,
		const std::vector<vertex_id_t> *neighbors, time_t timestamp,
		time_t time_interval, edge_type type)
{
	size_t num_local_edges = 0;
	page_byte_array::seq_const_iterator<vertex_id_t> it = get_ts_iterator(
			v, type, timestamp, time_interval);
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
				&& neigh_neighbor != this->get_id()) {
			if (std::binary_search(this_it, this_end, neigh_neighbor))
				num_local_edges++;
		}
	} PAGE_FOREACH_END
	return num_local_edges;
}

size_t scan_vertex::count_edges(const page_directed_vertex &v,
		const std::vector<vertex_id_t> *neighbors, time_t timestamp,
		time_t time_interval)
{
	assert(!neighbors->empty());
	return count_edges(v, neighbors, timestamp, time_interval, edge_type::IN_EDGE)
		+ count_edges(v, neighbors, timestamp, time_interval, edge_type::OUT_EDGE);
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
	page_byte_array::seq_const_iterator<vertex_id_t> it = get_ts_iterator(
			v, type, time_start, time_interval);
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

	long ret = num_working_vertices.inc(1);
	if (ret % 100000 == 0)
		printf("%ld working vertices\n", ret);
	std::vector<vertex_id_t> in_neighbors;
	std::vector<vertex_id_t> out_neighbors;
	get_neighbors(vertex, edge_type::IN_EDGE, timestamp, time_interval,
			in_neighbors);
	get_neighbors(vertex, edge_type::OUT_EDGE, timestamp, time_interval,
			out_neighbors);
	if (in_neighbors.size() + out_neighbors.size() == 0) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

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

	int num_neighbors = unique_merge(
			in_neighbors.begin(), in_neighbors.end(),
			out_neighbors.begin(), out_neighbors.end(),
			skip_self(vertex.get_id()),
			neighbors->begin());
	neighbors->resize(num_neighbors);

	local_scans->at(0) = in_neighbors.size() + out_neighbors.size();
	// If there is a self-loop in the in-edge list
	if (std::binary_search(in_neighbors.begin(), in_neighbors.end(), get_id()))
		local_scans->at(0)--;
	// If there is a self-loop in the out-edge list
	if (std::binary_search(out_neighbors.begin(), out_neighbors.end(), get_id()))
		local_scans->at(0)--;

	if (neighbors->size() == 0) {
		delete local_scans;
		delete neighbors;
		local_scans = NULL;
		neighbors = NULL;
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);
		return;
	}

	// Count the degree of the vertex in other time intervals.
	for (int ts_idx = 1; ts_idx < num_time_intervals
			&& timestamp >= ts_idx * time_interval; ts_idx++) {
		time_t timestamp2 = timestamp - ts_idx * time_interval;

		// For in-edges.
		page_byte_array::seq_const_iterator<vertex_id_t> it
			= get_ts_iterator(vertex, edge_type::IN_EDGE, timestamp2,
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
		size_t ret = count_edges(vertex, neighbors,
				timestamp - j * time_interval, time_interval);
		if (ret > 0)
			local_scans->at(j) += ret;
	}

	// If we have seen all required neighbors, we have complete
	// the computation. We can release the memory now.
	if (num_joined == (int) neighbors->size()) {
		long ret = num_completed_vertices.inc(1);
		if (ret % 100000 == 0)
			printf("%ld completed vertices\n", ret);

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

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"scan-statistics conf_file graph_file index_file [output_file]\n");
	fprintf(stderr, "-c conf\n");
	fprintf(stderr, "-n num: the number of time intervals\n");
	fprintf(stderr, "-m: use month for time unit\n");
	fprintf(stderr, "-d: use day for time unit\n");
	fprintf(stderr, "-h: use hour for time unit\n");
	fprintf(stderr, "-o output: the output file\n");
	fprintf(stderr, "-t time: the start time\n");
	fprintf(stderr, "-i time: the length of time interval\n");
	graph_conf.print_help();
	params.print_help();
	exit(-1);
}

int main(int argc, char *argv[])
{
	printf("time_t size: %ld\n", sizeof(time_t));
	int opt;
	int num_opts = 0;
	std::string confs;
	std::string output_file;
	int time_unit = 1;
	while ((opt = getopt(argc, argv, "c:n:mdho:t:i:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'n':
				num_time_intervals = atoi(optarg);
				num_opts++;
				break;
			case 'm':
				time_unit = MONTH_SECS;
				break;
			case 'd':
				time_unit = DAY_SECS;
				break;
			case 'h':
				time_unit = HOUR_SECS;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 't':
				timestamp = atol(optarg);
				num_opts++;
				break;
			case 'i':
				time_interval = atol(optarg);
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	timestamp *= time_unit;
	time_interval *= time_unit;
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 3) {
		print_usage();
		exit(1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];

	config_map configs(conf_file);
	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<scan_vertex>::create(
			index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	assert(graph->get_graph_header().get_graph_type() == graph_type::DIRECTED);
	assert(graph->get_graph_header().has_edge_data());
	printf("scan statistics starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
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
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("It takes %f seconds\n", time_diff(start, end));
	printf("There are %ld vertices\n", index->get_num_vertices());
	printf("process %ld vertices and complete %ld vertices\n",
			num_working_vertices.get(), num_completed_vertices.get());

	double max_res = LONG_MIN;
	vertex_id_t max_v = -1;
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	for (; it != end_it; ++it) {
		const scan_vertex &v = (const scan_vertex &) *it;
		if (max_res < v.get_result()) {
			max_v = v.get_id();
			max_res = v.get_result();
		}
	}
	printf("max value is on v%ld: %f\n", (unsigned long) max_v, max_res);

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			return -1;
		}
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			const scan_vertex &v = (const scan_vertex &) *it;
			fprintf(f, "\"%ld\" %f\n", (unsigned long) v.get_id(), v.get_result());
		}
		fclose(f);
	}
}
