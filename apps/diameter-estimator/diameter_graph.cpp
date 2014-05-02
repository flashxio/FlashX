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

#include <vector>
#include <unordered_map>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"

const int K = sizeof(uint64_t) * 8;
edge_type traverse_edge = edge_type::OUT_EDGE;

template<class T>
class embedded_bitmap
{
	T map;
public:
	embedded_bitmap() {
		map = 0;
	}

	void set(int idx) {
		assert((size_t) idx < sizeof(T) * 8);
		map |= ((T) 1) << idx;
	}

	void clear() {
		map = 0;
	}

	void merge(const embedded_bitmap<T> &map) {
		this->map |= map.map;
	}

	bool operator!=(const embedded_bitmap<T> &map) {
		return this->map != map.map;
	}

	T get_value() const {
		return map;
	}
};

class diameter_message: public vertex_message
{
	embedded_bitmap<uint64_t> bfs_ids;
public:
	diameter_message(embedded_bitmap<uint64_t> bfs_ids): vertex_message(
			sizeof(diameter_message), true) {
		this->bfs_ids = bfs_ids;
	}

	embedded_bitmap<uint64_t> get_bfs_ids() const {
		return bfs_ids;
	}
};

class diameter_vertex: public compute_directed_vertex
{
	embedded_bitmap<uint64_t> bfs_ids;
	embedded_bitmap<uint64_t> new_bfs_ids;
	// The largest distance from a start vertex among the BFS.
	short max_dist;
	bool updated;
public:
	diameter_vertex() {
		updated = false;
		max_dist = 0;
	}

	diameter_vertex(vertex_id_t id,
			const vertex_index &index): compute_directed_vertex(id, index) {
		updated = false;
		max_dist = 0;
	}

	short get_max_dist() const {
		return max_dist;
	}

	void init(int bfs_id) {
		max_dist = 0;
		bfs_ids.clear();
		bfs_ids.set(bfs_id);
		new_bfs_ids = bfs_ids;
		updated = true;
		printf("BFS %d starts at v%u\n", bfs_id, get_id());
	}

	void reset() {
		max_dist = 0;
		updated = false;
		bfs_ids.clear();
		new_bfs_ids.clear();
	}

	void run(vertex_program &prog) {
		if (updated) {
			updated = false;
			directed_vertex_request req(get_id(), traverse_edge);
			request_partial_vertices(&req, 1);
		}
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &vprog, const vertex_message &msg) {
		const diameter_message &dmsg = (const diameter_message &) msg;
		new_bfs_ids.merge(dmsg.get_bfs_ids());
		vprog.request_notify_iter_end(*this);
	}

	void notify_iteration_end(vertex_program &vprog);
};

class diameter_vertex_program: public vertex_program_impl<diameter_vertex>
{
	short max_dist;
public:
	typedef std::shared_ptr<diameter_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<diameter_vertex_program, vertex_program>(prog);
	}

	diameter_vertex_program() {
		max_dist = 0;
	}

	void set_max_dist(short max_dist) {
		if (this->max_dist < max_dist)
			this->max_dist = max_dist;
	}

	size_t get_max_dist() const {
		return max_dist;
	}
};

class diameter_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new diameter_vertex_program());
	}
};

void diameter_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(traverse_edge);
	if (num_dests == 0)
		return;

	// We need to add the neighbors of the vertex to the queue of
	// the next level.
	edge_seq_iterator it = vertex.get_neigh_seq_it(traverse_edge, 0, num_dests);
	diameter_message msg(bfs_ids);
	prog.multicast_msg(it, msg);
}

void diameter_vertex::notify_iteration_end(vertex_program &vprog)
{
	if (bfs_ids != new_bfs_ids) {
		updated = true;
		max_dist = max(vprog.get_graph().get_curr_level() + 1, max_dist);
		((diameter_vertex_program &) vprog).set_max_dist(max_dist);
	}
	bfs_ids = new_bfs_ids;
}

class diameter_initiator: public vertex_initiator
{
	std::unordered_map<vertex_id_t, int> start_vertices;
public:
	diameter_initiator(std::vector<vertex_id_t> &vertices) {
		for (size_t i = 0; i < vertices.size(); i++) {
			start_vertices.insert(std::pair<vertex_id_t, int>(vertices[i], i));
		}
	}

	void init(compute_vertex &v) {
		diameter_vertex &dv = (diameter_vertex &) v;
		std::unordered_map<vertex_id_t, int>::const_iterator it
			= start_vertices.find(v.get_id());
		assert(it != start_vertices.end());
		dv.init(it->second);
	}
};

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"diameter_estimator [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-n: the number of BFS\n");
	fprintf(stderr, "-o: the output file\n");
	fprintf(stderr, "-b: traverse with both in-edges and out-edges\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	int num_bfs = 1;
	bool preload = false;
	std::string output_file;
	while ((opt = getopt(argc, argv, "c:pn:o:b")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'n':
				num_bfs = atoi(optarg);
				num_bfs = min(num_bfs, K);
				num_opts++;
				break;
			case 'b':
				traverse_edge = edge_type::BOTH_EDGES;
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

	signal(SIGINT, int_handler);

	graph_index::ptr index = NUMA_graph_index<diameter_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("BFS starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	std::vector<vertex_id_t> start_vertices;
	while (start_vertices.size() < (size_t) num_bfs) {
		vertex_id_t id = random() % graph->get_max_vertex_id();
		// We should skip the empty vertices.
		if (graph->get_vertex_edges(id) == 0)
			continue;

		start_vertices.push_back(id);
	}

	short max_dist = 0;
	struct timeval start, end;
	gettimeofday(&start, NULL);
	graph->start(start_vertices.data(), start_vertices.size(),
			vertex_initiator::ptr(new diameter_initiator(start_vertices)),
			vertex_program_creater::ptr(new diameter_vertex_program_creater()));
	graph->wait4complete();
	gettimeofday(&end, NULL);

	std::vector<vertex_program::ptr> vprogs;
	graph->get_vertex_programs(vprogs);
	BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
		diameter_vertex_program::ptr diameter_vprog = diameter_vertex_program::cast2(vprog);
		max_dist = max(max_dist, diameter_vprog->get_max_dist());
	}
	printf("It takes %f seconds\n", time_diff(start, end));

#if 0
	for (size_t i = 0; i < start_vertices.size(); i++) {
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			diameter_vertex &v = (diameter_vertex &) *it;
			if (v.get_dist(i) != USHRT_MAX)
				fprintf(stderr, "v%u: %d\n", v.get_id(), v.get_dist(i));
		}
	}
#endif

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	printf("The max dist: %d\n", max_dist);

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		assert(f);
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			diameter_vertex &v = (diameter_vertex &) *it;
			fprintf(f, "%d\n", v.get_max_dist());
		}
		fclose(f);
	}
}
