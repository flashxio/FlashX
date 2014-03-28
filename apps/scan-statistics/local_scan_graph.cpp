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

#include <signal.h>
#include <google/profiler.h>

#include "scan_graph.h"

class count_msg: public vertex_message
{
	size_t num;
public:
	count_msg(size_t num): vertex_message(sizeof(count_msg), false) {
		this->num = num;
	}

	size_t get_num() const {
		return num;
	}
};

class extended_neighbor_list: public neighbor_list
{
	std::vector<uint32_t> count_list;
public:
	extended_neighbor_list(const page_vertex &vertex,
			const std::vector<attributed_neighbor> &neighbors): neighbor_list(
				vertex, neighbors) {
		count_list.resize(this->size());
	}

	uint32_t get_count(size_t idx) const {
		return count_list[idx];
	}

	size_t count_edges(const page_vertex *v);

	off_t find_idx(vertex_id_t id) {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return -1;
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return idx;
		}
	}

	attributed_neighbor find(vertex_id_t id) {
		edge_set_t::const_iterator it = neighbor_set->find(id);
		if (it == neighbor_set->end())
			return attributed_neighbor();
		else {
			size_t idx = (*it).get_idx();
			assert(idx < id_list.size());
			return at(idx);
		}
	}

	attributed_neighbor at(size_t idx) {
		return attributed_neighbor(id_list[idx], num_dup_list[idx]);
	}
};

size_t extended_neighbor_list::count_edges(const page_vertex *v)
{
	assert(!this->empty());
	if (v->get_num_edges(edge_type::BOTH_EDGES) == 0)
		return 0;

	std::vector<vertex_id_t> common_neighs1;
	std::vector<vertex_id_t> common_neighs2;
	size_t ret = neighbor_list::count_edges(v, edge_type::IN_EDGE, &common_neighs1)
		+ neighbor_list::count_edges(v, edge_type::OUT_EDGE, &common_neighs2);

	class skip_self {
	public:
		bool operator()(vertex_id_t id) {
			return false;
		}
	};

	class merge_edge {
	public:
		vertex_id_t operator()(vertex_id_t id1, vertex_id_t id2) {
			assert(id1 == id2);
			return id1;
		}
	};

	std::vector<vertex_id_t> common_neighs(common_neighs1.size()
			+ common_neighs2.size());
	size_t num_neighbors = unique_merge(
			common_neighs1.begin(), common_neighs1.end(),
			common_neighs2.begin(), common_neighs2.end(),
			skip_self(), merge_edge(), common_neighs.begin());
	common_neighs.resize(num_neighbors);

#ifdef PV_STAT
	rand_jumps += common_neighs.size() + 1;
#endif
	// The number of duplicated edges between v and this vertex.
	off_t neigh_off = this->find_idx(v->get_id());
	assert(neigh_off >= 0);
	attributed_neighbor neigh = this->at(neigh_off);
	int num_v_dups = neigh.get_num_dups();
	assert(num_v_dups > 0);
	size_t num_edges = 0;
	for (std::vector<vertex_id_t>::const_iterator it = common_neighs.begin();
			it != common_neighs.end(); it++) {
		off_t n_off = this->find_idx(*it);
		assert(n_off >= 0);
		attributed_neighbor n = this->at(n_off);
		num_edges += n.get_num_dups();
		count_list[n_off] += num_v_dups;
	}
	if (num_edges > 0)
		count_list[neigh_off] += num_edges;
	return ret;
}

class local_scan_vertex: public scan_vertex
{
public:
	local_scan_vertex() {
	}

	local_scan_vertex(vertex_id_t id,
			const vertex_index *index): scan_vertex(id, index) {
	}

	using scan_vertex::run;

	void run(graph_engine &graph) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	virtual void finding_triangles_end(graph_engine &graph) {
		extended_neighbor_list *neighbors
			= (extended_neighbor_list *) data->neighbors.get();
		// Inform all neighbors in the in-edges.
		for (size_t i = 0; i < data->neighbors->size(); i++) {
			size_t count = neighbors->get_count(i);
			if (count > 0) {
				count_msg msg(count);
				graph.send_msg(data->neighbors->get_neighbor_id(i), msg);
			}
		}
	}

	void run_on_messages(graph_engine &graph,
			const vertex_message *msgs[], int num) {
		for (int i = 0; i < num; i++) {
			const count_msg *msg = (const count_msg *) msgs[i];
			num_edges.inc(msg->get_num());
		}
	}

	virtual runtime_data_t *create_runtime(graph_engine &graph,
			const page_vertex *vertex);
};

class skip_larger {
	size_t size;
	vertex_id_t id;
	graph_engine &graph;
public:
	skip_larger(graph_engine &_graph, vertex_id_t id): graph(_graph) {
		this->size = graph.get_vertex_info(id).get_ext_mem_size();
		this->id = id;
	}

	bool operator()(attributed_neighbor &e) {
		return operator()(e.get_id());
	}

	/**
	 * We are going to count edges on the vertices with the most edges.
	 * If two vertices have the same number of edges, we compute
	 * on the vertices with the largest Id.
	 */
	bool operator()(vertex_id_t id) {
		const in_mem_vertex_info &info = graph.get_vertex_info(id);
		if (info.get_ext_mem_size() == size)
			return id >= this->id;
		return info.get_ext_mem_size() > size;
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

runtime_data_t *local_scan_vertex::create_runtime(graph_engine &graph,
		const page_vertex *vertex)
{
	merge_edge merge;
	std::vector<attributed_neighbor> neighbors(
			vertex->get_num_edges(edge_type::BOTH_EDGES));
	size_t num_neighbors = unique_merge(
			vertex->get_neigh_begin(edge_type::IN_EDGE),
			vertex->get_neigh_end(edge_type::IN_EDGE),
			vertex->get_neigh_begin(edge_type::OUT_EDGE),
			vertex->get_neigh_end(edge_type::OUT_EDGE),
			skip_larger(graph, vertex->get_id()), merge,
			neighbors.begin());
	neighbors.resize(num_neighbors);
	return new runtime_data_t(std::unique_ptr<neighbor_list>(
				new extended_neighbor_list(*vertex, neighbors)));
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"local-scan [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-o file: output local scan of each vertex to a file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string output_file;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "o:c:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 'c':
				confs = optarg;
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

	graph_index *index = NUMA_graph_index<local_scan_vertex>::create(
			index_file, graph_conf.get_num_threads(), params.get_num_nodes());
	graph_engine *graph = graph_engine::create(
			graph_conf.get_num_threads(), params.get_num_nodes(), graph_file,
			index);
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	struct timeval start, end;
	gettimeofday(&start, NULL);
	printf("Computing local scan\n");
	gettimeofday(&start, NULL);
	graph->start_all();
	graph->wait4complete();
	gettimeofday(&end, NULL);
	printf("It takes %f seconds\n", time_diff(start, end));
	printf("process %ld vertices and complete %ld vertices\n",
			num_working_vertices.get(), num_completed_vertices.get());

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
	graph_engine::destroy(graph);
	destroy_io_system();

#ifdef PV_STAT
	graph_index::const_iterator it = index->begin();
	graph_index::const_iterator end_it = index->end();
	size_t tot_scan_bytes = 0;
	size_t tot_rand_jumps = 0;
	for (; it != end_it; ++it) {
		const local_scan_vertex &v = (const local_scan_vertex &) *it;
		tot_scan_bytes += v.get_scan_bytes();
		tot_rand_jumps += v.get_rand_jumps();
	}
	printf("scan %ld bytes, %ld rand jumps\n",
			tot_scan_bytes, tot_rand_jumps);
#endif

	size_t max_scan = 0;
	vertex_id_t max_v = -1;
	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			return -1;
		}
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			const local_scan_vertex &v = (const local_scan_vertex &) *it;
			fprintf(f, "\"%ld\" %ld\n", (unsigned long) v.get_id(), v.get_local_scan());
			if (max_scan < v.get_local_scan()) {
				max_scan = v.get_local_scan();
				max_v = v.get_id();
			}
		}
		fclose(f);
	}

	if (max_scan == 0) {
		graph_index::const_iterator it = index->begin();
		graph_index::const_iterator end_it = index->end();
		for (; it != end_it; ++it) {
			const local_scan_vertex &v = (const local_scan_vertex &) *it;
			if (max_scan < v.get_local_scan()) {
				max_scan = v.get_local_scan();
				max_v = v.get_id();
			}
		}
	}
	printf("max scan: %ld, on v%u\n", max_scan, max_v);
}
