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

/**
 * This program computes local scan or # triangles on each vertex that belongs
 * to the largest weakly connected component.
 */

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "FGlib.h"

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
			"test_lcc [options] conf_file graph_file index_file algorithm\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "c:")) != -1) {
		num_opts++;
		switch (opt) {
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

	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	std::string alg = argv[3];

	config_map::ptr configs = config_map::create(conf_file);
	assert(configs);
	configs->add_options(confs);

	signal(SIGINT, int_handler);

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);
	FG_vector<size_t>::ptr counts;
	if (alg == "cycle_triangle")
		counts = compute_directed_triangles(graph, directed_triangle_type::CYCLE);
	else if (alg == "triangle")
		counts = compute_undirected_triangles(graph);
	else if (alg == "local_scan")
		counts = compute_local_scan(graph);

	FG_vector<vertex_id_t>::ptr comp_ids = compute_wcc(graph);
	count_map<vertex_id_t> comp_map;
	comp_ids->count_unique(comp_map);
	std::pair<vertex_id_t, size_t> lcc = comp_map.get_max_count();
	vertex_id_t lcc_id = lcc.first;
	size_t lcc_size = lcc.second;
	assert(comp_ids->get_size() == counts->get_size());

	count_map<size_t> map;
	for (size_t i = 0; i < comp_ids->get_size(); i++) {
		if (comp_ids->get(i) == lcc_id)
			map.add(counts->get(i));
	}

	static size_t tot_num = 0;

	class print
	{
	public:
		void operator()(const std::pair<size_t, size_t> &p) {
			tot_num += p.second;
			fprintf(stderr, "%ld %ld\n", p.first, p.second);
		}
	};
	map.apply(print());
	assert(tot_num == lcc_size);
	printf("There are %ld vertices in the lcc\n", tot_num);
}
