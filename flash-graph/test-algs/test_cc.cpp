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

#include "stat.h"
#include "FGlib.h"
#include "graph.h"

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
			"cc [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-s size: the output min component size\n");
	fprintf(stderr, "-o file: output the component size to the file\n");
	fprintf(stderr, "-t: the type of connected components\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	std::string output_file;
	size_t min_comp_size = 0;
	int num_opts = 0;
	std::string type = "wcc";
	while ((opt = getopt(argc, argv, "c:s:o:t:")) != -1) {
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
			case 't':
				type = optarg;
				num_opts++;
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

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);
	FG_vector<vertex_id_t>::ptr comp_ids;
	if (type == "wcc")
		comp_ids = compute_wcc(graph);
	else if (type == "scc")
		comp_ids = compute_scc(graph);
	else
		assert(0);

	count_map<vertex_id_t> map;
	comp_ids->count_unique(map);
	std::pair<vertex_id_t, size_t> max_comp = map.get_max_count();
	printf("There are %ld components, and largest comp has %ld vertices\n",
			map.get_size(), max_comp.second);

	struct large_comp_apply
	{
		std::set<vertex_id_t> &large_comps;
		vsize_t threshold;
		large_comp_apply(std::set<vertex_id_t> &_comps,
				vsize_t threshold): large_comps(_comps) {
			this->threshold = threshold;
		}

		void operator()(const std::pair<vertex_id_t, size_t> &v) {
			if (v.second >= threshold)
				large_comps.insert(v.first);
		}
	};

	std::set<vertex_id_t> wanted_comps;
	map.apply(large_comp_apply(wanted_comps, min_comp_size));
	// We ignore the largest component.
	wanted_comps.erase(max_comp.first);
	printf("There are %ld components of the size larger than %ld\n",
			wanted_comps.size() + 1, min_comp_size);

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			assert(0);
		}

		std::map<vertex_id_t, graph::ptr> clusters;
		fetch_subgraphs(graph, comp_ids, wanted_comps, clusters);

		typedef std::map<vertex_id_t, std::pair<size_t, size_t> > size_map_t;
		size_map_t cluster_sizes;
		compute_subgraph_sizes(graph, comp_ids, wanted_comps, cluster_sizes);

		BOOST_FOREACH(size_map_t::value_type v, cluster_sizes) {
			std::map<vertex_id_t, graph::ptr>::const_iterator it
				= clusters.find(v.first);
			assert(it != clusters.end());
			assert(v.second.first == it->second->get_num_vertices());
			assert(v.second.second == it->second->get_num_edges());
			fprintf(f, "comp %u: %ld, %ld\n", v.first, v.second.first,
					v.second.second);
		}

		fclose(f);

		typedef std::pair<vertex_id_t, directed_graph<>::ptr> id_graph_t;
		BOOST_FOREACH(id_graph_t v, clusters) {
			v.second->dump_as_edge_list(output_file
					+ "-" + ltoa(v.first));
		}
	}
}
