/*
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

#include "FGlib.h"
#include "ts_graph.h"
#include "matrix/FG_sparse_matrix.h"

using namespace fg;

void print_usage();

void int_handler(int sig_num)
{
#ifdef PROFILER
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
#endif
	exit(0);
}

void run_cycle_triangle(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	bool fast = false;

	while ((opt = getopt(argc, argv, "f")) != -1) {
		num_opts++;
		switch (opt) {
			case 'f':
				fast = true;
				break;
			default:
				print_usage();
				abort();
		}
	}

	FG_vector<size_t>::ptr triangles;
	if (fast)
		triangles = compute_directed_triangles_fast(graph,
				directed_triangle_type::CYCLE);
	else
		triangles = compute_directed_triangles(graph,
				directed_triangle_type::CYCLE);
	printf("There are %ld cycle triangles\n", triangles->sum());
}

void run_triangle(FG_graph::ptr graph, int argc, char *argv[])
{
	FG_vector<size_t>::ptr triangles;
	triangles = compute_undirected_triangles(graph);
	printf("There are %ld triangles\n", triangles->sum());
}

void run_local_scan(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	int num_hops = 1;

	while ((opt = getopt(argc, argv, "H:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'H':
				num_hops = atoi(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	FG_vector<size_t>::ptr scan;
	if (num_hops == 1)
		scan = compute_local_scan(graph);
	else if (num_hops == 2)
		scan = compute_local_scan2(graph);
	else {
		fprintf(stderr, "we don't support local scan of more than 2 hops\n");
		exit(1);
	}
	std::pair<size_t, off_t> ret = scan->max_val_loc();
	printf("Max local scan is %ld on v%ld\n", ret.first, ret.second);
}

void run_topK_scan(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	int topK = 1;

	while ((opt = getopt(argc, argv, "K:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'K':
				topK = atoi(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	FG_vector<std::pair<vertex_id_t, size_t> >::ptr scan
		= compute_topK_scan(graph, topK);
	printf("The top %d scans:\n", topK);
	for (int i = 0; i < topK; i++)
		printf("%u\t%ld\n", scan->get(i).first, scan->get(i).second);
}

void print_cc(FG_vector<vertex_id_t>::ptr comp_ids)
{
	count_map<vertex_id_t> map;
	comp_ids->count_unique(map);
	if (map.exists(INVALID_VERTEX_ID)) {
		printf("There are %ld empty vertices\n",
				map.get(INVALID_VERTEX_ID));
		map.remove(INVALID_VERTEX_ID);
	}
	std::pair<vertex_id_t, size_t> max_comp = map.get_max_count();
	printf("There are %ld components (exclude empty vertices), and largest comp has %ld vertices\n",
			map.get_size(), max_comp.second);
}

void run_cc(FG_graph::ptr graph, int argc, char *argv[])
{
	print_cc(compute_cc(graph));
}

void run_wcc(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	bool sync = false;
	std::string output_file;
	while ((opt = getopt(argc, argv, "so:")) != -1) {
		num_opts++;
		switch (opt) {
			case 's':
				sync = true;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}
	FG_vector<vertex_id_t>::ptr comp_ids;
	if (sync)
		comp_ids = compute_sync_wcc(graph);
	else
		comp_ids = compute_wcc(graph);
	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			return;
		}
		for (size_t i = 0; i < comp_ids->get_size(); i++)
			fprintf(f, "%ld %d\n", i, comp_ids->get(i));
		fclose(f);
	}
	print_cc(comp_ids);
}

void run_scc(FG_graph::ptr graph, int argc, char *argv[])
{
	print_cc(compute_scc(graph));
}

void run_diameter(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;

	int num_para_bfs = 1;
	bool directed = false;

	while ((opt = getopt(argc, argv, "p:ds:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'p':
				num_para_bfs = atoi(optarg);
				num_opts++;
				break;
			case 'd':
				directed = true;
				break;
			default:
				print_usage();
				abort();
		}
	}

	size_t diameter = estimate_diameter(graph, num_para_bfs, directed);
	printf("The estimated diameter is %ld\n", diameter);
}

void run_pagerank(FG_graph::ptr graph, int argc, char *argv[], int version)
{
	int opt;
	int num_opts = 0;

	int num_iters = 30;
	float damping_factor = 0.85;

	while ((opt = getopt(argc, argv, "i:D:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'i':
				num_iters = atoi(optarg);
				num_opts++;
				break;
			case 'D':
				damping_factor = atof(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	FG_vector<float>::ptr pr;
	switch (version) {
		case 1:
			pr = compute_pagerank(graph, num_iters, damping_factor);
			break;
		case 2:
			pr = compute_pagerank2(graph, num_iters, damping_factor);
			break;
		default:
			abort();
	}
	std::vector<std::pair<float, off_t> > val_locs;
	pr->max_val_locs(10, val_locs);
	printf("The sum of pagerank of all vertices: %f\n", pr->sum());
	for (size_t i = 0; i < val_locs.size(); i++)
		printf("v%ld: %f\n", val_locs[i].second, val_locs[i].first);
}

void run_sstsg(FG_graph::ptr graph, int argc, char *argv[])
{
	std::string start_time_str;
	std::string time_unit_str;
	std::string output_file;
	int num_time_intervals = 1;
	long time_interval = 1;
	bool compute_all = false;
	time_t start_time = -1;

	int opt;
	int num_opts = 0;

	while ((opt = getopt(argc, argv, "n:u:o:t:l:a")) != -1) {
		num_opts++;
		switch (opt) {
			case 'n':
				num_time_intervals = atoi(optarg);
				num_opts++;
				break;
			case 'u':
				time_unit_str = optarg;
				num_opts++;
				break;
			case 'o':
				output_file = optarg;
				num_opts++;
				break;
			case 't':
				start_time_str = optarg;
				if (is_time_str(start_time_str))
					start_time = conv_str_to_time(start_time_str);
				else
					start_time = atol(start_time_str.c_str());
				num_opts++;
				break;
			case 'l':
				time_interval = atol(optarg);
				num_opts++;
				break;
			case 'a':
				compute_all = true;
				break;
			default:
				print_usage();
				abort();
		}
	}

	if (time_unit_str == "hour")
		time_interval *= HOUR_SECS;
	else if (time_unit_str == "day")
		time_interval *= DAY_SECS;
	else if (time_unit_str == "month")
		time_interval *= MONTH_SECS;
	else
		fprintf(stderr, "a wrong time unit: %s\n", time_unit_str.c_str());

#if 0
	namespace bt = boost::posix_time;
	printf("start time: %s\n", start_time_str.c_str());
	bt::ptime pt = bt::time_from_string(start_time_str);
	struct tm tm = bt::to_tm(pt);
#endif
	if (compute_all) {
		std::pair<time_t, time_t> range = get_time_range(graph);
		std::string start_time_str = ctime(&range.first);
		start_time_str[start_time_str.length() - 1] = 0;
		std::string end_time_str = ctime(&range.second);
		end_time_str[end_time_str.length() - 1] = 0;
		printf("the time-series graph starts at %s, ends at %s\n",
				start_time_str.c_str(), end_time_str.c_str());
		time_t start_time = range.first;
		time_t end_time = range.second;
		for (time_t interval_start
				= start_time + num_time_intervals * time_interval;
				interval_start < end_time; interval_start += time_interval) {
			FG_vector<float>::ptr res = compute_sstsg(graph, interval_start,
					time_interval, num_time_intervals);
			std::pair<float, off_t> p = res->max_val_loc();
			printf("v%ld has max scan %f\n", p.second, p.first);
		}
	}
	else {
		printf("start time: %ld, interval: %ld\n", start_time, time_interval);
		FG_vector<float>::ptr res = compute_sstsg(graph, start_time,
				time_interval, num_time_intervals);

		std::pair<float, off_t> p = res->max_val_loc();
		printf("v%ld has max scan %f\n", p.second, p.first);
		if (!output_file.empty()) {
			FILE *f = fopen(output_file.c_str(), "w");
			if (f == NULL) {
				perror("fopen");
				return;
			}
			for (size_t i = 0; i < res->get_size(); i++)
				fprintf(f, "\"%ld\" %f\n", i, res->get(i));
			fclose(f);
		}
	}
}

void run_ts_wcc(FG_graph::ptr graph, int argc, char *argv[])
{
	std::string start_time_str;
	std::string time_unit_str;
	long time_interval = 1;

	int opt;
	int num_opts = 0;

	while ((opt = getopt(argc, argv, "u:t:l:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'u':
				time_unit_str = optarg;
				num_opts++;
				break;
			case 't':
				start_time_str = optarg;
				num_opts++;
				break;
			case 'l':
				time_interval = atol(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	if (time_unit_str == "hour")
		time_interval *= HOUR_SECS;
	else if (time_unit_str == "day")
		time_interval *= DAY_SECS;
	else if (time_unit_str == "month")
		time_interval *= MONTH_SECS;
	else
		fprintf(stderr, "a wrong time unit: %s\n", time_unit_str.c_str());

#if 0
	namespace bt = boost::posix_time;
	printf("start time: %s\n", start_time_str.c_str());
	bt::ptime pt = bt::time_from_string(start_time_str);
	struct tm tm = bt::to_tm(pt);
#endif
	time_t start_time = conv_str_to_time(start_time_str);
	printf("start time: %ld, interval: %ld\n", start_time, time_interval);
	FG_vector<vertex_id_t>::ptr comp_ids = compute_ts_wcc(graph, start_time,
			time_interval);
	print_cc(comp_ids);
}

void run_kcore(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	int num_opts = 0;
	size_t kmax = 0;
	size_t k = 0;
	std::string write_out = "";

	while ((opt = getopt(argc, argv, "k:m:w:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'k':
				k = atol(optarg);
				num_opts++;
				break;
			case 'm':
				kmax = atol(optarg);
				num_opts++;
				break;
			case 'w':
				write_out = optarg;
				break;
			default:
				print_usage();
				abort();
		}
	}

	if (k < 2) {
		fprintf(stderr, "[Error]: kmin cannot be < 2\n");
		exit(-1);
	}

	FG_vector<size_t>::ptr kcorev = compute_kcore(graph, k, kmax);
	if (!write_out.empty())
		kcorev->to_file(write_out);
}

int read_vertices(const std::string &file, std::vector<vertex_id_t> &vertices)
{
	FILE *f = fopen(file.c_str(), "r");
	assert(f);
	ssize_t ret;
	char *line = NULL;
	size_t line_size = 0;
	while ((ret = getline(&line, &line_size, f)) > 0) {
		if (line[ret - 1] == '\n')
			line[ret - 1] = 0;
		vertex_id_t id = atol(line);
		printf("%u\n", id);
		vertices.push_back(id);
	}
	fclose(f);
	return vertices.size();
}

void run_overlap(FG_graph::ptr graph, int argc, char* argv[])
{
	std::string output_file;
	std::string confs;

	if (argc < 2) {
		fprintf(stderr, "overlap requires vertex_file\n");
		exit(-1);
	}
	std::string vertex_file = argv[1];

	int opt;
	int num_opts = 0;
	std::string write_out;
	double threshold = 0;
	while ((opt = getopt(argc, argv, "o:t:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'o':
				write_out = optarg;
				num_opts++;
				break;
			case 't':
				threshold = atof(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	std::vector<vertex_id_t> overlap_vertices;
	read_vertices(vertex_file, overlap_vertices);
	std::vector<std::vector<double> > overlaps;
	std::sort(overlap_vertices.begin(), overlap_vertices.end());
	compute_overlap(graph, overlap_vertices, overlaps);

	if (!write_out.empty()) {
		FILE *fout = fopen(write_out.c_str(), "w");
		assert(fout);
		size_t num_vertices = overlap_vertices.size();
		assert(num_vertices == overlaps.size());
		for (size_t i = 0; i < num_vertices; i++) {
			assert(num_vertices == overlaps[i].size());
			for (size_t j = 0; j < num_vertices; j++) {
				double overlap = overlaps[i][j];
				if (overlap >= threshold)
					fprintf(fout, "%u %u %f\n", overlap_vertices[i],
							overlap_vertices[j], overlap);
			}
		}
		fclose(fout);
	}
}

void run_bfs(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	int num_opts = 0;
	edge_type edge = edge_type::OUT_EDGE;
	vertex_id_t start_vertex = 0;

	std::string edge_type_str;
	while ((opt = getopt(argc, argv, "e:s:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'e':
				edge_type_str = optarg;
				num_opts++;
				break;
			case 's':
				start_vertex = atol(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}
	if (!edge_type_str.empty()) {
		if (edge_type_str == "IN")
			edge = edge_type::IN_EDGE;
		else if (edge_type_str == "OUT")
			edge = edge_type::OUT_EDGE;
		else if (edge_type_str == "BOTH")
			edge = edge_type::BOTH_EDGES;
		else {
			fprintf(stderr, "wrong edge type");
			exit(1);
		}
	}

	size_t bfs(FG_graph::ptr fg, vertex_id_t start_vertex, edge_type);
	size_t num_vertices = bfs(graph, start_vertex, edge);
	printf("BFS from v%u traverses %ld vertices on edge type %d\n",
			start_vertex, num_vertices, edge);
}

void run_spmv(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	FG_adj_matrix::ptr m = FG_adj_matrix::create(graph);
	while ((opt = getopt(argc, argv, "t")) != -1) {
		switch (opt) {
			case 't':
				m = m->transpose();
				break;
			default:
				print_usage();
				abort();
		}
	}

	FG_vector<double>::ptr in = FG_vector<double>::create(m->get_num_cols());
	in->init(1);
	FG_vector<double>::ptr out = FG_vector<double>::create(m->get_num_rows());
	m->multiply<double>(*in, *out);
	printf("the graph has %ld vertices and %ld edges\n",
			graph->get_graph_header().get_num_vertices(),
			graph->get_graph_header().get_num_edges());
	printf("SpMV returns %lf\n", out->sum());
}

std::string supported_algs[] = {
	"cycle_triangle",
	"triangle",
	"local_scan",
	"topK_scan",
	"wcc",
	"scc",
	"diameter",
	"pagerank",
	"pagerank2",
	"sstsg",
	"ts_wcc",
	"kcore",
	"overlap",
	"bfs",
	"spmv",
};
int num_supported = sizeof(supported_algs) / sizeof(supported_algs[0]);

void print_usage()
{
	fprintf(stderr,
			"test_algs conf_file graph_file index_file algorithm [alg-options]\n");
	fprintf(stderr, "scan-statistics:\n");
	fprintf(stderr, "-K topK: topK vertices in topK scan\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "local scan\n");
	fprintf(stderr, "-H hops: local scan within the specified number of hops\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "diameter estimation:\n");
	fprintf(stderr, "-p num_para_bfs: the number of parallel bfs to estimate diameter\n");
	fprintf(stderr, "-d: whether we respect the direction of edges\n");
	fprintf(stderr, "-s num: the number of sweeps performed in diameter estimation\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "pagerank\n");
	fprintf(stderr, "-i num: the maximum number of iterations\n");
	fprintf(stderr, "-D v: damping factor\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "sstsg\n");
	fprintf(stderr, "-n num: the number of time intervals\n");
	fprintf(stderr, "-u unit: time unit (hour, day, month, etc)\n");
	fprintf(stderr, "-o output: the output file\n");
	fprintf(stderr, "-t time: the start time\n");
	fprintf(stderr, "-l time: the length of time interval\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "ts_wcc\n");
	fprintf(stderr, "-u unit: time unit (hour, day, month, etc)\n");
	fprintf(stderr, "-t time: the start time\n");
	fprintf(stderr, "-l time: the length of time interval\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "k-core:\n");
	fprintf(stderr, "-k k: the minimum k value to compute\n");
	fprintf(stderr, "-m kmax: the maximum k value to compute\n");
	fprintf(stderr, "-w output: the file name for a vector written to filen");
	fprintf(stderr, "\n");
	fprintf(stderr, "cycle_triangle\n");
	fprintf(stderr, "-f: run the fast implementation\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "wcc\n");
	fprintf(stderr, "-s: run wcc synchronously\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "overlap vertex_file\n");
	fprintf(stderr, "-o output: the output file\n");
	fprintf(stderr, "-t threshold: the threshold for printing the overlaps\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "bfs\n");
	fprintf(stderr, "-e edge type: the type of edge to traverse (IN, OUT, BOTH)\n");
	fprintf(stderr, "-s vertex id: the vertex where the BFS starts\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "spmv\n");
	fprintf(stderr, "-t: transpose the sparse matrix.\n");

	fprintf(stderr, "supported graph algorithms:\n");
	for (int i = 0; i < num_supported; i++)
		fprintf(stderr, "\t%s\n", supported_algs[i].c_str());
	graph_conf.print_help();
	safs::params.print_help();
}

int main(int argc, char *argv[])
{
	argv++;
	argc--;
	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];
	std::string alg = argv[3];
	// We should increase by 3 instead of 4. getopt() ignores the first
	// argument in the list.
	argv += 3;
	argc -= 3;

	config_map::ptr configs = config_map::create(conf_file);
	if (configs == NULL)
		configs = config_map::ptr();
	signal(SIGINT, int_handler);

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);
	if (alg == "cycle_triangle") {
		run_cycle_triangle(graph, argc, argv);
	}
	else if (alg == "triangle") {
		run_triangle(graph, argc, argv);
	}
	else if (alg == "local_scan") {
		run_local_scan(graph, argc, argv);
	}
	else if (alg == "topK_scan") {
		run_topK_scan(graph, argc, argv);
	}
	else if (alg == "diameter") {
		run_diameter(graph, argc, argv);
	}
	else if (alg == "pagerank") {
		run_pagerank(graph, argc, argv, 1);
	}
	else if (alg == "pagerank2") {
		run_pagerank(graph, argc, argv, 2);
	}
	else if (alg == "wcc") {
		run_wcc(graph, argc, argv);
	}
	else if (alg == "cc") {
		run_cc(graph, argc, argv);
	}
	else if (alg == "scc") {
		run_scc(graph, argc, argv);
	}
	else if (alg == "sstsg") {
		run_sstsg(graph, argc, argv);
	}
	else if (alg == "ts_wcc") {
		run_ts_wcc(graph, argc, argv);
	}
	else if (alg == "kcore") {
		run_kcore(graph, argc, argv);
	}
	else if (alg == "overlap") {
		run_overlap(graph, argc, argv);
	}
	else if (alg == "bfs") {
		run_bfs(graph, argc, argv);
	}
	else if (alg == "spmv") {
		run_spmv(graph, argc, argv);
	}
}
