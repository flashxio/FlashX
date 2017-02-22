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
#include "sparse_matrix.h"
#include "libgraph-algs/sem_kmeans.h"

#include "vector.h"
#include "col_vec.h"
#include "data_frame.h"

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

	fm::vector::ptr triangles;
	if (fast)
		triangles = compute_directed_triangles_fast(graph,
				directed_triangle_type::CYCLE);
	else
		triangles = compute_directed_triangles(graph,
				directed_triangle_type::CYCLE);
	if (triangles)
		printf("There are %ld cycle triangles\n", triangles->sum<size_t>());
}

void run_triangle(FG_graph::ptr graph, int argc, char *argv[])
{
	fm::vector::ptr triangles;
	triangles = compute_undirected_triangles(graph);
	if (triangles)
		printf("There are %ld triangles\n", triangles->sum<size_t>());
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

	fm::vector::ptr scan;
	if (num_hops == 1)
		scan = compute_local_scan(graph);
	else if (num_hops == 2)
		scan = compute_local_scan2(graph);
	else {
		fprintf(stderr, "we don't support local scan of more than 2 hops\n");
		exit(1);
	}
	if (scan) {
//		std::pair<size_t, off_t> ret = scan->max_val_loc();
//		printf("Max local scan is %ld on v%ld\n", ret.first, ret.second);
	}
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
	if (scan) {
		printf("The top %d scans:\n", topK);
		for (int i = 0; i < topK; i++)
			printf("%u\t%ld\n", scan->get(i).first, scan->get(i).second);
	}
}

fm::data_frame::ptr get_cluster_counts(fm::vector::ptr comp_ids)
{
	fm::col_vec::ptr vec = fm::col_vec::create(comp_ids);
	fm::bulk_operate::const_ptr add = fm::bulk_operate::conv2ptr(
			fm::get_scalar_type<size_t>().get_basic_ops().get_add());
	fm::agg_operate::const_ptr sum = fm::agg_operate::create(add);;
	return vec->groupby(sum, true);
}

std::pair<fg::vertex_id_t, size_t> get_max_cid(fm::data_frame::ptr counts)
{
	fm::detail::mem_vec_store::const_ptr ids
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				counts->get_vec(0));
	fm::detail::mem_vec_store::const_ptr cnts
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				counts->get_vec(1));
	assert(cnts->get_type() == fm::get_scalar_type<size_t>());
	size_t idx = 0;
	size_t max_size = cnts->get<size_t>(0);
	for (size_t i = 1; i < cnts->get_length(); i++) {
		// We should skip the clusters with only one vertex.
		if (ids->get<vertex_id_t>(i) == INVALID_VERTEX_ID)
			continue;

		if (cnts->get<size_t>(i) > max_size) {
			max_size = cnts->get<size_t>(i);
			idx = i;
		}
	}
	fg::vertex_id_t id = ids->get<fg::vertex_id_t>(idx);
	return std::pair<fg::vertex_id_t, size_t>(id, max_size);
}

void print_cc(fm::vector::ptr comp_ids)
{
	fm::data_frame::ptr cnts = get_cluster_counts(comp_ids);
	std::pair<vertex_id_t, size_t> max_comp = get_max_cid(cnts);
	printf("There are %ld components (exclude empty vertices), and largest comp has %ld vertices\n",
			cnts->get_num_entries(), max_comp.second);
}

void run_cc(FG_graph::ptr graph, int argc, char *argv[])
{
	fm::vector::ptr cc = compute_cc(graph);
	if (cc)
		print_cc(cc);
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
	fm::vector::ptr comp_ids;
	if (sync)
		comp_ids = compute_sync_wcc(graph);
	else
		comp_ids = compute_wcc(graph);
	if (comp_ids == NULL)
		return;

	if (!output_file.empty()) {
		FILE *f = fopen(output_file.c_str(), "w");
		if (f == NULL) {
			perror("fopen");
			return;
		}
		fm::detail::mem_vec_store::const_ptr mem_ids
			= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
					comp_ids->get_raw_store());
		for (size_t i = 0; i < comp_ids->get_length(); i++)
			fprintf(f, "%ld %d\n", i, mem_ids->get<vertex_id_t>(i));
		fclose(f);
	}
	print_cc(comp_ids);
}

void run_scc(FG_graph::ptr graph, int argc, char *argv[])
{
	fm::vector::ptr cc = compute_scc(graph);
	if (cc)
		print_cc(cc);
}

void run_diameter(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;

	int num_para_bfs = 1;
	int num_sweeps = std::numeric_limits<int>::max();
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
			case 's':
				num_sweeps = atoi(optarg);
				num_sweeps = num_sweeps; // STUB
				fprintf(stderr, "[Warning]: num_sweeps argument currently unused\n");
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

	fm::vector::ptr pr;
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
	if (pr == NULL)
		return;

	printf("The sum of pagerank of all vertices: %f\n", pr->sum<float>());

	// Get the top N pagerank values and the corresponding vertices.
	fm::detail::mem_vec_store::const_ptr pr_store
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				pr->get_raw_store());
	typedef std::pair<float, off_t> val_loc_t;
	struct comp_val {
		bool operator()(const val_loc_t &v1, const val_loc_t &v2) {
			return v1.first > v2.first;
		}
	};
	std::priority_queue<val_loc_t, std::vector<val_loc_t>, comp_val> queue;
	for (size_t i = 0; i < pr->get_length(); i++) {
		float val = pr_store->get<float>(i);
		queue.push(val_loc_t(val, i));
		if (queue.size() > 10)
			queue.pop();
	}
	while (!queue.empty()) {
		val_loc_t pair = queue.top();
		printf("v%ld: %f\n", pair.second, pair.first);
		queue.pop();
	}
}

template<class T>
std::pair<T, off_t> max_val_loc(fm::vector::ptr res)
{
	fm::detail::mem_vec_store::const_ptr res_store
		= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
				res->get_raw_store());
	T ret = std::numeric_limits<T>::min();
	off_t idx = 0;
	// TODO we should be able to parallelize it.
	for (size_t i = 0; i < res->get_length(); i++) {
		if (ret < res_store->get<T>(i)) {
			ret = res_store->get<T>(i);
			idx = i;
		}
	}
	return std::pair<T, off_t>(ret, idx);
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
			fm::vector::ptr res = compute_sstsg(graph, interval_start,
					time_interval, num_time_intervals);
			std::pair<float, off_t> p = max_val_loc<float>(res);
			printf("v%ld has max scan %f\n", p.second, p.first);
		}
	}
	else {
		printf("start time: %ld, interval: %ld\n", start_time, time_interval);
		fm::vector::ptr res = compute_sstsg(graph, start_time,
				time_interval, num_time_intervals);

		std::pair<float, off_t> p = max_val_loc<float>(res);
		printf("v%ld has max scan %f\n", p.second, p.first);
		if (!output_file.empty()) {
			FILE *f = fopen(output_file.c_str(), "w");
			if (f == NULL) {
				perror("fopen");
				return;
			}
			fm::detail::mem_vec_store::const_ptr res_store
				= std::dynamic_pointer_cast<const fm::detail::mem_vec_store>(
						res->get_raw_store());
			for (size_t i = 0; i < res->get_length(); i++)
				fprintf(f, "\"%ld\" %f\n", i, res_store->get<float>(i));
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
	fm::vector::ptr comp_ids = compute_ts_wcc(graph, start_time,
			time_interval);
	if (comp_ids)
		print_cc(comp_ids);
}

void run_kcore(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	int num_opts = 0;
	size_t kmax = 0;
	size_t k = 2;
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

	fm::vector::ptr kcorev = compute_kcore(graph, k, kmax);
//	if (!write_out.empty() && kcorev)
//		kcorev->to_file(write_out);
}

void run_betweenness_centrality(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	int num_opts = 0;
	std::string write_out = "";
	vertex_id_t id = INVALID_VERTEX_ID;

	while ((opt = getopt(argc, argv, "w:s:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'w':
				write_out = optarg;
				break;
			case 's':
				id = atol(optarg);
				break;
			default:
				print_usage();
				assert(0);
		}
	}

	std::vector<vertex_id_t> ids;

	if (id == INVALID_VERTEX_ID) {
		for (vertex_id_t id = 0; id < graph->get_graph_header().get_num_vertices(); id++) {
			ids.push_back(id);
		}
	} else {
		ids.push_back(id);
	}

	fm::vector::ptr btwn_v = compute_betweenness_centrality(graph, ids);
//	if (!write_out.empty() && btwn_v)
//		btwn_v->to_file(write_out);
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

#if 0
void run_louvain(FG_graph::ptr graph, int argc, char* argv[])
{
	int opt;
	int num_opts = 0;
	uint32_t levels = 1;

	while ((opt = getopt(argc, argv, "l:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'l':
				levels = atoi(optarg);
				break;
			default:
				print_usage();
				assert(0);
		}
	}

	compute_louvain(graph, levels);
}
#endif

void run_sem_kmeans(FG_graph::ptr graph, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
    unsigned k = 0;
    unsigned max_iters = std::numeric_limits<unsigned>::max();
    std::string init = "kmeanspp";
    double tolerance = -1;

	while ((opt = getopt(argc, argv, "k:i:t:l:c:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'k':
				k = atol(optarg);
				num_opts++;
				break;
			case 'i':
                max_iters = atol(optarg);
				break;
			case 't':
                init = optarg;
				break;
            case 'l':
                tolerance = atof(optarg);
                num_opts++;
                break;
			default:
				print_usage();
				abort();
		}
	}

    compute_sem_kmeans(graph, k, init, max_iters, tolerance);
}

void run_spmm(FG_graph::ptr fg, int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	std::string entry_type;
	int num_cols = 1;

	while ((opt = getopt(argc, argv, "e:c:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'e':
				entry_type = optarg;
				num_opts++;
				break;
			case 'c':
				num_cols = atoi(optarg);
				num_opts++;
				break;
			default:
				print_usage();
				abort();
		}
	}

	fm::sparse_matrix::ptr mat;
	if (entry_type.empty())
		mat = create_sparse_matrix(fg, NULL);
	else if (entry_type == "I")
		mat = create_sparse_matrix(fg, &fm::get_scalar_type<int>());
	else if (entry_type == "L")
		mat = create_sparse_matrix(fg, &fm::get_scalar_type<long>());
	else if (entry_type == "F")
		mat = create_sparse_matrix(fg, &fm::get_scalar_type<float>());
	else if (entry_type == "D")
		mat = create_sparse_matrix(fg, &fm::get_scalar_type<double>());
	else {
		fprintf(stderr, "unknown entry type\n");
		return;
	}
	int num_nodes = safs::params.get_num_nodes();
	fm::dense_matrix::ptr in = fm::dense_matrix::create_randu<double>(0, 1,
			mat->get_num_cols(), num_cols, fm::matrix_layout_t::L_ROW, num_nodes);
	fm::detail::matrix_store::ptr out = fm::detail::mem_matrix_store::create(
			mat->get_num_rows(), num_cols, fm::matrix_layout_t::L_ROW,
			fm::get_scalar_type<double>(), num_nodes);
	mat->multiply(in->get_raw_store(), out);
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
	"betweenness",
	"overlap",
	"bfs",
	"louvain",
    "sem_kmeans"
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
	fprintf(stderr, "kcore\n");
	fprintf(stderr, "-k k: the minimum k value to compute\n");
	fprintf(stderr, "-m kmax: the maximum k value to compute\n");
	fprintf(stderr, "-w output: the file name for a vector written to file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "betweenness\n");
	fprintf(stderr, "-w output: the file name for a vector written to file\n");
	fprintf(stderr, "-s vertex id: the vertex where BC starts. (Default runs all)\n");
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
	fprintf(stderr, "louvain\n");
	fprintf(stderr, "-l: how many levels in the hierarchy to compute\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "sem_kmeans\n");
	fprintf(stderr, "-k: the number of clusters to use\n");
	fprintf(stderr, "-i: max number of iterations\n");
	fprintf(stderr, "-t: init type [random, forgy, kmeanspp]. Default: kmeanspp\n");
	fprintf(stderr, "-l: convergence tolerance (defualt: -1 = no changes)\n");
	fprintf(stderr, "\n");

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

	graph_engine::init_flash_graph(configs);
	FG_graph::ptr graph;
	try {
		graph = FG_graph::create(graph_file, index_file, configs);
	} catch(std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		exit(-1);
	}

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
	else if (alg == "betweenness") {
		run_betweenness_centrality(graph, argc, argv);
	}
	else if (alg == "overlap") {
		run_overlap(graph, argc, argv);
	}
	else if (alg == "bfs") {
		run_bfs(graph, argc, argv);
	}
#if 0
	else if (alg == "louvain") {
		run_louvain(graph, argc, argv);
	}
#endif
	else if (alg == "sem_kmeans") {
		run_sem_kmeans(graph, argc, argv);
	}
	else if (alg == "spmm") {
		run_spmm(graph, argc, argv);
	}
	else {
		fprintf(stderr, "\n[ERROR]: Unknown algorithm '%s'!\n", alg.c_str());
	}
	graph_engine::destroy_flash_graph();
}
