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

#include <unordered_map>
#include <boost/filesystem.hpp>
#include <Rcpp.h>

#include "log.h"
#include "safs_file.h"
#include "matrix/FG_sparse_matrix.h"
#include "matrix/matrix_eigensolver.h"

#include "FGlib.h"
#include "utils.h"
#include "in_mem_storage.h"

#if 0
FG_vector<float>::ptr compute_sstsg(FG_graph::ptr fg, time_t start_time,
		time_t interval, int num_intervals);
FG_vector<float>::ptr compute_betweenness_centrality(FG_graph::ptr fg, vertex_id_t id);
FG_vector<vsize_t>::ptr get_ts_degree(FG_graph::ptr fg, edge_type type,
		time_t start_time, time_t time_interval);
#endif

/*
 * A global configuration of FlashGraphR.
 */
static config_map::ptr configs;

/*
 * This class maintains a reference to an in-memory graph.
 */
class graph_ref
{
	in_mem_graph::ptr g;
	vertex_index::ptr index;
	std::string name;
	int count;
public:
	graph_ref(in_mem_graph::ptr g, vertex_index::ptr index,
			const std::string &name) {
		this->g = g;
		this->index = index;
		this->name = name;
		count = 1;
	}

	int get_counts() {
		return count;
	}

	FG_graph::ptr get_graph() {
		return FG_graph::create(g, index, name, configs);
	}

	const std::string &get_name() const {
		return name;
	}

	void ref() {
		count++;
		printf("inc_ref: graph %s has %d refs\n", name.c_str(), count);
	}

	void deref() {
		count--;
		printf("dec_ref: graph %s has %d refs\n", name.c_str(), count);
	}
};

typedef std::unordered_map<std::string, graph_ref *> graph_map_t;
/*
 * This contains all in-memory graphs loaded to FlashGraphR.
 */
static graph_map_t graphs;

bool standalone = true;

static bool exist_cindex(const std::string &graph_name)
{
	if (standalone)
		return false;

	std::string cindex_name = graph_name + "-cindex-v" + itoa(CURR_VERSION);
	safs_file cindex_file(get_sys_RAID_conf(), cindex_name);
	return cindex_file.exist();
}

/*
 * Get a FG_graph for the specified graph.
 */
static FG_graph::ptr R_FG_get_graph(SEXP pgraph)
{
	Rcpp::List graph(pgraph);
	// If the pointer field is defined, we can get the FG_graph object
	// directly.
	if (graph.containsElementNamed("pointer")) {
		graph_ref *ref = (graph_ref *) R_ExternalPtrAddr(graph["pointer"]);
		return ref->get_graph();
	}
	else if (standalone) {
		fprintf(stderr, "Wrong state! Can't get a graph\n");
		return FG_graph::ptr();
	}
	else {
		std::string graph_name = graph["name"];
		std::string version = itoa(CURR_VERSION);
		std::string graph_file = graph_name + "-v" + version;
		std::string index_file;
		if (exist_cindex(graph_name))
			index_file = graph_name + "-cindex-v" + version;
		else
			index_file = graph_name + "-index-v" + version;
		return FG_graph::create(graph_file, index_file, configs);
	}
}

/**
 * Initialize FlashGraphR.
 */
RcppExport SEXP R_FG_init(SEXP pconf)
{
	set_log_level(c_log_level::warning);
	std::string conf_file = CHAR(STRING_ELT(pconf, 0));

	boost::filesystem::path p(conf_file);
	if (boost::filesystem::exists(p)) {
		configs = config_map::create(conf_file);
	}
	else {
		fprintf(stderr, "conf file %s doesn't exist.\n", conf_file.c_str());
		configs = config_map::create();
	}

	Rcpp::LogicalVector res(1);
	try {
		graph_engine::init_flash_graph(configs);
		standalone = false;
		res[0] = true;
	} catch (init_error &e) {
		fprintf(stderr, "init FlashGraphR: %s\n", e.what());
		res[0] = true;
	} catch (std::exception &e) {
		fprintf(stderr, "exception in init: %s\n", e.what());
		res[0] = false;
	}

	if (standalone)
		printf("Run FlashGraphR in standalone mode\n");
	else if (is_safs_init())
		printf("Run FlashGraphR\n");
	else {
		fprintf(stderr, "Can't enable the SAFS mode of FlashGraphR\n");
		res[0] = false;
	}
	return res;
}

/**
 * Destroy FlashGraphR
 */
RcppExport SEXP R_FG_destroy()
{
	int num_refs = 0;
	for (auto it = graphs.begin(); it != graphs.end(); it++) {
		if (it->second->get_counts() == 1) {
			delete it->second;
			it->second = NULL;
		}
		else {
			num_refs++;
			fprintf(stderr, "%s is still referenced\n",
					it->second->get_name().c_str());
		}
	}
	if (num_refs == 0)
		graphs.clear();
	graph_engine::destroy_flash_graph();
	return R_NilValue;
}

/*
 * This search in memory first, and then searches in SAFS.
 */
static bool exist_graph(std::string &graph_name)
{
	// Search in memory.
	graph_map_t::const_iterator it = graphs.find(graph_name);
	if (it != graphs.end())
		return true;

	// If FlashGraphR runs in the standalone mode, we can't search in SAFS.
	if (standalone)
		return false;

	std::string version = itoa(CURR_VERSION);
	std::string graph_file_name = graph_name + "-v" + version;
	safs_file graph_file(get_sys_RAID_conf(), graph_file_name);
	if (!graph_file.exist()) {
		fprintf(stderr, "The graph file of %s doesn't exist\n",
				graph_name.c_str());
		return false;
	}

	std::string graph_index_name = graph_name + "-index-v" + version;
	std::string graph_cindex_name = graph_name + "-cindex-v" + version;
	safs_file graph_index_file(get_sys_RAID_conf(), graph_index_name);
	safs_file graph_cindex_file(get_sys_RAID_conf(), graph_cindex_name);
	if (!graph_index_file.exist() && !graph_cindex_file.exist()) {
		fprintf(stderr, "The index file of %s doesn't exist\n",
				graph_name.c_str());
		return false;
	}

	return true;
}

static SEXP get_safs_params()
{
	Rcpp::List ret;
	ret["RAID_block_size"] = params.get_RAID_block_size();
	ret["SA_min_cell_size"] = params.get_SA_min_cell_size();
	ret["IO_dpeth"] = params.get_aio_depth_per_file();
	ret["cache_type"] = params.get_cache_type();
	ret["cache_size"] = params.get_cache_size();
	ret["RAID_mapping"] = params.get_RAID_mapping_option();
	ret["virtual_AIO"] = params.is_use_virt_aio();
	ret["use_flusher"] = params.is_use_flusher();
	ret["NUMA_num_process_threads"] = params.get_numa_num_process_threads();
	ret["num_nodes"] = params.get_num_nodes();
	ret["merge_requests"] = params.is_merge_reqs();
	ret["max_obj_alloc_size"] = params.get_max_obj_alloc_size();
	ret["writable"] = params.is_writable();
	ret["max_num_pending_IOs"] = params.get_max_num_pending_ios();
	ret["huge_page"] = params.is_huge_page_enabled();
	return ret;
}

static SEXP get_fg_params()
{
	Rcpp::List ret;
	ret["prof_file"] = graph_conf.get_prof_file();
	ret["num_threads"] = graph_conf.get_num_threads();
	ret["elevator"] = graph_conf.get_elevator_enabled();
	ret["max_processing_vertices"] = graph_conf.get_max_processing_vertices();
	ret["part_range_size_log"] = graph_conf.get_part_range_size_log();
	ret["preload"] = graph_conf.preload();
	ret["index_file_weight"] = graph_conf.get_index_file_weight();
	ret["in_mem_index"] = graph_conf.use_in_mem_index();
	ret["in_mem_graph"] = graph_conf.use_in_mem_graph();
	ret["serial_run"] = graph_conf.use_serial_run();
	ret["num_vertical_parts"] = graph_conf.get_num_vparts();
	ret["min_vpart_degree"] = graph_conf.get_min_vpart_degree();
	return ret;
}

/**
 * This returns all parameters of SAFS or FlashGraph.
 */
RcppExport SEXP R_FG_get_params(SEXP psys)
{
	std::string sys_name = CHAR(STRING_ELT(psys, 0));
	if (sys_name == "SAFS") {
		if (standalone) {
			fprintf(stderr,
					"Can't get SAFS parameters. FlashGraphR runs in standalone mode\n");
			return R_NilValue;
		}
		else
			return get_safs_params();
	}
	else if (sys_name == "FlashGraph")
		return get_fg_params();
	else {
		fprintf(stderr, "wrong system name\n");
		return R_NilValue;
	}
}

/**
 * This test whether a graph has been loaded to FlashGraphR.
 */
RcppExport SEXP R_FG_exist_graph(SEXP pgraph)
{
	std::string graph_name = CHAR(STRING_ELT(pgraph, 0));
	Rcpp::LogicalVector res(1);
	res[0] = exist_graph(graph_name);
	return res;
}

static std::string extract_graph_name(const std::string &file_name)
{
	std::string version = itoa(CURR_VERSION);
	size_t pos = file_name.rfind("-cindex-v" + version);
	if (pos == std::string::npos)
		pos = file_name.rfind("-index-v" + version);
	if (pos == std::string::npos)
		pos = file_name.rfind("-v" + version);
	if (pos == std::string::npos)
		return "";
	else
		return file_name.substr(0, pos);
}

/**
 * This lists all graphs that have been loaded to FlashGraphR.
 */
RcppExport SEXP R_FG_list_graphs()
{
	// Get all graphs in memory.
	std::map<std::string, bool> graph_names;
	for (graph_map_t::const_iterator it = graphs.begin();
			it != graphs.end(); it++) {
		graph_names.insert(std::pair<std::string, bool>(it->first, true));
	}

	// Get all graphs in SAFS.
	if (!standalone) {
		std::set<std::string> files;
		get_all_safs_files(files);

		BOOST_FOREACH(std::string file, files) {
			std::string graph_name = extract_graph_name(file);
			if (!graph_name.empty())
				graph_names.insert(std::pair<std::string, bool>(
							graph_name, false));
		}
	}

	Rcpp::CharacterVector names;
	Rcpp::LogicalVector in_mem;
	typedef std::pair<std::string, bool> graph_value_t;
	BOOST_FOREACH(graph_value_t graph, graph_names) {
		if (graph.second) {
			names.push_back(graph.first);
			in_mem.push_back(graph.second);
		}
		else if (exist_graph(graph.first)) {
			names.push_back(graph.first);
			in_mem.push_back(graph.second);
		}
	}
	return Rcpp::DataFrame::create(Rcpp::Named("name") = names,
			Rcpp::Named("in-mem") = in_mem);
}

RcppExport SEXP R_FG_set_log_level(SEXP plevel)
{
	std::string level = CHAR(STRING_ELT(plevel, 0));
	if (level == "debug") {
		set_log_level(c_log_level::debug);
	}
	else if (level == "info") {
		set_log_level(c_log_level::info);
	}
	else if (level == "warning") {
		set_log_level(c_log_level::warning);
	}
	else if (level == "error") {
		set_log_level(c_log_level::error);
	}
	else if (level == "fatal") {
		set_log_level(c_log_level::fatal);
	}
	else {
		fprintf(stderr, "unknown level %s\n", level.c_str());
	}
	return R_NilValue;
}

static void fg_clean_graph(SEXP p)
{
	graph_ref *ref = (graph_ref *) R_ExternalPtrAddr(p);
	ref->deref();
}

static SEXP create_FGR_obj(graph_ref *ref)
{
	std::string graph_name = ref->get_name();
	Rcpp::List ret;
	ret["name"] = Rcpp::String(graph_name);

	FG_graph::ptr graph = ref->get_graph();
	if (graph == NULL) {
		fprintf(stderr, "the graph table has inconsistent data.\n");
		return R_NilValue;
	}

	ref->ref();
	SEXP pointer = R_MakeExternalPtr(ref, R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(pointer, fg_clean_graph, FALSE);
	ret["pointer"] = pointer;

	graph_header header = graph->get_graph_header();
	Rcpp::LogicalVector directed(1);
	directed[0] = header.is_directed_graph();
	ret["directed"] = directed;

	Rcpp::NumericVector vcount(1);
	vcount[0] = header.get_num_vertices();
	ret["vcount"] = vcount;

	Rcpp::NumericVector ecount(1);
	ecount[0] = header.get_num_edges();
	ret["ecount"] = ecount;

	Rcpp::LogicalVector in_mem(1);
	in_mem[0] = graph->get_graph_data() != NULL;
	ret["in.mem"] = in_mem;
	return ret;
}

static SEXP create_FGR_obj(FG_graph::ptr graph, const std::string &graph_name)
{
	Rcpp::List ret;
	ret["name"] = Rcpp::String(graph_name);

	graph_header header = graph->get_graph_header();
	Rcpp::LogicalVector directed(1);
	directed[0] = header.is_directed_graph();
	ret["directed"] = directed;

	Rcpp::NumericVector vcount(1);
	vcount[0] = header.get_num_vertices();
	ret["vcount"] = vcount;

	Rcpp::NumericVector ecount(1);
	ecount[0] = header.get_num_edges();
	ret["ecount"] = ecount;

	Rcpp::LogicalVector in_mem(1);
	in_mem[0] = graph->get_graph_data() != NULL;
	ret["in.mem"] = in_mem;
	return ret;
}

graph_ref *register_in_mem_graph(FG_graph::ptr fg,
		const std::string &graph_name)
{
	if (!fg->is_in_mem())
		return NULL;

	graph_ref *ref = new graph_ref(fg->get_graph_data(), fg->get_index_data(),
			graph_name);
	auto ret = graphs.insert(std::pair<std::string, graph_ref *>(graph_name,
				ref));
	if (!ret.second)
		ret.first->second = ref;
	return ref;
}

RcppExport SEXP R_FG_load_graph_adj(SEXP pgraph_name, SEXP pgraph_file,
		SEXP pindex_file)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	std::string graph_file = CHAR(STRING_ELT(pgraph_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));
	FG_graph::ptr fg = FG_graph::create(graph_file, index_file, configs);
	graph_ref *ref = register_in_mem_graph(fg, graph_name);
	if (ref)
		return create_FGR_obj(ref);
	else
		return create_FGR_obj(fg, graph_name);
}

/*
 * Load a graph from edge lists in a data frame.
 */
RcppExport SEXP R_FG_load_graph_el_df(SEXP pgraph_name, SEXP pedge_lists,
		SEXP pdirected, SEXP pnthreads)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	Rcpp::DataFrame edge_lists = Rcpp::DataFrame(pedge_lists);
	bool directed = INTEGER(pdirected)[0];
	int num_threads = INTEGER(pnthreads)[0];

	Rcpp::IntegerVector from = edge_lists["from"];
	Rcpp::IntegerVector to = edge_lists["to"];
	std::vector<vertex_id_t> from_vec(from.begin(), from.end());
	std::vector<vertex_id_t> to_vec(to.begin(), to.end());

	std::pair<in_mem_graph::ptr, vertex_index::ptr> gpair = construct_mem_graph(
		from_vec, to_vec, graph_name, DEFAULT_TYPE, directed, num_threads);
	FG_graph::ptr fg = FG_graph::create(gpair.first, gpair.second, graph_name,
			configs);
	graph_ref *ref = register_in_mem_graph(fg, graph_name);
	if (ref)
		return create_FGR_obj(ref);
	else
		return create_FGR_obj(fg, graph_name);
}

/*
 * Load a graph from edge lists in a file.
 */
RcppExport SEXP R_FG_load_graph_el(SEXP pgraph_name, SEXP pgraph_file,
		SEXP pdirected, SEXP pnthreads)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	std::string graph_file = CHAR(STRING_ELT(pgraph_file, 0));
	bool directed = INTEGER(pdirected)[0];
	int num_threads = INTEGER(pnthreads)[0];

	native_file f(graph_file);
	if (!f.exist()) {
		fprintf(stderr, "%s doesn't exist\n", graph_file.c_str());
		return R_NilValue;
	}

	std::vector<std::string> edge_list_files(1);
	edge_list_files[0] = graph_file;
	std::pair<in_mem_graph::ptr, vertex_index::ptr> gpair = construct_mem_graph(
		edge_list_files, graph_name, DEFAULT_TYPE, directed, num_threads);
	FG_graph::ptr fg = FG_graph::create(gpair.first, gpair.second, graph_name,
			configs);
	graph_ref *ref = register_in_mem_graph(fg, graph_name);
	if (ref)
		return create_FGR_obj(ref);
	else
		return create_FGR_obj(fg, graph_name);
}

RcppExport SEXP R_FG_get_graph_obj(SEXP pgraph)
{
	std::string graph_name = CHAR(STRING_ELT(pgraph, 0));
	if (!exist_graph(graph_name)) {
		fprintf(stderr, "%s doesn't exist\n", graph_name.c_str());
		return R_NilValue;
	}

	auto it = graphs.find(graph_name);
	// If the graph exist, but it's not in the graph table. It's in SAFS.
	if (it == graphs.end()) {
		std::string version = itoa(CURR_VERSION);
		std::string graph_file = graph_name + "-v" + version;
		std::string index_file;

		std::string cindex_name = graph_name + "-cindex-v" + itoa(CURR_VERSION);
		safs_file cindex_file(get_sys_RAID_conf(), cindex_name);
		if (cindex_file.exist())
			index_file = graph_name + "-cindex-v" + version;
		else
			index_file = graph_name + "-index-v" + version;
		FG_graph::ptr fg = FG_graph::create(graph_file, index_file, configs);
		graph_ref *ref = register_in_mem_graph(fg, graph_name);
		if (ref)
			return create_FGR_obj(ref);
		else
			return create_FGR_obj(fg, graph_name);
	}
	else
		return create_FGR_obj(it->second);

}

///////////////////////////// graph algorithms ///////////////////////////

RcppExport SEXP R_FG_compute_cc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<vertex_id_t>::ptr fg_vec = compute_cc(fg);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_wcc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<vertex_id_t>::ptr fg_vec = compute_wcc(fg);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_scc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<vertex_id_t>::ptr fg_vec = compute_scc(fg);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_transitivity(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<float>::ptr fg_vec = compute_transitivity(fg);
	Rcpp::NumericVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_get_degree(SEXP graph, SEXP ptype)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);

	std::string type_str = CHAR(STRING_ELT(ptype, 0));
	edge_type type = edge_type::NONE;
	if (type_str == "in")
		type = edge_type::IN_EDGE;
	else if (type_str == "out")
		type = edge_type::OUT_EDGE;
	else if (type_str == "both")
		type = edge_type::BOTH_EDGES;
	else {
		fprintf(stderr, "wrong edge type\n");
		return R_NilValue;
	}

	FG_vector<vsize_t>::ptr fg_vec = get_degree(fg, type);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_pagerank(SEXP graph, SEXP piters, SEXP pdamping)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);

	int num_iters = REAL(piters)[0];
	float damping_factor = REAL(pdamping)[0];

	FG_vector<float>::ptr fg_vec = compute_pagerank2(fg, num_iters, damping_factor);
	Rcpp::NumericVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_undirected_triangles(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<size_t>::ptr fg_vec = compute_undirected_triangles(fg);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_directed_triangles(SEXP graph, SEXP ptype)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);

	std::string type_str = CHAR(STRING_ELT(ptype, 0));
	directed_triangle_type type = directed_triangle_type::CYCLE;
	if (type_str == "cycle")
		type = directed_triangle_type::CYCLE;
	else
		type = directed_triangle_type::ALL;

	FG_vector<size_t>::ptr fg_vec = compute_directed_triangles_fast(fg, type);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_local_scan(SEXP graph, SEXP order)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<size_t>::ptr fg_vec = compute_local_scan(fg);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_topK_scan(SEXP graph, SEXP order, SEXP K)
{
	size_t topK = REAL(K)[0];
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<std::pair<vertex_id_t, size_t> >::ptr fg_vec
		= compute_topK_scan(fg, topK);
	assert(fg_vec->get_size() == topK);
	Rcpp::IntegerVector vertices(fg_vec->get_size());
	Rcpp::IntegerVector scans(fg_vec->get_size());
	for (size_t i = 0; i < topK; i++) {
		std::pair<vertex_id_t, size_t> pair = fg_vec->get(i);
		vertices[i] = pair.first;
		scans[i] = pair.second;
	}
	return Rcpp::DataFrame::create(Named("vid", vertices), Named("scan", scans));
}

RcppExport SEXP R_FG_compute_kcore(SEXP graph, SEXP _k, SEXP _kmax)
{
	int k = REAL(_k)[0];
	int kmax = REAL(_kmax)[0];
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_vector<size_t>::ptr fg_vec = compute_kcore(fg, k, kmax);
	Rcpp::IntegerVector res(fg_vec->get_size());
	fg_vec->copy_to(res.begin(), fg_vec->get_size());
	return res;
}

RcppExport SEXP R_FG_compute_overlap(SEXP graph, SEXP _vids)
{
	Rcpp::IntegerVector Rvids(_vids);
	std::vector<vertex_id_t> vids(Rvids.begin(), Rvids.end());
	std::vector<std::vector<double> > overlap_matrix;
	size_t num_vertices = vids.size();

	FG_graph::ptr fg = R_FG_get_graph(graph);
	compute_overlap(fg, vids, overlap_matrix);
	assert(overlap_matrix.size() == num_vertices);

	Rcpp::NumericMatrix res(num_vertices, num_vertices);
	for (size_t i = 0; i < num_vertices; i++) {
		assert(overlap_matrix[i].size() == num_vertices);
		for (size_t j = 0; j < num_vertices; j++) {
			res(i, j) = overlap_matrix[i][j];
		}
	}
	return res;
}

/*
 * Generate an induced subgraph from a graph, given a list of vertices,
 * and convert it into an edge list.
 */
RcppExport SEXP R_FG_fetch_subgraph_el(SEXP graph, SEXP pvertices)
{
	Rcpp::IntegerVector vertices(pvertices);
	std::vector<vertex_id_t> vids(vertices.begin(), vertices.end());
	FG_graph::ptr fg = R_FG_get_graph(graph);
	in_mem_subgraph::ptr subg = fetch_subgraph(fg, vids);
	subg->compress();
	assert(subg->get_num_vertices() == vids.size());
	Rcpp::IntegerVector s_vs(subg->get_num_edges());
	Rcpp::IntegerVector d_vs(subg->get_num_edges());
	size_t num_tot_edges = 0;
	BOOST_FOREACH(vertex_id_t id, vids) {
		const in_mem_vertex &v = subg->get_vertex(id);
		if (v.has_edge_data())
			ABORT_MSG("we can't fetch a subgraph from a graph with attributes");
		if (subg->is_directed()) {
			const in_mem_directed_vertex<> &dv
				= (const in_mem_directed_vertex<> &) v;
			size_t num_edges = dv.get_num_out_edges();
			for (size_t i = 0; i < num_edges; i++) {
				edge<> e = dv.get_out_edge(i);
				s_vs[num_tot_edges] = e.get_from();
				d_vs[num_tot_edges] = e.get_to();
				num_tot_edges++;
			}
		}
		else {
			const in_mem_undirected_vertex<> &un_v
				= (const in_mem_undirected_vertex<> &) v;
			size_t num_edges = un_v.get_num_edges();
			for (size_t i = 0; i < num_edges; i++) {
				edge<> e = un_v.get_edge(i);
				// each edge appears twice in an undirected graph.
				// we only need to store one.
				if (e.get_from() <= e.get_to()) {
					s_vs[num_tot_edges] = e.get_from();
					d_vs[num_tot_edges] = e.get_to();
					num_tot_edges++;
				}
			}
		}
	}
	assert(s_vs.size() == num_tot_edges);
	Rcpp::List ret;
	ret["src"] = s_vs;
	ret["dst"] = d_vs;
	return ret;
}

RcppExport SEXP R_FG_fetch_subgraph(SEXP graph, SEXP pvertices, SEXP pname)
{
	std::string graph_name = CHAR(STRING_ELT(pname, 0));
	Rcpp::IntegerVector vertices(pvertices);
	std::vector<vertex_id_t> vids(vertices.begin(), vertices.end());
	FG_graph::ptr fg = R_FG_get_graph(graph);
	in_mem_subgraph::ptr subg = fetch_subgraph(fg, vids);
	assert(subg->get_num_vertices() == vids.size());
	std::pair<in_mem_graph::ptr, vertex_index::ptr> gpair
		= subg->compress(graph_name);
	FG_graph::ptr sub_fg = FG_graph::create(gpair.first, gpair.second,
			graph_name, configs);
	graph_ref *ref = register_in_mem_graph(sub_fg, graph_name);
	if (ref)
		return create_FGR_obj(ref);
	else
		return create_FGR_obj(sub_fg, graph_name);
}

RcppExport SEXP R_FG_estimate_diameter(SEXP graph, SEXP pdirected)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	bool directed = INTEGER(pdirected)[0];
	int diameter = estimate_diameter(fg, 1, directed);
	Rcpp::IntegerVector ret(1);
	ret[0] = diameter;
	return ret;
}

RcppExport SEXP R_FG_multiply_v(SEXP graph, SEXP pvec, SEXP ptranspose)
{
	bool transpose = INTEGER(ptranspose)[0];
	Rcpp::NumericVector vec(pvec);
	size_t length = vec.size();
	FG_vector<double>::ptr in_vec = FG_vector<double>::create(length);
	for (size_t i = 0; i < length; i++) {
		in_vec->get_data()[i] = vec[i];
	}
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_adj_matrix::ptr matrix = FG_adj_matrix::create(fg);
	if (transpose)
		matrix = matrix->transpose();
	assert(matrix->get_num_rows() == length);
	assert(matrix->get_num_cols() == length);

	FG_vector<double>::ptr out_vec = FG_vector<double>::create(length);
	matrix->multiply<double>(*in_vec, *out_vec);
	Rcpp::NumericVector ret(out_vec->get_data(), out_vec->get_data() + length);
	return ret;
}

#ifdef USE_EIGEN
/*
 * Compute eigen value/vector on an unweighted adjacency matrix.
 */
RcppExport SEXP R_FG_eigen_uw(SEXP graph, SEXP pwhich, SEXP pnev, SEXP pncv)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_adj_matrix::ptr matrix = FG_adj_matrix::create(fg);
	std::string which = CHAR(STRING_ELT(pwhich, 0));
	int nev = INTEGER(pnev)[0];
	int ncv = INTEGER(pncv)[0];
	std::vector<eigen_pair_t> eigen_pairs;
	compute_eigen<FG_adj_matrix>(matrix, ncv, nev, which, eigen_pairs);
	if (eigen_pairs.empty())
		return R_NilValue;

	size_t length = eigen_pairs[0].second->get_size();
	Rcpp::NumericVector eigen_values(nev);
	Rcpp::NumericMatrix eigen_matrix(length, nev);
	for (int i = 0; i < nev; i++) {
		eigen_values[i] = eigen_pairs[i].first;
		for (size_t j = 0; j < length; j++)
			eigen_matrix(j, i) = eigen_pairs[i].second->get(j);
	}
	Rcpp::List ret;
	ret["values"] = eigen_values;
	ret["vectors"] = eigen_matrix;
	return ret;
}

RcppExport SEXP R_FG_SVD_uw(SEXP graph, SEXP pwhich, SEXP pnev, SEXP pncv,
		SEXP ptype)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	FG_adj_matrix::ptr matrix = FG_adj_matrix::create(fg);
	std::string which = CHAR(STRING_ELT(pwhich, 0));
	std::string type = CHAR(STRING_ELT(ptype, 0));
	int nev = INTEGER(pnev)[0];
	int ncv = INTEGER(pncv)[0];
	std::vector<eigen_pair_t> eigen_pairs;
	compute_SVD<FG_adj_matrix>(matrix, ncv, nev, which, type, eigen_pairs);
	if (eigen_pairs.empty())
		return R_NilValue;

	size_t length = eigen_pairs[0].second->get_size();
	Rcpp::NumericVector eigen_values(nev);
	Rcpp::NumericMatrix eigen_matrix(length, nev);
	for (int i = 0; i < nev; i++) {
		eigen_values[i] = eigen_pairs[i].first;
		for (size_t j = 0; j < length; j++)
			eigen_matrix(j, i) = eigen_pairs[i].second->get(j);
	}
	Rcpp::List ret;
	ret["values"] = eigen_values;
	ret["vectors"] = eigen_matrix;
	return ret;
}
#endif
