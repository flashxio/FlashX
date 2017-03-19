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
#include <Rcpp.h>

#include "log.h"
#include "safs_file.h"
#include "io_interface.h"

#include "matrix/kmeans.h"
#include "libgraph-algs/sem_kmeans.h"

#include "FGlib.h"
#include "utils.h"
#include "in_mem_storage.h"
#include "fg_utils.h"

#include "mem_vec_store.h"
#include "data_frame.h"
#include "data_io.h"
#include "sparse_matrix.h"

#include "rutils.h"

using namespace safs;
using namespace fg;

/*
 * A global configuration of FlashGraph.
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
	}

	void deref() {
		count--;
	}
};

typedef std::unordered_map<std::string, graph_ref *> graph_map_t;
/*
 * This contains all in-memory graphs loaded to FlashGraphR.
 */
static graph_map_t graphs;

static bool standalone = true;

static std::pair<std::string, std::string> get_graph_files(
		const std::string &graph_name)
{
	std::string graph_file = graph_name + ".adj";
	std::string index_file = graph_name + ".index";
	return std::pair<std::string, std::string>(graph_file, index_file);
}

/*
 * Get a FG_graph for the specified graph.
 */
FG_graph::ptr R_FG_get_graph(SEXP pgraph)
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
		auto graph_files = get_graph_files(graph["name"]);
		return FG_graph::create(graph_files.first, graph_files.second, configs);
	}
}

/**
 * Initialize FlashGraph.
 */
RcppExport SEXP R_FG_init(SEXP pconf)
{
	set_log_level(c_log_level::warning);
	std::string conf_file;
	if (!R_is_null(pconf) && R_is_string(pconf))
		conf_file = CHAR(STRING_ELT(pconf, 0));

	if (!conf_file.empty() && safs::file_exist(conf_file)) {
		configs = config_map::create(conf_file);
		configs->add_options("writable=1");
	}
	else if (!conf_file.empty()) {
		fprintf(stderr, "conf file %s doesn't exist.\n", conf_file.c_str());
		configs = config_map::create();
	}
	// If there isn't a conf file, we just use the default settings.
	else
		configs = config_map::create();

	standalone = !is_safs_init();
	bool fg_success;
	try {
		graph_engine::init_flash_graph(configs);
		fg_success = true;
	} catch (std::exception &e) {
		fprintf(stderr, "exception in init: %s\n", e.what());
		fg_success = false;
	}

	Rcpp::LogicalVector res(1);
	res[0] = fg_success;
	if (standalone)
		printf("Run FlashGraphR in standalone mode\n");
	else
		printf("Run FlashGraphR\n");
	return res;
}

/**
 * Destroy FlashGraphR
 */
#if 0
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
	fm::destroy_flash_matrix();
	safs::destroy_io_system();
	return R_NilValue;
}
#endif

RcppExport SEXP R_FG_destroy()
{
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

	auto graph_files = get_graph_files(graph_name);
	safs_file graph_file(get_sys_RAID_conf(), graph_files.first);
	if (!graph_file.exist())
		return false;

	safs_file graph_index_file(get_sys_RAID_conf(), graph_files.second);
	if (!graph_index_file.exist())
		return false;
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
	size_t pos = file_name.rfind(".adj");
	if (pos == std::string::npos)
		pos = file_name.rfind(".index");
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
	if (ref->get_counts() > 1)
		return;

	auto it = graphs.find(ref->get_name());
	if (it == graphs.end()) {
		fprintf(stderr, "graph %s doesn't exist\n", ref->get_name().c_str());
	}
	else {
		// If the graph is still registered in the graph table.
		if (it->second == ref)
			graphs.erase(it);
	}

	// There are no R objects referencing this graph now.
	printf("delete graph %s\n", ref->get_name().c_str());
	delete ref;
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
	in_mem[0] = graph->is_in_mem();
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
	in_mem[0] = graph->is_in_mem();
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
	if (!ret.second) {
		// If the in-memory graph isn't referenced by any R objects, we should
		// delete it.
		if (ret.first->second->get_counts() == 1) {
			printf("delete the old graph registered with %s\n", graph_name.c_str());
			delete ret.first->second;
		}
		ret.first->second = ref;
	}
	return ref;
}

RcppExport SEXP R_FG_load_graph_adj(SEXP pgraph_name, SEXP pgraph_file,
		SEXP pindex_file)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	std::string graph_file = CHAR(STRING_ELT(pgraph_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));

	FG_graph::ptr fg;
	try {
		fg = FG_graph::create(graph_file, index_file, configs);
	} catch(std::exception &e) {
		fprintf(stderr, "%s\n", e.what());
		return R_NilValue;
	}
	graph_ref *ref = register_in_mem_graph(fg, graph_name);
	if (ref)
		return create_FGR_obj(ref);
	else
		return create_FGR_obj(fg, graph_name);
}

RcppExport SEXP R_FG_export_graph(SEXP pgraph, SEXP pgraph_file, SEXP pindex_file)
{
	FG_graph::ptr fg = R_FG_get_graph(pgraph);
	std::string graph_file = CHAR(STRING_ELT(pgraph_file, 0));
	std::string index_file = CHAR(STRING_ELT(pindex_file, 0));

	Rcpp::LogicalVector ret(1);
	if (!fg->is_in_mem()) {
		fprintf(stderr, "currently we only support exporting in-mem graphs\n");
		ret[0] = false;
	}
	else {
		in_mem_graph::ptr graph_data = fg->get_graph_data();
		vertex_index::ptr index_data = fg->get_index_data();
		graph_data->dump(graph_file);
		index_data->dump(index_file);
		ret[0] = true;
	}
	return ret;
}

/*
 * Load a graph from edge lists in a data frame.
 */
RcppExport SEXP R_FG_load_graph_el_df(SEXP pgraph_name, SEXP pedge_lists,
		SEXP pdirected)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	Rcpp::DataFrame edge_lists = Rcpp::DataFrame(pedge_lists);
	bool directed = INTEGER(pdirected)[0];

	Rcpp::IntegerVector from = edge_lists["from"];
	Rcpp::IntegerVector to = edge_lists["to"];
	std::vector<vertex_id_t> from_vec(from.begin(), from.end());
	std::vector<vertex_id_t> to_vec(to.begin(), to.end());
	fm::detail::mem_vec_store::ptr from_store;
	fm::detail::mem_vec_store::ptr to_store;
	if (directed) {
		from_store = fm::detail::mem_vec_store::create(from.size(),
				-1, fm::get_scalar_type<vertex_id_t>());
		to_store = fm::detail::mem_vec_store::create(to.size(), -1,
				fm::get_scalar_type<vertex_id_t>());
		from_store->copy_from((const char *) from_vec.data(),
				from_vec.size() * sizeof(vertex_id_t));
		to_store->copy_from((const char *) to_vec.data(),
				to_vec.size() * sizeof(vertex_id_t));
	}
	else {
		from_store = fm::detail::mem_vec_store::create(from.size() * 2,
				-1, fm::get_scalar_type<vertex_id_t>());
		to_store = fm::detail::mem_vec_store::create(to.size() * 2, -1,
				fm::get_scalar_type<vertex_id_t>());
		size_t num_bytes = from_vec.size() * sizeof(vertex_id_t);
		memcpy(from_store->get_raw_arr(), from_vec.data(), num_bytes);
		memcpy(from_store->get_raw_arr() + num_bytes, to_vec.data(),
				num_bytes);
		memcpy(to_store->get_raw_arr(), to_vec.data(), num_bytes);
		memcpy(to_store->get_raw_arr() + num_bytes, from_vec.data(),
				num_bytes);
	}

	fm::data_frame::ptr df = fm::data_frame::create();
	df->add_vec("source", from_store);
	df->add_vec("dest", to_store);
	edge_list::ptr el = edge_list::create(df, directed);
	FG_graph::ptr fg = create_fg_graph(graph_name, el);

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
		SEXP pdirected, SEXP pin_mem, SEXP pdelim, SEXP pattr_type)
{
	Rcpp::LogicalVector res(1);
	std::string graph_name = CHAR(STRING_ELT(pgraph_name, 0));
	std::string graph_file = CHAR(STRING_ELT(pgraph_file, 0));
	bool directed = LOGICAL(pdirected)[0];
	bool in_mem = LOGICAL(pin_mem)[0];
	std::string delim = CHAR(STRING_ELT(pdelim, 0));
	std::string attr_type = CHAR(STRING_ELT(pattr_type, 0));

	if (!in_mem && !is_safs_init()) {
		fprintf(stderr, "SAFS isn't initialized\n");
		return R_NilValue;
	}

	native_file f(graph_file);
	if (!f.exist()) {
		fprintf(stderr, "edge list file %s doesn't exist\n", graph_file.c_str());
		return R_NilValue;
	}

	std::vector<std::string> edge_list_files(1);
	edge_list_files[0] = graph_file;
	// TODO give more options when loading an edge list.
	fm::data_frame::ptr df = utils::read_edge_list(edge_list_files, in_mem,
			delim, attr_type, directed);
	if (df == NULL)
		return R_NilValue;
	edge_list::ptr el = edge_list::create(df, directed);
	FG_graph::ptr fg = create_fg_graph(graph_name, el);
	if (fg == NULL)
		return R_NilValue;

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
		fprintf(stderr, "graph %s doesn't exist\n", graph_name.c_str());
		return R_NilValue;
	}

	auto it = graphs.find(graph_name);
	// If the graph exist, but it's not in the graph table. It's in SAFS.
	if (it == graphs.end()) {
		try {
			auto graph_files = get_graph_files(graph_name);
			FG_graph::ptr fg = FG_graph::create(graph_files.first,
					graph_files.second, configs);
			graph_ref *ref = register_in_mem_graph(fg, graph_name);
			if (ref)
				return create_FGR_obj(ref);
			else
				return create_FGR_obj(fg, graph_name);
		} catch(wrong_format &e) {
			fprintf(stderr, "%s\n", e.what());
			return R_NilValue;
		}
	}
	else
		return create_FGR_obj(it->second);

}

///////////////////////////// graph algorithms ///////////////////////////

enum R_type
{
	R_LOGICAL,
	R_INT,
	R_REAL,
	R_NTYPES,
};

SEXP create_FMR_vector(fm::dense_matrix::ptr m, R_type type, const std::string &name);

SEXP create_FMR_vector(fm::dense_matrix::ptr m, const std::string &name)
{
	R_type type;
	if (m->get_type() == fm::get_scalar_type<double>())
		type = R_type::R_REAL;
	else if (m->get_type() == fm::get_scalar_type<int>())
		type = R_type::R_INT;
	else {
		fprintf(stderr, "unknown type\n");
		return R_NilValue;
	}
	return create_FMR_vector(m, type, name);
}

fm::dense_matrix::ptr get_vertex_ids(fm::vector::ptr vec)
{
	fm::dense_matrix::ptr mat = vec->conv2mat(vec->get_length(), 1, false);
	// unsigned int isn't supported by FlashR. let's cast them to double.
	mat = mat->cast_ele_type(fm::get_scalar_type<double>());
	return mat;
}

template<class T>
fm::dense_matrix::ptr cast_type(fm::vector::ptr vec)
{
	fm::dense_matrix::ptr mat = vec->conv2mat(vec->get_length(), 1, false);
	return mat->cast_ele_type(fm::get_scalar_type<T>());
}

RcppExport SEXP R_FG_compute_cc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	fm::vector::ptr fg_vec = compute_cc(fg);
	return create_FMR_vector(get_vertex_ids(fg_vec), "");
}

RcppExport SEXP R_FG_compute_wcc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	fm::vector::ptr fg_vec = compute_wcc(fg);
	return create_FMR_vector(get_vertex_ids(fg_vec), "");
}

RcppExport SEXP R_FG_compute_scc(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	fm::vector::ptr fg_vec = compute_scc(fg);
	return create_FMR_vector(get_vertex_ids(fg_vec), "");
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

	fm::vector::ptr fg_vec = get_degree(fg, type);
	return create_FMR_vector(get_vertex_ids(fg_vec), "");
}

RcppExport SEXP R_FG_compute_pagerank(SEXP graph, SEXP piters, SEXP pdamping)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);

	int num_iters = REAL(piters)[0];
	float damping_factor = REAL(pdamping)[0];

	fm::vector::ptr fg_vec = compute_pagerank2(fg, num_iters, damping_factor);
	return create_FMR_vector(cast_type<double>(fg_vec), "");
}

RcppExport SEXP R_FG_compute_undirected_triangles(SEXP graph)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	fm::vector::ptr fg_vec = compute_undirected_triangles(fg);
	return create_FMR_vector(cast_type<double>(fg_vec), "");
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

	fm::vector::ptr fg_vec = compute_directed_triangles_fast(fg, type);
	return create_FMR_vector(cast_type<double>(fg_vec), "");
}

RcppExport SEXP R_FG_compute_local_scan(SEXP graph, SEXP porder)
{
	FG_graph::ptr fg = R_FG_get_graph(graph);
	int order = INTEGER(porder)[0];
	if (order == 0) {
		fm::vector::ptr fg_vec = get_degree(fg, edge_type::BOTH_EDGES);
		return create_FMR_vector(get_vertex_ids(fg_vec), "");
	}
	else if (order == 1) {
		fm::vector::ptr fg_vec = compute_local_scan(fg);
		return create_FMR_vector(cast_type<double>(fg_vec), "");
	}
	else if (order == 2) {
		fm::vector::ptr fg_vec = compute_local_scan2(fg);
		return create_FMR_vector(cast_type<double>(fg_vec), "");
	}
	else
		return R_NilValue;
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
	fm::vector::ptr fg_vec = compute_kcore(fg, k, kmax);
	return create_FMR_vector(cast_type<double>(fg_vec), "");
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

RcppExport SEXP R_FG_fetch_subgraph(SEXP graph, SEXP pvertices, SEXP pname,
		SEXP pcompress)
{
	bool compress = LOGICAL(pcompress)[0];
	std::string graph_name = CHAR(STRING_ELT(pname, 0));
	Rcpp::IntegerVector vertices(pvertices);
	if (vertices.length() == 0) {
		fprintf(stderr, "There aren't vertices to fetch\n");
		return R_NilValue;
	}
	std::vector<vertex_id_t> vids(vertices.begin(), vertices.end());

	FG_graph::ptr fg = R_FG_get_graph(graph);
	vertex_id_t max_vid = fg->get_graph_header().get_num_vertices() - 1;
	BOOST_FOREACH(vertex_id_t vid, vids) {
		if (vid > max_vid) {
			fprintf(stderr, "invalid vertex id: %d\n", vid);
			return R_NilValue;
		}
	}

	FG_graph::ptr sub_fg = fetch_subgraph(fg, vids, graph_name, compress);
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

RcppExport SEXP R_FG_sem_kmeans(SEXP graph, SEXP pk, SEXP pinit,
        SEXP pmax_iters, SEXP ptolerance)
{
    // Argparse
	FG_graph::ptr fg = R_FG_get_graph(graph);
	vsize_t k = INTEGER(pk)[0];
	std::string init = CHAR(STRING_ELT(pinit,0));
	vsize_t max_iters = INTEGER(pmax_iters)[0];
	double tolerance = REAL(ptolerance)[0];

    Rcpp::List ret;
    sem_kmeans_ret::ptr fg_ret = compute_sem_kmeans(fg, k, init, max_iters, tolerance);

	fm::vector::ptr clusters = fg_ret->get_cluster_assignments();
    ret["cluster"] = create_FMR_vector(get_vertex_ids(clusters), "");
    ret["iter"] = fg_ret->get_iters();

    Rcpp::IntegerVector res1(fg_ret->get_size().begin(), fg_ret->get_size().end());
    ret["size"] = res1;

    const unsigned NUM_COLS = fg_ret->get_centers()[0].size();
	Rcpp::NumericMatrix centers = Rcpp::NumericMatrix(k, NUM_COLS);
#pragma omp parallel for firstprivate(fg_ret) shared(centers)
	for (unsigned row = 0; row < k; row++) {
		for (unsigned col = 0; col < NUM_COLS; col++) {
			centers(row, col) =  fg_ret->get_centers()[row][col];
		}
	}
    ret["centers"] = centers;
	return ret;
}

RcppExport SEXP R_FG_compute_betweenness(SEXP graph, SEXP _vids)
{
	Rcpp::IntegerVector Rvids(_vids);
	std::vector<vertex_id_t> vids(Rvids.begin(), Rvids.end());
	FG_graph::ptr fg = R_FG_get_graph(graph);

	fm::vector::ptr fg_vec = compute_betweenness_centrality(fg, vids);
	return create_FMR_vector(cast_type<double>(fg_vec), "");
}

SEXP create_FMR_matrix(fm::sparse_matrix::ptr m, R_type type, const std::string &name);

namespace fg
{
fm::sparse_matrix::ptr create_sparse_matrix(FG_graph::ptr fg,
		const fm::scalar_type *entry_type);
}

RcppExport SEXP R_FG_get_matrix_fg(SEXP pgraph)
{
	Rcpp::List graph = Rcpp::List(pgraph);
	Rcpp::LogicalVector res(1);
	fg::FG_graph::ptr fg = R_FG_get_graph(pgraph);
	// TODO does this work if this isn't a binary matrix?
	fm::sparse_matrix::ptr m = fg::create_sparse_matrix(fg, NULL);
	std::string name = graph["name"];
	// TODO change it later for non-binary matrix.
	return create_FMR_matrix(m, R_type::R_LOGICAL, name);
}

RcppExport SEXP R_FG_print_graph(SEXP pgraph, SEXP pfile, SEXP pdelim,
		SEXP ptype)
{
	fg::FG_graph::ptr fg = R_FG_get_graph(pgraph);
	std::string file_name = CHAR(STRING_ELT(pfile, 0));
	std::string delim = CHAR(STRING_ELT(pdelim, 0));
	std::string type = CHAR(STRING_ELT(ptype, 0));
	Rcpp::LogicalVector res(1);

	FILE *f = fopen(file_name.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		res[0] = false;
		return res;
	}
	print_graph_el(fg, delim, type, f);
	fclose(f);
	res[0] = true;
	return res;
}
