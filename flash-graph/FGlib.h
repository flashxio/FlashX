#ifndef __FGLIB_H__
#define __FGLIB_H__

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

#include "graph_engine.h"
#include "FG_vector.h"
#include "graph_file_header.h"

namespace safs
{
	class file_io_factory;
};

namespace fg
{

/**
  * \brief A user-friendly wrapper for FlashGraph's raw graph type.
  *         Very usefule when when utilizing FlashGraph
  *         pre-written/library algorithms.
  *
*/
class FG_graph
{
	graph_header header;
	std::string graph_file;
	std::string index_file;
	std::shared_ptr<in_mem_graph> graph_data;
	std::shared_ptr<vertex_index> index_data;
	config_map::ptr configs;

	// In this case, the graph file is kept in SAFS and the index is read to
	// memory.
	FG_graph(const std::string &graph_file,
			std::shared_ptr<vertex_index> index_data, config_map::ptr configs);
	// In this case, both the graph and the index are read to memory.
	FG_graph(std::shared_ptr<in_mem_graph> graph_data,
			std::shared_ptr<vertex_index> index_data,
			const std::string &graph_name, config_map::ptr configs);
public:
	typedef std::shared_ptr<FG_graph> ptr; /**Smart pointer through which object is accessed*/

	~FG_graph() {
	}

	/**
	 * \brief  Method to instantiate a graph object.
	 *         This method is used in lieu of explicitly calling a ctor.
	 *
	 * \param graph_file Path to the graph file in SAFS or in Linux filesystem.
	 * \param index_file Path to the graph index file in SAFS or
	 *        in Linux filesystem.
	 * \param configs Configuration in configuration file.
	 */
	static ptr create(const std::string &graph_file,
			const std::string &index_file, config_map::ptr configs);

	/**
	 * \brief  Method to instantiate a graph object.
	 *         This method is used in lieu of explicitly calling a ctor.
	 *
	 * \param graph_data The adjacency lists of the graph stored in memory.
	 * \param index_data The index of the graph stored in memory.
	 * \param graph_name The name of the graph.
	 * \param configs Configuration in configuration file.
	 */
	static ptr create(std::shared_ptr<in_mem_graph> graph_data,
			std::shared_ptr<vertex_index> index_data,
			const std::string &graph_name, config_map::ptr configs) {
		return ptr(new FG_graph(graph_data, index_data, graph_name, configs));
	}

	std::shared_ptr<safs::file_io_factory> get_graph_io_factory(int access_option);

/**
  * \brief Get the map that contains the runtime configurations
  *        for FlashGraph.
  *
  * \return The config_map that contains all FlashGraph configurations.
  *
*/
	config_map::ptr get_configs() const {
		return configs;
	}

	bool is_in_mem() const {
		return graph_data != NULL;
	}

	std::shared_ptr<in_mem_graph> get_graph_data() const;
	std::shared_ptr<vertex_index> get_index_data() const;

	graph_engine::ptr create_engine(graph_index::ptr index);

	/**
	 * \brief Get the header of the graph that contains basic
     *      information of the graph.
	 * \return The graph header.
	 */
	const graph_header &get_graph_header() const {
		return header;
	}

	bool is_directed() const {
		return get_graph_header().is_directed_graph();
	}

	size_t get_num_vertices() const {
		return get_graph_header().get_num_vertices();
	}

	size_t get_num_edges() const {
		return get_graph_header().get_num_edges();
	}
};

/**
  * \brief Triangle computation type.
  *
  * - CYCLE triangles are defined for directed graphs and
  *         depend on the direction of each edge. All edges must be
  *         head to tail connections.
  *     E.g A -----> B
  *         ^     /
  *         |   /
  *         | v
  *         C
  *
  * - ALL triangles. Edge direction is disregarded.
  *     E.g A ----- B
  *         |     /
  *         |   /
  *         | /
  *         C
  */
enum directed_triangle_type
{
	CYCLE,
	ALL,
};

std::vector<vertex_id_t> compute_cc(FG_graph::ptr fg);

/**
  * \brief Compute all weakly connectected components of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
  *
*/
std::vector<vertex_id_t> compute_wcc(FG_graph::ptr fg);

/**
  * \brief Compute all weakly connectected components of a graph synchronously.
  * The reason of having this implementation is to understand the performance
  * of synchronous wcc and asynchronous wcc.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
 */
std::vector<vertex_id_t> compute_sync_wcc(FG_graph::ptr fg);

/**
 * \brief Compute all weakly connectected components of a time-series graph
 *        in a specified time interval.
 *
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param start_time The start time of the time interval.
 * \param time_interval The length of the time interval.
 * \return A vector with a component ID for each vertex in the graph.
 *
 */
std::vector<vertex_id_t> compute_ts_wcc(FG_graph::ptr fg,
		time_t start_time, time_t time_interval);

/**
  * \brief Compute all strongly connected components of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
  *
*/
std::vector<vertex_id_t>compute_scc(FG_graph::ptr fg);

/**
  * \brief Compute the directed triangle count for each each vertex.
  *        Currently, it only counts the number of cycle triangles.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param type The type of triangles you wish to count.
  * \return A vector that contains the number of triangles associated with
  *         each vertex in the graph.
  *
*/
std::vector<size_t> compute_directed_triangles(FG_graph::ptr fg,
		directed_triangle_type type);
std::vector<size_t> compute_directed_triangles_fast(FG_graph::ptr fg,
		directed_triangle_type type);
/**
  * \brief Compute undirected triangle count for each vertex.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector that contains the number of triangles associated with
  *         each vertex in the graph.
  *
*/
std::vector<size_t> compute_undirected_triangles(FG_graph::ptr fg);

/**
  * \brief Compute the per-vertex local Scan Statistic
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with an entry for each vertex in the graph's
  *          local scan value.
  *
*/
std::vector<size_t> compute_local_scan(FG_graph::ptr);
std::vector<size_t> compute_local_scan2(FG_graph::ptr fg);

/**
  * \brief Obtain the top K vertices with the largest local Scan
  *     Statistic value.
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param topK The value for K used for the `top K` vertices.
  * \return A vector of std::pair with an entry for each vertex
  *         in the top K together with its value.
  *
*/
FG_vector<std::pair<vertex_id_t, size_t> >::ptr compute_topK_scan(
		FG_graph::ptr, size_t topK);

/**
  * \brief Compute the diameter estimation for a graph.
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return The diameter estimate value.
  *
*/
size_t estimate_diameter(FG_graph::ptr fg, int num_bfs, bool directed);

/**
  * \brief Compute the PageRank of a graph using the pull method
  *       where vertices request the data from all their neighbors
  *       each iteration. Tends to converge to stable values.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param num_iters The maximum number of iterations for PageRank.
  * \param damping_factor The damping factor. Originally .85.
  *
  * \return A vector with an entry for each vertex in the graph's
  *         PageRank value.
  *
*/
std::vector<float> compute_pagerank(FG_graph::ptr fg, int num_iters,
		float damping_factor);

/**
  * \brief Compute the PageRank of a graph using the push method
  *       where vertices send deltas of their PageRank to neighbors
  *       in the event their own PageRank changes.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param num_iters The maximum number of iterations for PageRank.
  * \param damping_factor The damping factor. Originally .85.
  *
  * \return A vector with an entry for each vertex in the graph's
  *         PageRank value.
  *
*/
std::vector<float> compute_pagerank2(FG_graph::ptr, int num_iters,
		float damping_factor);

std::vector<float> compute_sstsg(FG_graph::ptr fg, time_t start_time,
		time_t interval, int num_intervals);

/**
 * \brief Compute the k-core/coreness of a graph. The algorithm will
 *        determine which vertices are between core `k` and `kmax` --
 *        all other vertices will be assigned to core 0.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param k The core value to be computed.
 * \param kmax (Optional) The kmax value. If omitted then all cores are
 *        computed i.e., coreness. *This is not recommended for very large
 *        graphs.*
 * \return An `FG_vector` containing the core of each vertex between `k`
 *         and `kmax`. All other vertices are assigned to core 0.
 */
std::vector<size_t> compute_kcore(FG_graph::ptr fg, size_t k, size_t kmax=0,
        bool skip=true);

/**
 * \brief Get the degree of all vertices in the graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param type The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
 * \return A vector with an entry for each vertex degree.
 */
std::vector<vertex_id_t> get_degree(FG_graph::ptr fg, edge_type type);

/**
  * \brief Compute the transitivity/clustering coefficient of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with an entry for each vertex in the graph's
  *         transitivity value.
*/
// TODO
std::vector<float> compute_transitivity(FG_graph::ptr fg);// {
    //throw std::runtime_error("Transitivity not yet implemented");
//}

/**
  * \brief Compute the betweenness centrality of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param vids The vertex IDs for which BC should be computed
  * \return A vector with an entry for each vertex in the graph's
  *         betweenness centrality value.
*/
std::vector<float> compute_betweenness_centrality(FG_graph::ptr fg,
		const std::vector<vertex_id_t>& vids);

/**
 * \brief Get the degree of all vertices in a specified time interval in
 *        a time-series graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param type The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
 * \param start_time The start time of the time interval.
 * \param time_interval length of the time interval.
 * \return A vector with an entry for each vertex degree.
 */
std::vector<vertex_id_t> get_ts_degree(FG_graph::ptr fg, edge_type type,
		time_t start_time, time_t time_interval);

/**
 * \brief Get the time range in which the time-series graph is.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \return A pair of timestamp that defines the time range of the time-series
 * graph.
 */
std::pair<time_t, time_t> get_time_range(FG_graph::ptr fg);

/**
 * \brief Get the neighborhood overlap of each pair of vertices in `vids'.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param vids The vertices whose neighborhood overlap is computed.
 * \param overlap_matrix A dense matrix that stores the overlap of each pair of vertices.
 */
void compute_overlap(FG_graph::ptr fg, const std::vector<vertex_id_t> &vids,
		std::vector<std::vector<double> > &overlap_matrix);

/**
 * \brief Compute louvain clustering for a graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param levels The number of levels of the hierarchy to do.
 */
std::vector<unsigned> compute_louvain(FG_graph::ptr fg, const uint32_t levels);

/**
  * \brief Run breadth first search and return the number of vertices visited.
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param start_vertex The source vertex for the BFS
  * \param type The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
  *
  */
size_t bfs(FG_graph::ptr fg, vertex_id_t start_vertex,
        edge_type type=BOTH_EDGES);

/**
  * \brief Compute the closeness centrality of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param ids The vertex IDs for which closeness should be computed
  * \param traverse_e The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
  * \return A vector with an entry for each vertex in the graph's
  *         closeness centrality value.
*/
std::vector<double> compute_closeness_centrality(FG_graph::ptr fg,
        std::vector<vsize_t>& ids, edge_type traverse_e=BOTH_EDGES);

/**
  * \brief The diversity of a vertex is defined as the (scaled)
  *     Shannon entropy of the weights of its incident edges.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \param traverse_e The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
  * \param memopt Optimize for minimal memory rather than performance.
  * \return A vector with an entry for each vertex in the graph's
  *         diversity value.
*/
std::vector<float> compute_diversity(FG_graph::ptr fg,
        edge_type traverse_e=OUT_EDGE, bool memopt=false);

/**
  * \brief Compute one (of possibly many) topological sortings of the graph.
  *
  * \return A vector with a topological sorting of the vertices.
*/
std::vector<vertex_id_t> compute_topo_sort(FG_graph::ptr fg, bool approx=true);

}
#endif
