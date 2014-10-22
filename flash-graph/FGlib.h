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
#include "graph.h"
#include "FG_vector.h"
#include "graph_file_header.h"

/**
  * \brief A user-friendly wrapper for FlashGraph's raw graph type.
  *         Very usefule when when utilizing FlashGraph 
  *         pre-written/library algorithms.
  *
*/
class FG_graph
{
	std::string graph_file;
	std::string index_file;
	std::shared_ptr<in_mem_graph> graph_data;
	std::shared_ptr<vertex_index> index_data;
	config_map::ptr configs;

	FG_graph(const std::string &graph_file,
			const std::string &index_file, config_map::ptr configs);
public:
	typedef std::shared_ptr<FG_graph> ptr; /**Smart pointer through which object is accessed*/

	~FG_graph() {
		graph_engine::destroy_flash_graph();
	}

/**
  * \brief  Method to instantiate a graph object.
  *         This method is used in lieu of explicitly calling a ctor.
  *    
  * \param graph_file Path to the graph file on disk.
  * \param index_file Path to the graph index file on disk.
  * \param configs Configuration in configuration file.
*/
	static ptr create(const std::string &graph_file,
			const std::string &index_file, config_map::ptr configs) {
		return ptr(new FG_graph(graph_file, index_file, configs));
	}

/**
  * \brief Get the path to the graph file.
  *
  * \return The path to the graph file on disk.
  *
*/
	const std::string &get_graph_file() const {
		return graph_file;
	}

/**
  * \brief Get the graph index file path.
  *
  * \return The path to the graph index file on disk.
*/
	const std::string &get_index_file() const {
		return index_file;
	}

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

	std::shared_ptr<in_mem_graph> get_graph_data() const {
		return graph_data;
	}

	std::shared_ptr<vertex_index> get_index_data() const {
		return index_data;
	}

	graph_engine::ptr create_engine(graph_index::ptr index);
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

/**
  * \brief Compute all weakly connectected components of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
  *
*/
FG_vector<vertex_id_t>::ptr compute_wcc(FG_graph::ptr fg);

/**
  * \brief Compute all weakly connectected components of a graph synchronously.
  * The reason of having this implementation is to understand the performance
  * of synchronous wcc and asynchronous wcc.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
 */
FG_vector<vertex_id_t>::ptr compute_sync_wcc(FG_graph::ptr fg);

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
FG_vector<vertex_id_t>::ptr compute_ts_wcc(FG_graph::ptr fg,
		time_t start_time, time_t time_interval);

/**
  * \brief Compute all strongly connected components of a graph.
  *
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with a component ID for each vertex in the graph.
  *
*/
FG_vector<vertex_id_t>::ptr compute_scc(FG_graph::ptr fg);

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
FG_vector<size_t>::ptr compute_directed_triangles(FG_graph::ptr fg,
		directed_triangle_type type);
FG_vector<size_t>::ptr compute_directed_triangles_fast(FG_graph::ptr fg,
		directed_triangle_type type);
/**
  * \brief Compute undirected triangle count for each vertex.
  * 
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector that contains the number of triangles associated with
  *         each vertex in the graph.
  *
*/
FG_vector<size_t>::ptr compute_undirected_triangles(FG_graph::ptr fg);

/**
  * \brief Compute the per-vertex local Scan Statistic 
  * \param fg The FlashGraph graph object for which you want to compute.
  * \return A vector with an entry for each vertex in the graph's 
  *          local scan value.
  *
*/
FG_vector<size_t>::ptr compute_local_scan(FG_graph::ptr);

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
FG_vector<float>::ptr compute_pagerank(FG_graph::ptr fg, int num_iters,
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
FG_vector<float>::ptr compute_pagerank2(FG_graph::ptr, int num_iters,
		float damping_factor);

FG_vector<float>::ptr compute_sstsg(FG_graph::ptr fg, time_t start_time,
		time_t interval, int num_intervals);

/**
 * \brief Fetch the clusters with the wanted cluster IDs.
 *  
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param vertices The vertices that the induced subgraph has.
 * \return A subgraph.
 */
in_mem_subgraph::ptr fetch_subgraph(FG_graph::ptr graph,
		const std::vector<vertex_id_t> &vertices);

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
FG_vector<size_t>::ptr compute_kcore(FG_graph::ptr fg,
		                size_t k, size_t kmax=0);

/**
 * \brief Get the degree of all vertices in the graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param type The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
 * \return A vector with an entry for each vertex degree.
 */
FG_vector<vsize_t>::ptr get_degree(FG_graph::ptr fg, edge_type type);

/**
 * \brief Get the degree of all vertices in a specified time interval in
 *        a time-series graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \param type The edge type: IN_EDGE, OUT_EDGE, BOTH_EDGES.
 * \param start_time The start time of the time interval.
 * \param time_interval length of the time interval.
 * \return A vector with an entry for each vertex degree.
 */
FG_vector<vsize_t>::ptr get_ts_degree(FG_graph::ptr fg, edge_type type,
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
 * \brief Compute transitivity of all vertices in the graph.
 * \param fg The FlashGraph graph object for which you want to compute.
 * \return A vector with an transitivity value for each vertex.
 */
FG_vector<float>::ptr compute_transitivity(FG_graph::ptr fg);

/**
 * \brief Get the header of the graph that contains basic information of the graph.
 * \return The graph header.
 */
graph_header get_graph_header(FG_graph::ptr fg);

#endif
