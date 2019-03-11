/*
 * Copyright 2019
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

#ifndef __FG_CGRAPH__H__
#define __FG_CGRAPH__H__

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


namespace fg {

class CGraph {
    private:
        bool isC;
        std::string graph_file;
        std::string index_file;
        config_map::ptr configs;
        FG_graph::ptr graph;

    public:
        CGraph() { }
        CGraph(const std::string& graph_file, const std::string& index_file,
                const std::string& config_file) : graph_file(graph_file),
                index_file(index_file) {

                    config_map::ptr configs = config_map::create(config_file);
                    if (configs == NULL)
                        configs = config_map::ptr();

                    graph_engine::init_flash_graph(configs);
                    graph = FG_graph::create(graph_file, index_file, configs);
                }

        // Utils
        const vertex_id_t vcount() {
            return graph->get_graph_header().get_num_vertices();
        }

        // Coreness
        std::vector<size_t> coreness(const size_t kmax=0, const size_t kmin=0) {
            return compute_kcore(graph, kmin, kmax, true)->conv2std<size_t>();
        }

        // Betweenness Centrality
        std::vector<float> betweenness(std::vector<vertex_id_t>& ids) {

            if (ids.empty()) {
                for (vertex_id_t id = 0; id < vcount(); id++) {
                    ids.push_back(id);
                }
            }

            return compute_betweenness_centrality(graph, ids)->conv2std<float>();
        }

        // Degree
        std::vector<vertex_id_t> degree(const std::string& etype="both") {
            edge_type type;
            if (etype == "in")
                type = IN_EDGE;
            else if (etype == "out")
                type = OUT_EDGE;
            else
                type = BOTH_EDGES;

            return get_degree(graph, type)->conv2std<vertex_id_t>();
        }

        // Triangle counting
        std::vector<vertex_id_t> triangles(bool cycles_only=false) {
            if (graph->is_directed()) {
                if (cycles_only)
                    return compute_directed_triangles_fast(graph,
                            directed_triangle_type::CYCLE)->
                                conv2std<vertex_id_t>();
                else
                    return compute_directed_triangles(graph,
                            directed_triangle_type::CYCLE)->
                                conv2std<vertex_id_t>();
            } else {
                return compute_undirected_triangles(graph)->
                            conv2std<vertex_id_t>();
            }
        }

        // Local scan
        std::vector<size_t> local_scan(const short num_hops=1) {
            if (num_hops == 1)
                return compute_local_scan(graph)->conv2std<size_t>();
            else if (num_hops == 2)
                return compute_local_scan2(graph)->conv2std<size_t>();
            else
                throw std::runtime_error(
                        "local_scan: num_hops must be 1 or 2\n");
        }

        // Top k scan
        typedef std::pair<vertex_id_t, size_t> tk_t;
        std::vector<tk_t> topk_scan(const size_t k) {

            FG_vector<tk_t>::ptr scan = compute_topK_scan(graph, k);

            std::vector<tk_t> ret(scan->get_size());
            scan->copy_to<tk_t>(&ret[0], ret.size());
            return ret;
        }

        // Diameter
        size_t diameter(const unsigned num_para_bfs=2, bool directed=true) {
	        return estimate_diameter(graph, num_para_bfs, directed);
        }

        // Pagerank
        std::vector<float> pagerank(const int num_iters=30,
                const float damping_factor=0.85,
                const std::string& algo="push") {

            if (algo == "pull") {
                return compute_pagerank(graph, num_iters,
                        damping_factor)->conv2std<float>();
            } else { /* Default */
                return compute_pagerank2(graph, num_iters,
                        damping_factor)->conv2std<float>();
            }
        }

        // Weakly connected componentes
        std::vector<vertex_id_t> weakly_connected_components(
                const bool sync=false) {
            if (sync) {
                return compute_sync_wcc(graph)->conv2std<vertex_id_t>();
            } else {
                return compute_wcc(graph)->conv2std<vertex_id_t>();
            }
        }

        // Connected Components
        std::vector<vertex_id_t> connected_components() {
            return compute_cc(graph)->conv2std<vertex_id_t>();
        }

        // Strongly Connected Components
        std::vector<vertex_id_t> strongly_connected_components() {
            return compute_scc(graph)->conv2std<vertex_id_t>();
        }

        // bfs
        size_t bfs_vcount(vertex_id_t start_vertex=INVALID_VERTEX_ID,
                const std::string edge_type_str="both") {
            auto edge = edge_type::BOTH_EDGES;

            if (edge_type_str == "in")
                edge = edge_type::IN_EDGE;
            else if (edge_type_str == "out")
                edge = edge_type::OUT_EDGE;

            if (start_vertex == INVALID_VERTEX_ID ||
                    start_vertex >= vcount())
                start_vertex = random() % vcount() - 1;

            size_t bfs(FG_graph::ptr fg, vertex_id_t start_vertex, edge_type);
            return bfs(graph, start_vertex, edge);
        }

        ~CGraph() {
            graph_engine::destroy_flash_graph();
        }
};
}

#endif
