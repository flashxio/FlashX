/*
 * Copyright 2019
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

namespace fg {

class CGraph {
    private:
        std::string graph_file;
        std::string index_file;
        config_map::ptr configs;
        FG_graph::ptr graph;
        bool min_vertex_id_set;
        vertex_id_t min_vertex_id;

    public:
        CGraph() { }
        CGraph(std::string graph_file, std::string index_file,
                std::string config_file) : graph_file(graph_file),
                index_file(index_file), min_vertex_id_set(false),
                min_vertex_id(0) {

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

        const size_t ecount() const {
            return graph->get_num_edges();
        }

        const bool is_directed() const {
            return graph->is_directed();
        }

        const bool is_in_mem() const {
            return graph->is_in_mem();
        }
        // End Utils

        // Begin Algs
        // Coreness
        std::vector<size_t> coreness(const size_t kmax=0, const size_t kmin=0) {
            return compute_kcore(graph, kmin, kmax, true);
        }

        // Betweenness Centrality
        std::vector<float> betweenness(std::vector<vertex_id_t>& ids) {

            if (ids.empty()) {
                for (vertex_id_t id = 0; id < vcount(); id++) {
                    ids.push_back(id);
                }
            }

            return compute_betweenness_centrality(graph, ids);
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

            return get_degree(graph, type);
        }

        // Triangle counting
        std::vector<size_t> triangles(bool cycles_only=false) {
            if (graph->is_directed()) {
                if (cycles_only)
                    return compute_directed_triangles_fast(graph,
                            directed_triangle_type::CYCLE);
                else
                    return compute_directed_triangles(graph,
                            directed_triangle_type::CYCLE);
            } else {
                return compute_undirected_triangles(graph);
            }
        }

        // Local scan
        std::vector<size_t> local_scan(const short num_hops=1) {
            if (num_hops == 1)
                return compute_local_scan(graph);
            else if (num_hops == 2)
                return compute_local_scan2(graph);
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
                        damping_factor);
            } else { /* Default */
                return compute_pagerank2(graph, num_iters,
                        damping_factor);
            }
        }

        // Weakly connected componentes
        std::vector<vertex_id_t> weakly_connected_components(
                const bool sync=false) {
            if (sync) {
                return compute_sync_wcc(graph);
            } else {
                return compute_wcc(graph);
            }
        }

        // Connected Components
        std::vector<vertex_id_t> connected_components() {
            return compute_cc(graph);
        }

        // Strongly Connected Components
        std::vector<vertex_id_t> strongly_connected_components() {
            return compute_scc(graph);
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

            return bfs(graph, start_vertex, edge);
        }

        std::string to_str() {
			return std::string("Graphyti Graph:\n") +
				std::string("v: ") + std::to_string(vcount()) +
                std::string(", e: ") + std::to_string(ecount()) +
                (is_directed() ? std::string("\nDirected"):
                std::string("\nUndirected")) +
                (is_in_mem() ? std::string("\nIn-memory\n"):
                 std::string("\nOn-Disk\n"));
        }

        const vertex_id_t min_id() {
            if (!min_vertex_id_set) {
                class qvertex : public compute_vertex {
                    public:
                    qvertex(vertex_id_t id): compute_vertex(id) { }
                    void run(vertex_program&) {}
                    void run(vertex_program&, const page_vertex&) {}
                    void run_on_message(vertex_program&, const vertex_message&) {}
                };

                graph_index::ptr index = NUMA_graph_index<qvertex>::create(
                        graph->get_graph_header());
                min_vertex_id = index->get_min_vertex_id();
                min_vertex_id_set = true; // set it forever
            }
            return min_vertex_id;
        }

        const vertex_id_t max_id() {
            return graph->get_num_vertices() - 1;
        }

        // TODO: tocsr
        // TODO: toigraph

        ~CGraph() {
            graph_engine::destroy_flash_graph();
        }
        // End Algs
};
}
#endif
