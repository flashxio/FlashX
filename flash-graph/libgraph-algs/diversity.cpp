/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY CURRENT_KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <signal.h>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <cmath>
#include "graph_engine.h"
#include "FGlib.h"

using namespace fg;

namespace {

edge_type g_edge_type = edge_type::OUT_EDGE;
typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
std::vector<float> diversity_v;
bool g_memopt = false;

class diversity_vertex: public compute_directed_vertex
{
	public:
	diversity_vertex(vertex_id_t id): compute_directed_vertex(id) {
	}

	void run(vertex_program &prog) {
        directed_vertex_request req(prog.get_vertex_id(*this), g_edge_type);
        request_partial_vertices(&req, 1);
    }

	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg) { };
};

void diversity_vertex::run(vertex_program &prog, const page_vertex &vertex) {
    // Get the total_degree k
    vsize_t id = prog.get_vertex_id(*this);
    vsize_t tot_deg = prog.get_graph().get_num_edges(id, g_edge_type);

    if (tot_deg < 2)
        return; // We leave div as 0

    // Get the weights: w
    data_seq_iterator weight_it =
        ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(
            g_edge_type);

    float H = 0;

    if (!g_memopt) {
        std::vector<float> p;
        size_t weight_sum = 0;

        while (weight_it.has_next()) {
            edge_count e = weight_it.next();
            weight_sum += e.get_count();
            p.push_back(e.get_count());
        }

        // Normalize p and compute H and D
        for (auto& _ : p) {
            float tmp = _ / (float)weight_sum;
            H += (tmp * std::log2(tmp));
        }

    } else {// Memory Saver
        size_t weight_sum = 0;
        while (weight_it.has_next()) {
            edge_count e = weight_it.next();
            weight_sum += e.get_count();
        }

        weight_it =
            ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(
                g_edge_type);

        while (weight_it.has_next()) {
            edge_count e = weight_it.next();
            float p_ij = e.get_count()/(float)weight_sum;
            H += p_ij * std::log2(p_ij);
        }
    }
    diversity_v[id] = (-H) / std::log2(tot_deg);
}
}


namespace fg
{

std::vector<float> compute_diversity(FG_graph::ptr fg, edge_type traverse_e,
        bool memopt)
{
	graph_index::ptr index = NUMA_graph_index<diversity_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);
    g_edge_type = traverse_e;
    g_memopt = memopt;

    if (g_edge_type == BOTH_EDGES)
        throw std::runtime_error("Only IN and OUT edges usable\n");

	struct timeval start, end;
	gettimeofday(&start, NULL);

    // Run the computation
    printf("Starting diversity ...\n");
    diversity_v.assign(fg->get_num_vertices(), 0.0);
    graph->start_all();
    graph->wait4complete();

	gettimeofday(&end, NULL);
    printf("diversity took %.5f sec to complete\n", time_diff(start, end));

#if 1
    printf("diversity values: \n");
    printf("[ ");
    for (auto const& v : diversity_v)
        printf("%.3f ", v);
    printf("]\n");
#endif

	return diversity_v;
}
}
