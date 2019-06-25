/**
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <limits>
#include <atomic>

#include "FGlib.h"
#include "save_result.h"
#include "helper.hpp"

using namespace fg;

namespace {
	// Global map from cluster_id : cluster(volume, weight)
	//std::vector<atomicwrapper<uint32_t>> g_weight_vec;
	//std::vector<atomicwrapper<uint32_t>> g_volume_vec;

	enum stage_t
	{
		INIT, /* Compute m, ki, ki_in*/
		PC_RUN, /* Compute sigma_in, sigma_tot */
		RUN, /* Compute moudularity diffrerence */
		REBUILD, /* Move nodes to clusters */
	};

	typedef safs::page_byte_array::seq_const_iterator
                        <edge_count> data_seq_iterator;
	typedef vertex_id_t cluster_id_t;

    // Globals
	stage_t stage;
    std::vector<cluster_id_t> membership;
    barrier::ptr stage_barrier;
	size_t m = 0;
	bool g_changed = false; // If no change then end
    unsigned max_levels = 0;
    unsigned cur_level = 0;

	// If anything changes cluster we cannot converge
	void set_changed(bool changed) {
		if (!changed) // Only once per iter so no prob with this
			g_changed = changed;
		if (changed && (!g_changed))
			g_changed = changed;
	}

    // This message is sent from vertex a -> b when a wants to merge with b
	class merge_message: public vertex_message
	{
        vertex_id_t sender_id; // will include `a` too

		public:
		merge_message(vertex_id_t& sender_id):
			vertex_message(sizeof(merge_message), true), sender_id(sender_id) {
			}

		const vertex_id_t& get_sender_id() const {
			return sender_id;
		}
	};

	class louvain_vertex: public compute_vertex
	{
		cluster_id_t next_cluster_id; // current cluster
		size_t ki; // sum of weights incident to this vertex
        size_t ki_in; // -- computed per cluster // TODO: Reset at iteration end
        std::vector<vertex_id_t> members; // Worst case O(n) total

		public:
		louvain_vertex(vertex_id_t id): compute_vertex(id),
                next_cluster_id(id), ki(0) {
                        members.push_back(id);
		}

		void run(vertex_program &prog) {
            vertex_id_t id = prog.get_vertex_id(*this);
            switch (stage) {
                case INIT:
                    // Compute global edge weight
                    request_vertices(&id, 1);
                    break;
                case RUN:
                    assert(membership.size() > id);
                    if (membership[id] == id) {
                        request_vertices(&members[0], members.size());
                    } // else
                        // I've been merged into another vertex
                    break;
                default:
                    throw std::runtime_error("Unknown louvain stage");
			}
		}

		void run(vertex_program &prog, const page_vertex &vertex);
		void run_on_message(vertex_program &prog, const vertex_message &msg1);

		void compute_modularity(cluster_id_t neigh_cluster_id,
                vertex_id_t my_id, vertex_program& prog);

		void compute_statics(data_seq_iterator& weight_it,
                edge_seq_iterator& id_it, vertex_program &prog);

		//// Remove vertex from old cluster
		//void switch_clusters () {
			//std::cout << " Switching! Moving from c" <<
                //cluster_id << " to c" << next_cluster_id << std::endl;
			//g_weight_vec[cluster_id].minus_eq(weight);
			//g_volume_vec[cluster_id].minus_eq(volume);

			//cluster_id = next_cluster_id;
			//g_weight_vec[cluster_id].plus_eq(weight);
			//g_volume_vec[cluster_id].plus_eq(volume);
		//}

		//void notify_iteration_end(vertex_program &vprog);

	};

		void louvain_vertex::run_on_message(vertex_program &prog,
                const vertex_message &msg1) {
            // TODO: Notify neighbor that you're joining their cluster
        }

	void louvain_vertex::run(vertex_program &prog, const page_vertex &vertex) {
		switch (stage) {
			case INIT: /* INIT just accums the global edge_count */
				{
					// Out edges
                    if (vertex.get_num_edges(BOTH_EDGES) == 0) return;

                    // IN
                    auto num_dests = vertex.get_num_edges(IN_EDGE);
					data_seq_iterator weight_it =
						((const page_directed_vertex&)vertex).
                            get_data_seq_it<edge_count>(IN_EDGE, 0, num_dests);
					edge_seq_iterator id_it = vertex.get_neigh_seq_it(IN_EDGE,
                            0, num_dests);
					compute_statics(weight_it, id_it, prog);

                    // OUT
                    num_dests = vertex.get_num_edges(OUT_EDGE);
					weight_it =
						((const page_directed_vertex&)vertex).
                            get_data_seq_it<edge_count>(OUT_EDGE, 0, num_dests);
					id_it = vertex.get_neigh_seq_it(OUT_EDGE,
                            0, num_dests);
					compute_statics(weight_it, id_it, prog);
				}
				break;
			case RUN:
				{
                    // Maximum modularity value and cluster id
                    float max_mod = std::numeric_limits<float>::min();
                    cluster_id_t cid = INVALID_VERTEX_ID;

                    // TODO

				}
				break;
			default:
				assert(0);
		}
	}

	// Iterate through neighbors to see if we switch membership, instead of
    //  the min computation which is every cluster.
	void louvain_vertex::compute_modularity(cluster_id_t neigh_cluster_id,
            vertex_id_t my_id, vertex_program& prog) {

#if 0
		// TODO: Remove divisions & hope for no overflow
		float delta_mod = ((int)(g_weight_vec[neigh_cluster_id].get() -
                (g_weight_vec[this->next_cluster_id].get() - this->weight))
                / (float) m) +
			(((int)((g_volume_vec[this->next_cluster_id].get() - this->volume)
                    - g_volume_vec[neigh_cluster_id].get())
              * (int)this->volume) / (float)(2*(m*m)));

		if (delta_mod > this->max_modularity) {
			/** DEBUG **/
			std::cout << "v" << my_id << " tentative move from c" <<
                cluster_id << " to c" << neigh_cluster_id << ", delta_mod=" <<
                delta_mod << " > max_mod=" << max_modularity << "\n";
			/** GUBED **/

			max_modularity = delta_mod;
			this->next_cluster_id = neigh_cluster_id;
			//std::cout << "Activating v" << my_id << " for the switch iteration";
			set_changed(true); // Global notification of cluster change

			//prog.request_notify_iter_end(*this);
            // Activate myself so I can change my clusters metadata
			//prog.activate_vertices(&my_id,1);
		}
#endif
	}

#if 0
	void louvain_vertex::notify_iteration_end(vertex_program &prog) {
		assert (next_cluster_id != cluster_id);
		vertex_id_t my_id = prog.get_vertex_id(*this);
        // Activate myself so I can change my clusters metadata in NEXT iteration
		prog.activate_vertices(&my_id,1);
	}

#endif

	class louvain_vertex_program: public vertex_program_impl<louvain_vertex>
	{
		size_t local_m; // For the global edge count
        // TODO: Local clusters

		public:
		louvain_vertex_program() : local_m(0) {
		}

		typedef std::shared_ptr<louvain_vertex_program> ptr;

		static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<louvain_vertex_program,
                   vertex_program>(prog);
		}

		void pp_m(const size_t weight) {
			this->local_m += weight;
		}

		const size_t get_local_m() const {
			return local_m;
		}
	};

	class louvain_vertex_program_creater: public vertex_program_creater
	{
		public:
			vertex_program::ptr create() const {
				return vertex_program::ptr(new louvain_vertex_program());
			}
	};

	// Only need to do this once per vertex ever
	void louvain_vertex::compute_statics(
            data_seq_iterator& weight_it, edge_seq_iterator& id_it,
			vertex_program &prog) {

        auto vid = prog.get_vertex_id(*this);

		while (weight_it.has_next()) {
			edge_count e = weight_it.next();
			vertex_id_t nid = id_it.next();

            auto cnt = e.get_count();
			ki += cnt;
			if (nid == vid)
				ki_in += cnt;
		}

		((louvain_vertex_program&)prog).pp_m(ki);
	}
}

namespace fg
{
    std::vector<unsigned> compute_louvain(FG_graph::ptr fg,
            const unsigned levels) {

        if (!fg->get_graph_header().is_directed_graph()) {
            throw std::runtime_error(
                    "This algorithm works on a directed graph\n");
        }

		graph_index::ptr index = NUMA_graph_index<louvain_vertex>::create(
				fg->get_graph_header());
		graph_engine::ptr graph = fg->create_engine(index);

		printf("Starting Louvain with %u\n ..", levels);
        printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
#ifdef PROFILER
		if (!graph_conf.get_prof_file().empty())
			ProfilerStart(graph_conf.get_prof_file().c_str());
#endif
		struct timeval start, end;
		gettimeofday(&start, NULL);

        // Initialize state
        stage = INIT;
        max_levels = levels;
        assert(max_levels);
        membership.resize(fg->get_num_vertices());
        unsigned nthreads =
            atoi(fg->get_configs()->get_option("threads").c_str());
        stage_barrier = barrier::create(nthreads);
        printf("Louvain runs with %u threads ...\n", nthreads);

		graph->start_all(vertex_initializer::ptr(),
            vertex_program_creater::ptr(new louvain_vertex_program_creater()));
		graph->wait4complete();

		// Aggregate the global edge-weight
		std::vector<vertex_program::ptr> progs;
		graph->get_vertex_programs(progs);
		for (vertex_program::ptr vprog : progs) {
			louvain_vertex_program::ptr lvp = louvain_vertex_program::cast2(vprog);
			m += lvp->get_local_m();
		}

        printf("\n\n\x1B[31mThe graph's total edge weight is: %lu ..."
                "\x1B[0m\n\n", m);

        stage = RUN;
		do {
			set_changed(false);

			printf("\n\n\x1B[31m****************** LEVEL ITERATION: %u"
                " ********************************\x1B[0m\n\n", ++cur_level);

            graph->start_all(vertex_initializer::ptr(), vertex_program_creater::
                    ptr(new louvain_vertex_program_creater()));
			graph->wait4complete();
		} while (g_changed && cur_level < max_levels);

		gettimeofday(&end, NULL);

#ifdef PROFILER
		if (!graph_conf.get_prof_file().empty())
			ProfilerStop();
#endif
		printf("Louvain modularity takes %.5f seconds\n", time_diff(start, end));
		return membership;
	}
}
