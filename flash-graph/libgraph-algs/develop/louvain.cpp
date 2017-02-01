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
#include <vector>
#include <map>
#include <algorithm>
#include <atomic>

#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"
#include "FG_vector.h"
#include "save_result.h"

using namespace fg;

// NOTE: This routine is only meant for undirected graphs!!!
namespace {
	// INIT lets us accumulate
	enum stage_t
	{
		INIT,
		RUN,
		REBUILD,
	};

	template <typename T>
		class atomicwrapper
		{
			private:
				std::atomic<T> _a;

			public:
				atomicwrapper() :_a() {}

				atomicwrapper(const std::atomic<T> &a) :_a(a.load()) {}

				atomicwrapper(const atomicwrapper &other) :_a(other._a.load()) {}

				atomicwrapper &operator=(const atomicwrapper &other) {
					_a.store(other._a.load());
					return *this;
				}

				T get() {
					return _a;
				}

				void assign(T val) {
					_a = val;
				}

				void minus_eq(T val) {
					_a -= val;
				}

				void plus_eq(T val) {
					_a += val;
				}
		};

	template <typename T>
		static void print_atomicwrapper_v(std::vector<atomicwrapper<T>>& v ) {

			std::cout << "[ ";
			for (typename std::vector<atomicwrapper<T>>::iterator it = v.begin(); it != v.end(); ++it)  {
				std::cout << it->get() << " ";
			}
			std::cout << "]\n";
		}

	// Global map from cluster_id : cluster(volume, weight)
	std::vector<atomicwrapper<uint32_t>> g_weight_vec;
	std::vector<atomicwrapper<uint32_t>> g_volume_vec;

	typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
	typedef vertex_id_t cluster_id_t;


	uint64_t g_edge_weight = 0;
	bool g_changed = false;
	stage_t louvain_stage = INIT; // Init the stage

	// If anything changes cluster we cannot converge
	static void set_changed(bool changed) {
		if (!changed) // Only once per iter so no prob with this
			g_changed = changed;
		if (changed && (!g_changed))
			g_changed = changed;
	}

	// Easy printer for cluster
	template<typename T>
		static void print_vector(typename std::vector<T>& v) {
			std::cout << "[ ";
			for (typename std::vector<T>::iterator it = v.begin(); it != v.end(); ++it) {
				std::cout << *it << " ";
			}
			std::cout << "]\n";
		}

	class louvain_vertex: public compute_vertex
	{
		cluster_id_t cluster_id; // current cluster
		cluster_id_t next_cluster_id; // current cluster
		float max_modularity; // max modularity to any cluster
		uint32_t weight;
		uint32_t volume;
		bool update_globals;

		public:
		louvain_vertex(vertex_id_t id): compute_vertex(id) {
			cluster_id = id;
			next_cluster_id = id;
			max_modularity = 0;
			weight = 0;
			volume = 0;
			update_globals = false;
		}

		void run(vertex_program &prog) {
			if (update_globals) {
				if (cluster_id != next_cluster_id)
					switch_clusters();
				else
					BOOST_LOG_TRIVIAL(info) << "ignoring ...";
				update_globals = false;
			} else {
				vertex_id_t id = prog.get_vertex_id(*this);
				request_vertices(&id, 1);
			}
		}

		void run(vertex_program &prog, const page_vertex &vertex);
		void run_on_message(vertex_program &prog, const vertex_message &msg1);

		void compute_modularity(cluster_id_t neigh_cluster_id, vertex_id_t my_id, vertex_program& prog);
		void compute_per_vertex_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
				vertex_program &prog);
		void compute_per_cluster_vol_weight(vertex_program &prog);

		// Used for save_query join
		const cluster_id_t get_result() const {
			return cluster_id;
		}

		// Remove vertex from old cluster 
		void switch_clusters () {
			BOOST_LOG_TRIVIAL(info) << " Switching! Moving from c" << cluster_id << " to c" << next_cluster_id;
			g_weight_vec[cluster_id].minus_eq(weight);
			g_volume_vec[cluster_id].minus_eq(volume);

			cluster_id = next_cluster_id;
			g_weight_vec[cluster_id].plus_eq(weight);
			g_volume_vec[cluster_id].plus_eq(volume);
		}

		//void notify_iteration_end(vertex_program &vprog);

	};

	/* We need this to get the total edge_weight of the graph */
	class init_vertex_program: public vertex_program_impl<louvain_vertex>
	{
		// Thread local
		edge_count th_local_edge_count; // For the global edge count

		public:
		init_vertex_program() {
			th_local_edge_count = 0;
		}

		typedef std::shared_ptr<init_vertex_program> ptr;

		static ptr cast2(vertex_program::ptr prog) {
			return std::static_pointer_cast<init_vertex_program, vertex_program>(prog);
		}

		void pp_ec(edge_count weight) {
			this->th_local_edge_count += weight;
		}

		const uint32_t get_local_ec() const {
			return th_local_edge_count.get_count();
		}
	};

	/* We need this to do the INIT phase and pass vertex_program_impl to the start method */
	class init_vertex_program_creater: public vertex_program_creater
	{
		public:
			vertex_program::ptr create() const {
				return vertex_program::ptr(new init_vertex_program());
			}
	};

	class cluster_message: public vertex_message
	{
		cluster_id_t sender_cluster_id;

		public:
		cluster_message(cluster_id_t cluster_id):
			vertex_message(sizeof(cluster_message), true)	{
				this->sender_cluster_id = cluster_id;
			}

		const cluster_id_t get_sender_cluster_id() const {
			return sender_cluster_id;
		}
	};

	void louvain_vertex::run(vertex_program &prog, const page_vertex &vertex) {
		switch (louvain_stage) {
			case INIT: /* INIT just accums the global edge_count */
				{
					// Out edges
					data_seq_iterator weight_it = 
						((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);
					edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
					compute_per_vertex_vol_weight(weight_it, id_it, prog);

					// In edges
					weight_it = ((const page_directed_vertex&)vertex).
						get_data_seq_it<edge_count>(IN_EDGE);
					id_it = vertex.get_neigh_seq_it(IN_EDGE);
					compute_per_vertex_vol_weight(weight_it, id_it, prog);
				}
				break;
			case RUN:
				{
					// Broadcast cluster_id to all my neighbors then rely on run_on_message
					// Out edges
					edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
					cluster_message out_msg(this->cluster_id);
					prog.multicast_msg(id_it, out_msg);

					// In edges
					id_it = vertex.get_neigh_seq_it(IN_EDGE);
					cluster_message in_msg(this->cluster_id);
					prog.multicast_msg(id_it, in_msg);
				}
				break;
			default:
				assert(0);
		}
	}

	void louvain_vertex::run_on_message(vertex_program &prog, const vertex_message &msg1) {

		if (!update_globals) 
			update_globals = true;

		const cluster_message& msg = (const cluster_message& ) msg1;
		if (msg.get_sender_cluster_id() == next_cluster_id) { 
			return; // Ignore since my cluster will not change
		} else {
			// This may update the next_cluster_id and max_modularity values
			compute_modularity(msg.get_sender_cluster_id(), prog.get_vertex_id(*this), prog);

		}
	}

	// Tradeoff by iterating through neighbors, instead of the min computation which is every cluster.
	void louvain_vertex::compute_modularity(cluster_id_t neigh_cluster_id, vertex_id_t my_id, vertex_program& prog) {

		// TODO: Remove divisions & hope for no overflow
		float delta_mod = ((int)(g_weight_vec[neigh_cluster_id].get() - 
					(g_weight_vec[this->next_cluster_id].get() - this->weight)) / (float) g_edge_weight) +
			(((int)((g_volume_vec[this->next_cluster_id].get() - this->volume) - g_volume_vec[neigh_cluster_id].get()) 
			  * (int)this->volume) / (float)(2*(g_edge_weight*g_edge_weight)));

		if (delta_mod > this->max_modularity) {
#if 1
			/** DEBUG **/
			BOOST_LOG_TRIVIAL(info) << "v" << my_id << " tentative move from c" <<  cluster_id << " to c"
				<< neigh_cluster_id << ", delta_mod="<< delta_mod << " > max_mod=" << max_modularity << "\n";
			/** GUBED **/
#endif
			max_modularity = delta_mod;
			this->next_cluster_id = neigh_cluster_id;
			//BOOST_LOG_TRIVIAL(info) << "Activating v" << my_id << " for the switch iteration";
			set_changed(true); // Global notification of cluster change

			//prog.request_notify_iter_end(*this);
			//prog.activate_vertices(&my_id,1); // Activate myself so I can change my clusters metadata
		}
	}

#if 0
	void louvain_vertex::notify_iteration_end(vertex_program &prog) {
		assert (next_cluster_id != cluster_id);
		vertex_id_t my_id = prog.get_vertex_id(*this);
		prog.activate_vertices(&my_id,1); // Activate myself so I can change my clusters metadata in NEXT iteration
	}

#endif

	// Only need to do this once per vertex ever
	void louvain_vertex::compute_per_vertex_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
			vertex_program &prog) {
		edge_count local_edge_weight = 0;
		uint32_t self_edge_weight = 0;

		while (weight_it.has_next()) {
			edge_count e = weight_it.next();
			vertex_id_t nid = id_it.next();

			local_edge_weight += e.get_count();
			if (nid == prog.get_vertex_id(*this)) {
				self_edge_weight += e.get_count();
			}
		}

		((init_vertex_program&)prog).
			pp_ec(local_edge_weight);
		this->volume += local_edge_weight.get_count() + (2*self_edge_weight);
		this->weight += local_edge_weight.get_count() + self_edge_weight;

		// Aggregate for the cluster INIT
		g_weight_vec[cluster_id].plus_eq(local_edge_weight.get_count() + self_edge_weight);
		g_volume_vec[cluster_id].plus_eq(local_edge_weight.get_count() + (2*self_edge_weight));
	}
}

namespace fg 
{
	void compute_louvain(FG_graph::ptr fg, const uint32_t levels)
	{
		graph_index::ptr index = NUMA_graph_index<louvain_vertex>::create(
				fg->get_graph_header());
		graph_engine::ptr graph = fg->create_engine(index);

		BOOST_LOG_TRIVIAL(info) << "Starting Louvain with " << levels << " levels";
		BOOST_LOG_TRIVIAL(info) << "prof_file: " << graph_conf.get_prof_file().c_str();
#ifdef PROFILER
		if (!graph_conf.get_prof_file().empty())
			ProfilerStart(graph_conf.get_prof_file().c_str());
#endif

		struct timeval start, end;
		gettimeofday(&start, NULL);
		// Resize weight and volume vectors
		BOOST_LOG_TRIVIAL(info) << "Resizing vectors to " << graph->get_num_vertices() << "\n";
		g_weight_vec.resize(graph->get_num_vertices(), atomicwrapper<uint32_t>(0));
		g_volume_vec.resize(graph->get_num_vertices(), atomicwrapper<uint32_t>(0));

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute Vol & Weight ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		graph->start_all(vertex_initializer::ptr(), 
				vertex_program_creater::ptr(new init_vertex_program_creater()));
		graph->wait4complete();

		// Aggregate the global edge-weight
		std::vector<vertex_program::ptr> ec_progs;
		graph->get_vertex_programs(ec_progs);
		BOOST_FOREACH(vertex_program::ptr vprog, ec_progs) {
			init_vertex_program::ptr lvp = init_vertex_program::cast2(vprog);
			g_edge_weight += lvp->get_local_ec(); 
		}
		BOOST_LOG_TRIVIAL(info) << "The graph's total edge weight is " << g_edge_weight << "\n";

		int iter = 0;
		do {
			set_changed(false);
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute modularity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			louvain_stage = RUN;
			BOOST_LOG_TRIVIAL(info) << "\n\n\x1B[31m****************** LEVEL ITERATION: " << iter++ 
				<< " ********************************\x1B[0m\n\n";

			graph->start_all();
			graph->wait4complete();

#if 1
			/** DEBUG **/
			FG_vector<cluster_id_t>::ptr ret = FG_vector<cluster_id_t>::create(
					graph->get_num_vertices());
			graph->query_on_all(vertex_query::ptr(new save_query<cluster_id_t, louvain_vertex>(ret)));
			BOOST_LOG_TRIVIAL(info) << "Printing vertex clusters:";
			ret->print(); 
			/** GUBED **/
#endif
#if 1
			/** DEBUG **/
			BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
			BOOST_LOG_TRIVIAL(info) << "Global volume vector:";
			print_atomicwrapper_v(g_volume_vec);

			BOOST_LOG_TRIVIAL(info) << "Global weight vector:";
			print_atomicwrapper_v(g_weight_vec);
			BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
			/** GUBED **/
			if (g_changed) BOOST_LOG_TRIVIAL(info) << "Loop continues!";
			else BOOST_LOG_TRIVIAL(info) << "Loop ending!";

			//if (iter > 2) { fprintf(stderr, "Premature kill"); exit(-1); } 
			/** GUBED **/
#endif
		} while (g_changed);

#if 1
		/** DEBUG **/
		FG_vector<cluster_id_t>::ptr ret = FG_vector<cluster_id_t>::create(
				graph->get_num_vertices());
		graph->query_on_all(vertex_query::ptr(new save_query<cluster_id_t, louvain_vertex>(ret)));
		BOOST_LOG_TRIVIAL(info) << "\nVertex clusters @ end of Level1:";
		ret->print(); 
		/** GUBED **/
#endif

		louvain_stage = REBUILD; // TODO: Do something here
		BOOST_LOG_TRIVIAL(info) << "\n Reached rebuild graph stage\n"; 
		gettimeofday(&end, NULL);

#ifdef PROFILER
		if (!graph_conf.get_prof_file().empty())
			ProfilerStop();
#endif

		BOOST_LOG_TRIVIAL(info) << boost::format("It takes %1% seconds")
			% time_diff(start, end);

		return;
	}
}
