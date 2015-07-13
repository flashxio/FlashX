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

// In iterations > 1:
//		If I'm a vertex that joins a cluster with higher_id I send no messages

using namespace fg;

// NOTE: This routine is only meant for undirected graphs!!!
namespace {
// INIT lets us accumulate
enum stage_t
{
	INIT,
	LEVEL1,
	RUN,
	UPDATE,
};

class cluster
{
	private:
	friend std::ostream& operator<< (std::ostream& ost, const cluster& clust);

	uint32_t weight;
	uint32_t volume;

	public:
	cluster(uint32_t volume=0, uint32_t weight=0) {
		this->weight = weight;
		this->volume = volume;
	}

	void weight_pe(uint32_t weight) {
		this->weight += weight;
	}

	uint32_t get_weight() const {
		return weight;
	}

	void volume_pe(uint32_t volume) {
		this->volume += volume;
	}

	uint32_t get_volume() const {
		return volume;
	}
};

// Easy printer for cluster
std::ostream& operator<< (std::ostream& ost, const cluster& clust) {
	ost << "(vol:" << clust.get_volume() << ", wgt:" 
		<< clust.get_weight() << "), ";
	return ost;
}

typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
typedef vertex_id_t cluster_id_t;
typedef std::map<cluster_id_t, cluster> cluster_map; // Keep track of every cluster & it's metadata

uint64_t g_edge_weight = 0;
bool g_changed = false;
stage_t louvain_stage = INIT; // Init the stage

cluster_map g_cluster_map; // Global map from cluster_id : cluster(volume, weight)

static void print_cluster_map(cluster_map cl_map) {

	for (cluster_map::iterator it = cl_map.begin(); it != cl_map.end(); ++it) {
		BOOST_LOG_TRIVIAL(info) << "Cluster: " << it->first << ", " << it->second;
	}
}

class louvain_vertex: public compute_vertex
{
	cluster_id_t cluster_id; // current cluster
	float modularity;
	uint32_t weight;
	uint32_t volume;

	public:
	louvain_vertex(vertex_id_t id): compute_vertex(id) {
		cluster_id = id;
		modularity = std::numeric_limits<float>::min();
		weight = 0;
		volume = 0;
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &, const vertex_message &msg1) {}

	void compute_modularity(edge_seq_iterator& id_it,
		cluster_id_t& max_cluster, float& max_mod, vertex_id_t my_id);
	void compute_per_vertex_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog);
	void compute_per_cluster_vol_weight(vertex_program &prog);

	// Used for save_query join
	cluster_id_t get_result() const {
		return cluster_id;
	}
};

// Each vertex sends a message to its neighbors` 
class cluster_message: public vertex_message
{
	cluster_id_t sender_cluster_id;

	public:
	cluster_message(cluster_id_t id): 
		vertex_message(sizeof(cluster_message), true) {
			this->sender_cluster_id = id;
		}

	const vertex_id_t get_sender_cluster_id() const {
		return sender_cluster_id;
	}
};

/* We need this to get the total edge_weight of the graph */
class louvain_vertex_program: public vertex_program_impl<louvain_vertex>
{
	// Thread local
	edge_count th_local_edge_count; // For the global edge count
	cluster_map th_local_cluster_map;

	public:
	louvain_vertex_program() {
		th_local_edge_count = 0;
	}

	typedef std::shared_ptr<louvain_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<louvain_vertex_program, vertex_program>(prog);
	}

	void pp_ec(edge_count weight) {
		this->th_local_edge_count += weight;
	}

	uint32_t get_local_ec() {
		return th_local_edge_count.get_count();
	}

	void update(cluster_id_t id, float vol, uint32_t weight) {
		// TODO: Case of INIT we always just add since it's one cluster per vertex
		cluster_map::iterator it = th_local_cluster_map.find(id);

		if (it == th_local_cluster_map.end()) {
			th_local_cluster_map[id] = cluster(vol, weight);
		} else {
			it->second.volume_pe(vol);
			it->second.weight_pe(weight);
		}
	}

	cluster_map& get_cluster_map() {
		return th_local_cluster_map;
	}
};

/* We need this to do the INIT phase and pass vertex_program_impl to the start method */
class louvain_vertex_program_creater: public vertex_program_creater
{
	public:
		vertex_program::ptr create() const {
			return vertex_program::ptr(new louvain_vertex_program());
		}
};

void louvain_vertex::run(vertex_program &prog) {
	
	switch (louvain_stage) {
		case INIT:
		case LEVEL1:
		case RUN:
			{
				vertex_id_t id = prog.get_vertex_id(*this);
				request_vertices(&id, 1);
			}
			break;
		case UPDATE:
			break;
		default:
			BOOST_LOG_TRIVIAL(fatal) << "Unknown case louvain stage!";
	}
}

// Tradeoff by iterating through neighbors, instead of the min computation which is every cluster.
void louvain_vertex::compute_modularity(edge_seq_iterator& id_it, 
		cluster_id_t& max_cluster, float& max_mod, vertex_id_t my_id) {

	float delta_mod = 0;
	cluster curr_cluster = g_cluster_map[my_id];

	BOOST_LOG_TRIVIAL(info) << "Vertex " << my_id << " ==> cluster:" << curr_cluster;

	while(id_it.has_next()) {
		vertex_id_t	nid = id_it.next();
		cluster neigh_cluster = g_cluster_map[nid];

		BOOST_LOG_TRIVIAL(info) << "Neighbor " << nid << " ==> cluster: " << neigh_cluster; 

		if (nid == this->cluster_id) { // If I'm in the same cluster as my neigh
			delta_mod = (((neigh_cluster.get_weight() - this->weight) - 
						(curr_cluster.get_weight()-this->weight)) / (float) g_edge_weight) +
						(((int)((curr_cluster.get_volume() - this->volume) - (neigh_cluster.get_volume()-this->volume)) 
						 * (int)this->volume) / (float)(2*(g_edge_weight*g_edge_weight)));
		} else {
			delta_mod = ((neigh_cluster.get_weight() - (curr_cluster.get_weight()-this->weight))
						/ (float) g_edge_weight) +
						(((int)((curr_cluster.get_volume() - this->volume) - neigh_cluster.get_volume()) 
						 * (int)this->volume) / (float)(2*(g_edge_weight*g_edge_weight)));
		}

		BOOST_LOG_TRIVIAL(info) << "v" << my_id << " delta_mod for v" << nid << " = " << delta_mod << "\n";

		if (delta_mod > max_mod) {
			max_mod = delta_mod;
			max_cluster = nid;
		}
	}
}

// If anything changes cluster we cannot converge
void set_changed(bool changed) {
	if (!changed) // Only once per iter so no prob with this
		g_changed = changed;
	if (changed && (!g_changed))
		g_changed = changed;
	else 
		BOOST_LOG_TRIVIAL(fatal) << "Unknown option for g_changed";
}

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

	((louvain_vertex_program&)prog).
		pp_ec(local_edge_weight);
	this->volume += local_edge_weight.get_count() + (2*self_edge_weight);
	this->weight += local_edge_weight.get_count() + self_edge_weight;

	// Add this for the cluster
	((louvain_vertex_program&)prog).update(this->cluster_id, local_edge_weight.get_count() + 
		(2*self_edge_weight), local_edge_weight.get_count() + self_edge_weight);
}

//  Does not require edgelist
void louvain_vertex::compute_per_cluster_vol_weight(vertex_program &prog) {
	((louvain_vertex_program&)prog).update(this->cluster_id, this->volume, this->weight);
}

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
		case LEVEL1:
			{
				// Compute the new cluster based on modularity
				float max_mod = this->modularity;
				vertex_id_t max_cluster = this->cluster_id;

				// Out edges
				edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
				compute_modularity(id_it, max_cluster, max_mod, prog.get_vertex_id(*this));
				
				// In edges
				id_it = vertex.get_neigh_seq_it(IN_EDGE);
				compute_modularity(id_it, max_cluster, max_mod, prog.get_vertex_id(*this));
				
				if (this->cluster_id != max_cluster) {
					BOOST_LOG_TRIVIAL(info) << "Vertex " << prog.get_vertex_id(*this) << " with mod = " 
						<< this->modularity << " < " << max_mod <<
						", moved from cluster " << this->cluster_id << " ==> " << max_cluster << "\n";
					set_changed(true);
				} else {
					BOOST_LOG_TRIVIAL(info) << "Vertex " << prog.get_vertex_id(*this) << " with mod = " << this->modularity
						<< " ,stayed in cluster " << this->cluster_id << "\n";
				}

				this->cluster_id = max_cluster;
				this->modularity = max_mod;	
			}
			break;
		case RUN:
			BOOST_LOG_TRIVIAL(fatal) << "Run unimplemented!";
			break;

		default:
			assert(0);
	}
}

// General merge function for cluster maps
cluster merge_cluster(cluster& c1, cluster& c2) {
	return cluster(c1.get_volume()+c2.get_volume(), c2.get_weight()+c2.get_weight());
}

template<typename T, typename U >
void print_hash(std::map<T,U>& map) {

	std::cout << "{ ";
	for (typename std::map<T,U>::iterator it=map.begin(); it != map.end(); ++it)
		std::cout << it->first << ":" << it->second << ", ";

	std::cout << "}\n";
}

template <typename T, typename U>
std::map<T, U> build_merge_map (std::map<T, U>& add_map, std::map<T, U>& agg_map,
		 U (*merge_func) (U&, U&)) {

	if (agg_map.size() == 0 || add_map.size() == 0) {
		return add_map;
	}

	std::map<T,U> new_map;

	typename std::map<T,U>::iterator add_it = add_map.begin(); // Always iterate over the add maps keys
	typename std::map<T,U>::iterator agg_it = agg_map.begin();

	for (; add_it != add_map.end(); ++add_it) {

		// skip the keys we don't care about
		while (agg_it->first < add_it->first) {
			if (++agg_it == agg_map.end()) {  // There are no more keys in the agg_map
				new_map.insert(add_it, add_map.end()); // Get the rest of the map
				break;
			}
		}

		if (add_it->first == agg_it->first) {
			new_map[add_it->first] = merge_func(add_it->second, agg_it->second);

			BOOST_LOG_TRIVIAL(info) << "new_map add:"; print_cluster_map(new_map);
		} else {
			BOOST_LOG_TRIVIAL(info) << "new_map new:"; print_cluster_map(new_map);
			new_map[add_it->first] = add_it->second;
		}
	}

	return new_map;
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

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute Vol & Weight ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		graph->start_all(vertex_initializer::ptr(), 
				vertex_program_creater::ptr(new louvain_vertex_program_creater()));
		graph->wait4complete();
		std::vector<vertex_program::ptr> ec_progs;
		graph->get_vertex_programs(ec_progs);
		BOOST_FOREACH(vertex_program::ptr vprog, ec_progs) {
			louvain_vertex_program::ptr lvp = louvain_vertex_program::cast2(vprog);
			g_edge_weight += lvp->get_local_ec();
			
			// Merge the volume maps
			
			cluster (*merge_func) (cluster&, cluster&); // Function pointer to merge a cluster map
			merge_func = &merge_cluster;

			cluster_map merge_map = build_merge_map(lvp->get_cluster_map(), g_cluster_map, merge_func);
			g_cluster_map.insert(merge_map.begin(), merge_map.end());
		}
#if 1
		BOOST_LOG_TRIVIAL(info) << "The graph's total edge weight is " << g_edge_weight << "\n";

		BOOST_LOG_TRIVIAL(info) << "Global cluster map: ";
		print_cluster_map(g_cluster_map);

		BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
#endif
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute modularity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#if 1
		louvain_stage = LEVEL1;

		int iter = 0;

		do {
			/** DEBUG **/
			BOOST_LOG_TRIVIAL(info) << "\n\n******************** ITERATION: " << iter++ 
												<< " ********************************\n\n";
#if 1
			FG_vector<cluster_id_t>::ptr ret = FG_vector<cluster_id_t>::create(
					graph->get_num_vertices());
			graph->query_on_all(vertex_query::ptr(new save_query<cluster_id_t, louvain_vertex>(ret)));
			BOOST_LOG_TRIVIAL(info) << "Printing vertex clusters:";
			ret->print(); 

			if (iter > 5) { fprintf(stderr, "Premature kill"); exit(-1); } 
#endif
			/** DEBUG **/

			set_changed(false);
			graph->start_all(); // TODO: Init correctly so we get per vertex weight & volume
			graph->wait4complete();
		} while (g_changed);

		louvain_stage = RUN;
		BOOST_LOG_TRIVIAL(info) << "\n Reached running stage\n"; 
#endif

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
