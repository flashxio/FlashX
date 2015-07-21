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

using namespace fg;

// NOTE: This routine is only meant for undirected graphs!!!
namespace {
// INIT lets us accumulate
enum stage_t
{
	INIT,
	RUN,
	UPDATE_CLUSTERS,
	REBUILD,
};

class cluster
{
	private:
	friend std::ostream& operator<< (std::ostream& ost, const cluster& clust);

	uint32_t weight;
	uint32_t volume;

	public:
	cluster() {

	}

	cluster(uint32_t volume, uint32_t weight) {
		this->weight = weight;
		this->volume = volume;
	}

	void weight_pe(uint32_t weight) {
		this->weight += weight;
	}

	const uint32_t get_weight() const {
		return weight;
	}

	void volume_pe(uint32_t volume) {
		this->volume += volume;
	}

	const uint32_t get_volume() const {
		return volume;
	}
};

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

// If anything changes cluster we cannot converge
static void set_changed(bool changed) {
	if (!changed) // Only once per iter so no prob with this
		g_changed = changed;
	if (changed && (!g_changed))
		g_changed = changed;
}


// Easy printer for cluster
static void print_cluster_map(cluster_map& cl_map) {

	if (cl_map.size() == 0) {
		BOOST_LOG_TRIVIAL(info) << "EMPTY";
		return;
	}
	for (cluster_map::iterator it = cl_map.begin(); it != cl_map.end(); ++it) {
		BOOST_LOG_TRIVIAL(info) << "Cluster: " << it->first << ", " << it->second;
	}
}

static void print_cluster_map_sum(cluster_map& cl_map) {
	uint32_t weight_sum = 0;
	uint32_t vol_sum = 0;

	for (cluster_map::iterator it = cl_map.begin(); it != cl_map.end(); ++it) {
		weight_sum += it->second.get_weight();
		vol_sum += it->second.get_volume();
	}
	BOOST_LOG_TRIVIAL(info) << "Cluster map weight sum: " << weight_sum;
	BOOST_LOG_TRIVIAL(info) << "Cluster map volume sum: " << vol_sum;
}

class louvain_vertex: public compute_vertex
{
	cluster_id_t cluster_id; // current cluster
	float max_modularity; // max modularity to any cluster
	uint32_t weight;
	uint32_t volume;

	public:
	louvain_vertex(vertex_id_t id): compute_vertex(id) {
		cluster_id = id;
		max_modularity = 0;
		weight = 0;
		volume = 0;
	}

	void run(vertex_program &prog) {
		vertex_id_t id = prog.get_vertex_id(*this);
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg1);

	void compute_modularity(cluster_id_t neigh_cluster_id, vertex_id_t my_id);
	void compute_per_vertex_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog);
	void compute_per_cluster_vol_weight(vertex_program &prog);

	// Used for save_query join
	const cluster_id_t get_result() const {
		return cluster_id;
	}
};

// TODO: Opt -- create 2 of these one for INIT and other for others
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

	const uint32_t get_local_ec() const {
		return th_local_edge_count.get_count();
	}

	void update(cluster_id_t id, float vol, uint32_t weight) {

#if 0
		BOOST_LOG_TRIVIAL(info) << "Cluster before update:";
		print_cluster_map(th_local_cluster_map);
		BOOST_LOG_TRIVIAL(info) << "Adding c" << id << ", vol=" << vol << ", wgt=" << weight;
		printf("\n");
#endif

		cluster_map::iterator it = th_local_cluster_map.find(id);

		if (it == th_local_cluster_map.end()) {
			th_local_cluster_map[id] = cluster(vol, weight);
		} else {
			it->second.volume_pe(vol);
			it->second.weight_pe(weight);
		}

#if 0
		BOOST_LOG_TRIVIAL(info) << "Cluster after update:";
		print_cluster_map(th_local_cluster_map);
		printf("\n");
#endif
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

class cluster_message: public vertex_message
{
	cluster_id_t sender_cluster_id;

	public:
	cluster_message(cluster_id_t cluster_id):
		vertex_message(sizeof(cluster_message), false)	{
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
		case UPDATE_CLUSTERS:
			compute_per_cluster_vol_weight(prog);
			break;
		default:
			assert(0);
	}
}

void louvain_vertex::run_on_message(vertex_program &prog, const vertex_message &msg1) {
	
	const cluster_message& msg = (const cluster_message& ) msg1;
	if (msg.get_sender_cluster_id() == cluster_id) { 
		return; // Ignore since my cluster will not change
	} else {
		// This may update the cluster and modularity values
		compute_modularity(msg.get_sender_cluster_id(), prog.get_vertex_id(*this));
	}
}

// Tradeoff by iterating through neighbors, instead of the min computation which is every cluster.
void louvain_vertex::compute_modularity(cluster_id_t neigh_cluster_id, vertex_id_t my_id) {

	// TODO: Optimize these lookups
	cluster curr_cluster = g_cluster_map[cluster_id];
	cluster neigh_cluster = g_cluster_map[neigh_cluster_id];

#if 0
	/** DEBUG **/
	BOOST_LOG_TRIVIAL(info) << "v" << my_id << " in c" << cluster_id << 
		" ==> " << curr_cluster << "processing neigh in cluster: " << neigh_cluster; 
	/** GUBED **/
#endif
	
	// TODO: Remove divisions & hope for no overflow
	float delta_mod = ((neigh_cluster.get_weight() - (curr_cluster.get_weight()-this->weight))
			/ (float) g_edge_weight) +
		(((int)((curr_cluster.get_volume() - this->volume) - neigh_cluster.get_volume()) 
		  * (int)this->volume) / (float)(2*(g_edge_weight*g_edge_weight)));

	if (delta_mod > this->max_modularity) {
#if 0
		/** DEBUG **/
		BOOST_LOG_TRIVIAL(info) << "v" << my_id << " moved from c" <<  cluster_id << " to c"
		   << neigh_cluster_id << ", delta_mod = "<< delta_mod << "\n";
		/** GUBED **/
#endif

		max_modularity = delta_mod;
		this->cluster_id = neigh_cluster_id;
		set_changed(true); // Global notification of cluster change
	}
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
	//BOOST_LOG_TRIVIAL(info) << "v"<< prog.get_vertex_id(*this) << ", " << "c" << cluster_id 
		//<< " += vol: " << volume << ", += wgt: " << weight; 
	((louvain_vertex_program&)prog).update(this->cluster_id, this->volume, this->weight);
}

// General merge function for cluster maps
cluster merge_cluster(cluster& c1, cluster& c2) {
	return cluster(c1.get_volume()+c2.get_volume(), c1.get_weight()+c2.get_weight());
}

template<typename T, typename U >
void print_hash(std::map<T,U>& map) {

	std::cout << "{ ";
	for (typename std::map<T,U>::iterator it = map.begin(); it != map.end(); ++it)
		std::cout << it->first << ":" << it->second << ", ";

	std::cout << "}\n";
}

template <typename T, typename U>
void prune_merge_map (std::map<T, U>& add_map, std::map<T, U>& agg_map,
		U (*merge_func) (U&, U&)) {

	if (agg_map.size() == 0 || add_map.size() == 0) {
		return;
	}

	typename std::map<T,U>::iterator add_it = add_map.begin(); // Always iterate over the add maps keys
	typename std::map<T,U>::iterator agg_it = agg_map.begin();

	for (; add_it != add_map.end(); ++add_it) {
		// skip the keys we don't care about
		while (agg_it->first < add_it->first) {
			if (++agg_it == agg_map.end()) {  // There are no more keys in the agg_map
				break;
			}
		}

		if (add_it->first == agg_it->first) {
			agg_map[agg_it->first] = merge_func(add_it->second, agg_it->second);
		}
	}
}

// FIXME: Building this is serial and likely already slow to
// FIXME: Remove accum_edges eventually 
void build_global_cluster_map (graph_engine::ptr graph, bool accum_edges) {

	std::vector<vertex_program::ptr> ec_progs;
	graph->get_vertex_programs(ec_progs);
	BOOST_FOREACH(vertex_program::ptr vprog, ec_progs) {
		louvain_vertex_program::ptr lvp = louvain_vertex_program::cast2(vprog);

		if (accum_edges) { g_edge_weight += lvp->get_local_ec(); }

		// Merge the volume maps
		cluster (*merge_func) (cluster&, cluster&); // Function pointer to merge a cluster map
		merge_func = &merge_cluster;

#if 0
		BOOST_LOG_TRIVIAL(info) << "Before th_local_cluster_map:";
		print_cluster_map(lvp->get_cluster_map());
		print_cluster_map_sum(lvp->get_cluster_map());
		printf("\n");

		BOOST_LOG_TRIVIAL(info) << "Before g_cluster_map";
		print_cluster_map(g_cluster_map);
		print_cluster_map_sum(g_cluster_map);
		printf("\n");
#endif
		prune_merge_map(lvp->get_cluster_map(), g_cluster_map, merge_func);

		// Merging this clustermap 
		g_cluster_map.insert(lvp->get_cluster_map().begin(), lvp->get_cluster_map().end());

#if 0
		BOOST_LOG_TRIVIAL(info) << "After g_cluster_map";
		print_cluster_map(g_cluster_map);
		print_cluster_map_sum(g_cluster_map);
		printf("\n");
#endif
	}
}

// TODO: Opt -- eliminate memory copies
// TODO: Opt -- omp opts - schedule, loop unroll etc.. 
#if 0
void par_build_global_cluster_map (graph_engine::ptr graph, bool accum_edges) {
	std::vector<vertex_program::ptr> ec_progs;
	graph->get_vertex_programs(ec_progs);


#if 1
	if (accum_edges) {
		BOOST_FOREACH(vertex_program::ptr vprog, ec_progs) {
			louvain_vertex_program::ptr lvp = louvain_vertex_program::cast2(vprog);
			g_edge_weight += lvp->get_local_ec();
		}
	}
#endif

	std::vector<cluster_map> merged_maps;
	merged_maps.resize(ec_progs.size() / 2);

	cluster (*merge_func) (cluster&, cluster&); // Function pointer to merge a cluster map
	merge_func = &merge_cluster;

	//#pragma omp parallel for
	for (uint32_t i = 0; i < ec_progs.size()-1; i+=2) {

		merged_maps[i/2] = prune_merge_map(louvain_vertex_program::cast2(ec_progs[i])->get_cluster_map(), 
								louvain_vertex_program::cast2(ec_progs[i+1])->get_cluster_map(), merge_func);
		if (i == 0 && (ec_progs.size() % 2 != 0)) { // Handle odd # of threads in the first merge iteration
			merged_maps[i/2] = prune_merge_map(louvain_vertex_program::cast2(ec_progs.back())->get_cluster_map(),
														merged_maps[0], merge_func);
		}
	}

	while (merged_maps.size() > 1) {
		//#pragma omp parallel for
		for(uint32_t i = 0; i < merged_maps.size()-1; i+=2) {
			merged_maps[i/2] = prune_merge_map(merged_maps[i], 
					merged_maps[i+1], merge_func);
		}

		if ((merged_maps.size() / 2) % 2 != 0) {
			merged_maps[0] = prune_merge_map(merged_maps[0], 
					merged_maps[(merged_maps.size()/2)-1], merge_func);
		}
		merged_maps.erase(merged_maps.begin() + (merged_maps.size()/2), 
				merged_maps.end()); // Get rid of the back end of the vector
	}
}
#endif

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
		
		build_global_cluster_map(graph, true);

		BOOST_LOG_TRIVIAL(info) << "The graph's total edge weight is " << g_edge_weight << "\n";
#if 0
		/** DEBUG **/
		BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
		BOOST_LOG_TRIVIAL(info) << "Global cluster map: ";
		print_cluster_map(g_cluster_map);
		BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
		/** GUBED **/
#endif
		int iter = 0;

		do {
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute modularity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			louvain_stage = RUN;
			BOOST_LOG_TRIVIAL(info) << "\n\n\x1B[31m****************** LEVEL ITERATION: " << iter++ 
												<< " ********************************\x1B[0m\n\n";

#if 1
			/** DEBUG **/
			FG_vector<cluster_id_t>::ptr ret = FG_vector<cluster_id_t>::create(
					graph->get_num_vertices());
			graph->query_on_all(vertex_query::ptr(new save_query<cluster_id_t, louvain_vertex>(ret)));
			BOOST_LOG_TRIVIAL(info) << "Printing vertex clusters:";
			ret->print(); 
#endif
#if 1
			BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";
			BOOST_LOG_TRIVIAL(info) << "Global cluster map: ";
			print_cluster_map(g_cluster_map);
			print_cluster_map_sum(g_cluster_map);
			BOOST_LOG_TRIVIAL(info) << "\x1B[31m===========================================\x1B[0m\n";

			//if (iter > 2) { fprintf(stderr, "Premature kill"); exit(-1); } 
			/** GUBED **/
#endif

			set_changed(false);
			graph->start_all(vertex_initializer::ptr(), 
				vertex_program_creater::ptr(new louvain_vertex_program_creater()));
			graph->wait4complete();
			g_cluster_map.clear(); // TODO: Opt
		/*~~~~~~~~~~~~~~~~~~~ Update per thread then global cluster maps ~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			louvain_stage = UPDATE_CLUSTERS;
			graph->start_all(vertex_initializer::ptr(), 
				vertex_program_creater::ptr(new louvain_vertex_program_creater()));
			graph->wait4complete();

			// Now aggregate the cluster maps
			build_global_cluster_map(graph, false);
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
