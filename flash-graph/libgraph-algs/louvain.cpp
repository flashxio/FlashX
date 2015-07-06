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

using namespace fg;

// NOTE: This script is only meant for undirected graphs!!!
namespace {
typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;
uint64_t g_edge_weight = 0;
bool g_changed = false;

// INIT lets us accumulate
enum stage_t
{
	INIT,
	LEVEL1,
	RUN,
};

stage_t louvain_stage = INIT; // Init the stage
edge_count tot_edge_weight = 0; // Edge weight of the entire graph

class cluster
{
	uint32_t weight;
	uint32_t volume;
	std::vector<vertex_id_t> members;
	pthread_spinlock_t lock;

	public:
	cluster() {
		weight = 0;
		volume = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	void weight_pe(uint32_t weight) {
		this->weight += weight;
	}

	uint32_t get_weight() {
		return weight;
	}

	void volume_pe(uint32_t volume) {
		this->volume += volume;
	}

	uint32_t get_volume() {
		return volume;
	}

	void add_member(vertex_id_t member, uint32_t volume) {
		pthread_spin_lock(&lock); // FIXME: Locking :(
		members.push_back(member);
		this->volume += volume; 
		pthread_spin_unlock(&lock);
	}	

	void remove_member(vertex_id_t member_id, uint32_t volume) {
		pthread_spin_lock(&lock); // FIXME: Locking :(
		std::vector<vertex_id_t>::iterator it = std::find(members.begin(), members.end(), member_id);
		assert(it == members.end()); // TODO: rm -- Should be impossible
		members.erase(it);
		this->volume -= volume;
		pthread_spin_unlock(&lock);
	}

	std::vector<vertex_id_t>& get_members() {
		return this->members;
	}
};

std::map<vertex_id_t, float> g_volume_map; // Per-vertex volume 
std::map<vertex_id_t, cluster> cluster_map; // Keep track of what cluster each vertex is in

class louvain_vertex: public compute_vertex
{
	vertex_id_t cluster; // current cluster
	float modularity;
	uint32_t self_edge_weight;
	uint32_t volume;

	public:
	louvain_vertex(vertex_id_t id): compute_vertex(id) {
		cluster = id;
		modularity = std::numeric_limits<float>::min();
		self_edge_weight = 0;
		volume = 0;
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &, const vertex_message &msg1) {}

	void compute_modularity(edge_seq_iterator& id_it, data_seq_iterator& weight_it, 
		vertex_id_t& max_cluster, float& max_mod, vertex_id_t my_id);
	void compute_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog);

};

/* We need this to get the total edge_weight of the graph */
class louvain_vertex_program: public vertex_program_impl<louvain_vertex>
{
	// Thread local
	edge_count th_local_edge_count;
	std::map<vertex_id_t, float> th_local_volume_map;

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

	void update_volume(vertex_id_t vid, float vol) {
		std::map<vertex_id_t, float>::iterator it = th_local_volume_map.find(vid);

		if (it == th_local_volume_map.end()) {
			th_local_volume_map[vid] = vol;
		} else {
			it->second += vol;
		}
	}

	std::map<vertex_id_t, float>& get_vol_map() {
		return th_local_volume_map;
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
	// Need to recompute the modularity
	vertex_id_t id = prog.get_vertex_id(*this);
	request_vertices(&id, 1);
}

// FIXME: This will only work for the first level
void louvain_vertex::compute_modularity(edge_seq_iterator& id_it, data_seq_iterator& weight_it, 
		vertex_id_t& max_cluster, float& max_mod, vertex_id_t my_id) {
	float delta_mod;

	// Iterate through all vertices in the g_volume_map & if one is my 
	//	neighbor then I need to modify it's volume to not include me.
	std::map<vertex_id_t, float>::iterator graph_it = g_volume_map.begin(); // Iterator for all graph vertices

	while(id_it.has_next()) {
		vertex_id_t	nid = id_it.next();
		edge_count e = weight_it.next();

		delta_mod = ((e.get_count() - 0) / g_edge_weight) + 
			(((g_volume_map[my_id] - this->self_edge_weight) - 
			  (g_volume_map[nid] - e.get_count())) * g_volume_map[my_id])
			/ (float)(2*(g_edge_weight^2));

		BOOST_LOG_TRIVIAL(info) << "v" << my_id << " delta_mod for v" << graph_it->first << " = " << delta_mod;

		if (delta_mod > max_mod) {
			max_mod = delta_mod;
			max_cluster = graph_it->first;
		}
	}
}

// If anything changes cluster we cannot converge
void set_changed(bool changed) {
	if (!g_changed)
		g_changed = changed;
}

void louvain_vertex::compute_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog) {
	edge_count local_edge_weight = 0;

	while (weight_it.has_next()) {
		edge_count e = weight_it.next();
		vertex_id_t nid = id_it.next();

		local_edge_weight += e.get_count();
		if (nid == prog.get_vertex_id(*this)) {
			this->self_edge_weight += e.get_count();
		}
	}

	((louvain_vertex_program&)prog).
		pp_ec(local_edge_weight);
	this->volume = local_edge_weight.get_count() + (2*this->self_edge_weight);

	((louvain_vertex_program&)prog).update_volume(this->cluster,
		local_edge_weight.get_count() + (2*this->self_edge_weight));
}

// NOTE: If I know your ID I know what cluster you're in
void louvain_vertex::run(vertex_program &prog, const page_vertex &vertex) {

	switch (louvain_stage) {
		case INIT: /* INIT just accums the global edge_count. I chose OUT_EDGE at random */
			{
				data_seq_iterator weight_it = 
					((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);
				edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
				compute_vol_weight(weight_it, id_it, prog);
			}
			break;
		case LEVEL1:
			{
				// Compute the new cluster based on modularity
				float max_mod = this->modularity;
				vertex_id_t max_cluster = this->cluster;

				// Ignore's all not connected to this vertex
				edge_seq_iterator id_it = vertex.get_neigh_seq_it(OUT_EDGE);
				data_seq_iterator weight_it = 
					((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(OUT_EDGE);

				compute_modularity(id_it, weight_it, max_cluster, max_mod, prog.get_vertex_id(*this));

				if (this->cluster != max_cluster) {
					BOOST_LOG_TRIVIAL(info) << "Vertex " << prog.get_vertex_id(*this) << " with mod = " 
						<< this->modularity << " < " << max_mod <<
						" ,moved from cluster " << this->cluster << " ==> " << max_cluster << "\n";
					set_changed(true);
				} else {
					BOOST_LOG_TRIVIAL(info) << "Vertex " << prog.get_vertex_id(*this) << " with mod = " << this->modularity
						<< " ,stayed in cluster " << this->cluster << "\n";
				}

				this->cluster = max_cluster;
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

// General Addition function for merging maps
template <typename T>
T add(T arg1, T arg2) {
	  return arg1 + arg2;
}

template <typename T, typename U>
void build_merge_map (std::map<T, U>& add_map, std::map<T, U>& agg_map,
		std::map<T, U>& new_map, U (*merge_func) (U, U)) {
	if (agg_map.size() == 0 || add_map.size() == 0) { return; }

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
		} else {
			new_map[add_it->first] = add_it->second;
		}
	}
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
			std::map<vertex_id_t, float> merge_map;
			float (*add_func) (float, float); // Function pointer to add map
			add_func = &add;
			build_merge_map(lvp->get_vol_map(), g_volume_map, merge_map, add_func);

			g_volume_map.insert(merge_map.begin(), merge_map.end());
		}
#if 1
		BOOST_LOG_TRIVIAL(info) << "The graph's total edge weight is " << g_edge_weight << "\n";
		for (std::map<vertex_id_t, float>::iterator it = g_volume_map.begin();
				it != g_volume_map.end(); it++) {
			BOOST_LOG_TRIVIAL(info) << "Vertex: " << it->first << ", Volume: " << it->second;
		}
		BOOST_LOG_TRIVIAL(info) << "\x1B[31m====================================================\x1B[0m\n";
#endif
		/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute modularity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#if 1
		louvain_stage = LEVEL1;

		do {
			set_changed(false);
			graph->start_all();
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
