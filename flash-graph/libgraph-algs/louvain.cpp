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

// We will ignore the race conditions first
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
edge_count tot_edge_weight = 0; // E

class louvain_vertex: public compute_vertex
{
	vertex_id_t cluster; // current cluster
	float modularity;
	uint64_t volume;

	public:
	louvain_vertex(vertex_id_t id): compute_vertex(id) {
		cluster = id;
		modularity = 0;
		volume = 0;
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &, const vertex_message &msg1) {}

	void compute_modularity(edge_seq_iterator& id_it, data_seq_iterator& weight_it, 
		vertex_id_t& max_cluster, float& max_mod);
	void compute_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog);

	// void move();
};

/* We need this to get the total edge_weight of the graph */
class louvain_vertex_program: public vertex_program_impl<louvain_vertex>
{
	edge_count th_local_edge_count; // Local

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
		BOOST_LOG_TRIVIAL(info) << "Current global thread count " << this->th_local_edge_count << "\n";
	}

	uint32_t get_local_ec() {
		return th_local_edge_count.get_count();
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

// NOTE: This will only work for the first level
void louvain_vertex::compute_modularity(edge_seq_iterator& id_it, data_seq_iterator& weight_it, 
		 vertex_id_t& max_cluster, float& max_mod) {
	float delta_mod;

	while (weight_it.has_next()) {
		vertex_id_t nid = id_it.next();
		edge_count e = weight_it.next();

		delta_mod = (e.get_count()/g_edge_weight) + (0.0000 - (0.0000 * this->volume))/((2*g_edge_weight)^2); // FIXME
		if (delta_mod > max_mod) {
			max_mod = delta_mod;
			max_cluster = nid;
		}
	}
}

// If anything changes cluster we cannot converge
void toggle_changed() {
	if (!g_changed)
		g_changed = true;
}

void louvain_vertex::compute_vol_weight(data_seq_iterator& weight_it, edge_seq_iterator& id_it, 
		vertex_program &prog) {
	edge_count local_edge_weight = 0;
	edge_count self_edge_weight = 0;

	while (weight_it.has_next()) {
		edge_count e = weight_it.next();
		vertex_id_t nid = id_it.next();

		local_edge_weight += e.get_count();
		if (nid  == prog.get_vertex_id(*this)) {
			self_edge_weight += e.get_count();
		}
	}

	((louvain_vertex_program&)prog).
		pp_ec(local_edge_weight);
	this->volume = local_edge_weight.get_count() + (2*self_edge_weight.get_count());

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

				compute_modularity(id_it, weight_it, max_cluster, max_mod);

				// TODO: Flag one has changed so we do not converge
				if (this->cluster != max_cluster) {
					toggle_changed();
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

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
		graph->start_all(vertex_initializer::ptr(), 
				vertex_program_creater::ptr(new louvain_vertex_program_creater()));
		graph->wait4complete();
		std::vector<vertex_program::ptr> ec_progs;
		graph->get_vertex_programs(ec_progs);
		BOOST_FOREACH(vertex_program::ptr vprog, ec_progs) {
			louvain_vertex_program::ptr lvp = louvain_vertex_program::cast2(vprog);
			g_edge_weight += lvp->get_local_ec();
		}
		BOOST_LOG_TRIVIAL(info) << "The graph's total edge weight is " << g_edge_weight << "\n";

		/*~~~~~~~~~~~~~~~~~~~~~~~~~~ Compute E ~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

		louvain_stage = LEVEL1;
		graph->start_all();
		graph->wait4complete();

		


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
