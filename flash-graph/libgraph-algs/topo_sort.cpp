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

#include <vector>
#include <algorithm>

#include "graph_engine.h"
//#include "graph_config.h"
#include "FGlib.h"
//#include "save_result.h"

using namespace fg;

namespace {
enum stage_t
{
	INIT,
	TRANSMIT,
};
stage_t stage;

class topo_vertex: public compute_vertex
{
	vsize_t out_degree;
    vsize_t pos; // TODO: Stub

	public:
	topo_vertex(vertex_id_t id): compute_vertex(id), out_degree(0),
        pos(INVALID_VERTEX_ID) {

    }

	vsize_t get_result() const {
        return pos;
	}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg);

	void run_on_vertex_header(vertex_program &prog,
            const vertex_header &header) {
		out_degree = ((const directed_vertex_header&) header).get_num_out_edges();
	}

    void multicast_degree_msg(vertex_program &prog,
            const page_vertex &vertex);
};

class degree_message: public vertex_message
{
    vsize_t degree;

	public:
    degree_message(vsize_t od):
        vertex_message(sizeof(degree_message), true), degree(od) {
    }

    const vsize_t get_degree() const {
        return degree;
    }
};

void topo_vertex::multicast_degree_msg(vertex_program &prog,
		const page_vertex &vertex)
{
    // TODO: We can choose which vertices receive messages. We can do those with
    //  vid > mine
	int num_dests = vertex.get_num_edges(OUT_EDGE);
	edge_seq_iterator it = vertex.get_neigh_seq_it(OUT_EDGE, 0, num_dests);

	// Doesn't matter who sent it, just --degree on reception
	degree_message msg(out_degree);
	prog.multicast_msg(it, msg);
}

void topo_vertex::run(vertex_program &prog) {
	vertex_id_t id = prog.get_vertex_id(*this);
	if (stage == INIT) {
		request_vertex_headers(&id, 1);
		return;
	} else if (stage == TRANSMIT) {
		request_vertices(&id, 1); // put my edgelist in page cache
    } else {
        throw std::runtime_error("Unknown topological sort phase\n");
    }
}

void topo_vertex::run(vertex_program &prog, const page_vertex &vertex) {
	if (stage == TRANSMIT) {
        multicast_degree_msg(prog, vertex);
    }
}

void topo_vertex::run_on_message(vertex_program &prog,
        const vertex_message &msg1) {
	if (stage == TRANSMIT) {
        const degree_message &msg = (const degree_message &) msg1;
        if (msg.get_degree() == out_degree) {
            // FIXME: This is where we tie break using whether the vertex is an
            //  out-neighbor or not
        }
    }
}
}

namespace fg
{

std::vector<vsize_t> compute_topo_sort(FG_graph::ptr fg)
{
	graph_index::ptr index = NUMA_graph_index<topo_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	struct timeval start, end;
	gettimeofday(&start, NULL);

	printf("Running the init degree stage ...\n");
	stage = INIT;
	graph->start_all();
	graph->wait4complete();

    printf("Running transmit degree stage ...\n");
	stage = TRANSMIT;
	graph->start_all();
	graph->wait4complete();


	gettimeofday(&end, NULL);
    printf("Topological sort took %.5f sec to complete.s\n",
            time_diff(start, end));

    std::vector<vsize_t> res(fg->get_num_vertices());
    // FIXME: Get total ordering
	return res;
}
}
