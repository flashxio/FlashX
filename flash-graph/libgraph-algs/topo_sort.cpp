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
#include <set>
#include <algorithm>

#include "graph_engine.h"
#include "FGlib.h"
#include "topological.hpp"
#include "helper.hpp"

using namespace fg;
using namespace monya;

namespace {
enum stage_t
{
	INIT,
	TRANSMIT,
};
stage_t stage;

typedef monya::IndexVector<IndexVal<vertex_id_t, vsize_t>, vsize_t, vsize_t>
    TopoVector;
TopoVector topo_v;
std::vector<vsize_t> part_index;

class topo_vertex: public compute_vertex
{
	vsize_t out_degree;
    //// A neighbor with whom you have a collision
    //std::vector<vertex_id_t> collision;

	public:
	topo_vertex(vertex_id_t id): compute_vertex(id), out_degree(0) {
    }

	//vsize_t get_result() const {
        //return pos;
	//}

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg);

	void run_on_vertex_header(vertex_program &prog,
            const vertex_header &header) {
        out_degree = ((const directed_vertex_header&) header).
                                                get_num_out_edges();
        auto id = header.get_id();
        topo_v[id].set(id, out_degree); // No conflicts here
	}

    void multicast_degree_msg(vertex_program &prog,
            const page_vertex &vertex);
};

class degree_message: public vertex_message
{
    vsize_t degree;
    vertex_id_t sender_id;

	public:
    degree_message(const vsize_t& od, const vertex_id_t& sender_id):
        vertex_message(sizeof(degree_message), false),
        degree(od), sender_id(sender_id) {
    }

    const vertex_id_t& get_sender_id() const {
        return sender_id;
    }

    const vsize_t& get_degree() const {
        return degree;
    }
};

void topo_vertex::multicast_degree_msg(vertex_program &prog,
		const page_vertex &vertex)
{
    // TODO: Choose which vertices receive messages. Maybe vid > mine ?
	int num_dests = vertex.get_num_edges(IN_EDGE);
	edge_seq_iterator it = vertex.get_neigh_seq_it(IN_EDGE, 0, num_dests);

	degree_message msg(out_degree, vertex.get_id());
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

// TODO: Test me for correctness
void swap(TopoVector& tv, const size_t& part_id,
        const vertex_id_t& id1, const vertex_id_t& id2) {

    auto start_idx = tv.get_ranges()[part_id];
    auto end_offset = part_id == tv.get_max_part_id() ?
                        tv.locks_size() : tv.get_ranges()[part_id+1];

    //vertex_id_t find_id1 = id1 < id2 ? id1 : id2;
    //vertex_id_t find_id2 = find_id1 == id1 ? id2 : id1;

    tv.lock(part_id);
    // This finds each and swaps them with a locked partition
    tv.swap(id1, id2, start_idx, end_offset);
    tv.unlock(part_id);
}

void topo_vertex::run_on_message(vertex_program &prog,
        const vertex_message &msg1) {
	if (stage == TRANSMIT) {
        const degree_message &msg = (const degree_message &) msg1;
        const vertex_id_t id = prog.get_vertex_id(*this);

        if (msg.get_degree() == out_degree) {

#if 0
            printf("Running vid: %u (%u), found vid: %u (%u) ...\n",
                    id, out_degree, msg.get_sender_id(), msg.get_degree());
#endif

            if (msg.get_sender_id() > id && msg.get_sender_id() != id) {
                auto part_id = part_index[id]; // What partition are you in

                assert(part_id == part_index[msg.get_sender_id()]);
#if 0
                printf("\t ===> Swapping vid: %u with vid: %u in part: %u!\n",
                        id, msg.get_sender_id(), part_id);
#endif
                swap(topo_v, part_id, id, msg.get_sender_id());
            }
        }
    }
}

}

namespace fg
{

std::vector<vertex_id_t> compute_topo_sort(FG_graph::ptr fg, bool approx)
{
	graph_index::ptr index = NUMA_graph_index<topo_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	struct timeval start, end;
	gettimeofday(&start, NULL);

    topo_v.resize(fg->get_num_vertices());

	printf("Running the init degree stage ...\n");
	stage = INIT;
	graph->start_all();
	graph->wait4complete();
    topo_v.sort();

#if 0
    printf("\n\nApproximate sort:\n");
    topo_v.print();
#endif

    if (!approx) {
        part_index.resize(fg->get_num_vertices());

        printf("\nPartitioning the topological vector ...\n");
        topo_v.partition(part_index);

        printf("\nWe have %lu partitions ...\n", topo_v.get_ranges().size());
#if 1
        printf("Partition ranges:\n");
        print(topo_v.get_ranges());
        printf("\nPartition index:\n");
        print(part_index);
        printf("\n");
        std::set<size_t> parts(part_index.begin(), part_index.end());
        printf("parts.size() = %lu , topo_v.get_ranges.size() = %lu\n",
                parts.size(), topo_v.get_ranges().size());
        assert(parts.size() == topo_v.get_ranges().size());
#endif
        printf("Running transmit degree stage ...\n");
        stage = TRANSMIT;
        topo_v.set_sorted(false); // This stage immediately unsorts
        graph->start_all();
        graph->wait4complete();
    }

	gettimeofday(&end, NULL);
    printf("Topological sort took %.5f sec to complete\n",
            time_diff(start, end));

    topo_v.print();

    std::vector<vertex_id_t> res;
    topo_v.get_indexes(res);
	return res;
}
}
