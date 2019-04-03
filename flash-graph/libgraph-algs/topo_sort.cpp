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
#include <mutex>

#include "graph_engine.h"
#include "FGlib.h"
#include "topological.hpp"

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
std::vector<size_t> part_ranges;
size_t max_part_id;

// This class contains a bunch of locks for each partition
class TopoVectorLocks {
public:
    std::vector<std::mutex*> locks;
    TopoVectorLocks() { }

    void resize(const size_t nlocks) {
        for (size_t i = 0; i < nlocks; i++) {
            locks.push_back(new std::mutex());
        }
    }

    void lock(const size_t lock_num) {
        locks[lock_num]->lock();
    }

    void unlock(const size_t lock_num) {
        locks[lock_num]->unlock();
    }

    const size_t size() const {
        return locks.size();
    }

    ~TopoVectorLocks() {
        for (auto& lock : locks)
            delete lock;
    }
};

TopoVectorLocks topo_v_locks;

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
void swap(TopoVector& tv, const size_t& part_id, TopoVectorLocks& locks,
        std::vector<size_t> ranges,
        const vertex_id_t& id1, const vertex_id_t& id2) {

    auto start_idx = ranges[part_id];
    auto end_offset = part_id == max_part_id ? locks.size() : ranges[part_id+1];
    vertex_id_t find_id1 = id1 < id2 ? id1 : id2;
    vertex_id_t find_id2 = find_id1 == id1 ? id2 : id1;

    locks.lock(part_id);
    //auto find_it1 = tv.find_index(find_id1, tv.begin()+start_idx,
            //tv.begin()+end_offset);
    //auto find_it2 = tv.find_index(find_id2, find_it1, tv.begin()+end_offset);
    //tv.swap(*find_it1, *find_it2);
    tv.swap(0, 1);
    locks.unlock(part_id);
}

void topo_vertex::run_on_message(vertex_program &prog,
        const vertex_message &msg1) {
	if (stage == TRANSMIT) {
        const degree_message &msg = (const degree_message &) msg1;
        const vertex_id_t id = prog.get_vertex_id(*this);

        if (msg.get_degree() == out_degree) {

            printf("Running vid: %u (%u), found vid: %u (%u) ...\n",
                    id, out_degree, msg.get_sender_id(), msg.get_degree());

            if (msg.get_sender_id() > id && msg.get_sender_id() != id) {
                auto part_id = part_index[id]; // What partition are you in

                assert(part_id == part_index[msg.get_sender_id()]);
                printf("\t ===> Swapping vid: %u with vid: %u in part: %u!\n",
                        id, msg.get_sender_id(), part_id);
                swap(topo_v, part_id, topo_v_locks, part_ranges,
                        id, msg.get_sender_id());
            }
        }
    }
}

// Obtain the partition mapping
void partition_topo_v(TopoVector& v, std::vector<size_t>& ranges,
        std::vector<vertex_id_t>& part_idx) {
    if (v.size())
        ranges.push_back(0);
    else
        return;

    vsize_t part = 0;

    auto prev_val = v[0].get_val();
    part_idx[v[0].get_index()] = part;

    for (size_t i = 1; i < v.size(); i++) {
        if (v[i].get_val() != prev_val) {
            part++;

            ranges.push_back(i);
            prev_val = v[i].get_val();
        }
        part_idx[v[i].get_index()] = part;
    }
    max_part_id = part;
}

// Helpers
template <typename T>
void p(std::vector<T>& v) {
    std::cout << "[ ";
    for (const auto& _ : v)
        std::cout << _ << " ";
    std::cout << "]\n";
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

    if (!approx) {
        part_index.resize(fg->get_num_vertices());

        printf("\n\nPartitioning the topological vector ...\n");
        partition_topo_v(topo_v, part_ranges, part_index);

        printf("\n\nWe have %lu partitions ...\n", part_ranges.size());
        topo_v_locks.resize(part_ranges.size());
#if 0
        printf("Partition ranges:\n");
        p(part_ranges);
        printf("\nPartition index:\n");
        p(part_index);
        printf("\n");
        std::set<size_t> parts(part_index.begin(), part_index.end());
        std::cout << "parts.size() = " << parts.size() <<
            ", part_ranges.size() = " << part_ranges.size() << std::endl;
        assert(parts.size() == part_ranges.size());
#endif
        printf("Running transmit degree stage ...\n");
        stage = TRANSMIT;
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
