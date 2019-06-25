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

#include "graph_engine.h"
#include "FGlib.h"

using namespace fg;

namespace {

edge_type g_edge_type = edge_type::BOTH_EDGES;

class closeness_vertex: public compute_directed_vertex
{
	vsize_t dist;

	public:
	closeness_vertex(vertex_id_t id): compute_directed_vertex(id) {
        // Note we start with all vertices at max dist
		dist = std::numeric_limits<vsize_t>::max();
	}

	size_t get_result() const {
		return dist;
	}

    const vsize_t get_dist() const {
        return dist;
    }

    void reset() {
		dist = std::numeric_limits<vsize_t>::max();
    }

	void run(vertex_program &prog);
	void run(vertex_program &prog, const page_vertex &vertex);
	void run_on_message(vertex_program &prog, const vertex_message &msg) {};
    void multicast_dist_msg(vertex_program &prog,
            const page_vertex &vertex);
};

// An activateion message for neighbors
class dist_message: public vertex_message
{
    public:
        // auto-activate
        dist_message(): vertex_message(sizeof(dist_message), true) {
        }
};

void closeness_vertex::multicast_dist_msg(vertex_program &prog,
		const page_vertex &vertex)
{
	int num_dests = vertex.get_num_edges(g_edge_type);
	edge_seq_iterator it = vertex.get_neigh_seq_it(g_edge_type, 0, num_dests);

	dist_message msg;
	prog.multicast_msg(it, msg);
}

void closeness_vertex::run(vertex_program &prog, const page_vertex &vertex) {
    //auto vid = prog.get_vertex_id(*this);
    //printf("Vertex: %u active in run(arg,arg) with dist: %u\n", vid, dist);
    multicast_dist_msg(prog, vertex);
}

// We can accumulate at the very end ...
class accuml_query: public vertex_query
{
    size_t dist;
    public:
    accuml_query() {
        dist = 0;
    }

     void run(graph_engine &graph, compute_vertex &v) override {
        closeness_vertex &closeness_v = (closeness_vertex &) v;
            dist += closeness_v.get_dist();
    }

    void merge (graph_engine &graph, vertex_query::ptr q) override {
        accuml_query *a = (accuml_query *) q.get();
        dist += a->get_dist();
    }

    ptr clone() override {
        return vertex_query::ptr(new accuml_query());
    }

    size_t get_dist() const {
        return dist;
    }
};

// We can also accumulate on the fly
class closeness_vertex_program : public vertex_program_impl<closeness_vertex>
{
    size_t cuml_dist;

    public:
    typedef std::shared_ptr<closeness_vertex_program> ptr;

    closeness_vertex_program() : cuml_dist(0) { }

    void cuml_dist_peq(const vsize_t dist) {
        cuml_dist += dist;
    }

    const size_t get_cuml_dist() const {
        return cuml_dist;
    }

    static ptr cast2(vertex_program::ptr prog) {
        return std::static_pointer_cast<closeness_vertex_program,
               vertex_program>(prog);
    }
};

class closeness_vertex_program_creater: public vertex_program_creater
{
    public:
        vertex_program::ptr create() const {
            return vertex_program::ptr(new closeness_vertex_program());
        }
};

void closeness_vertex::run(vertex_program &prog) {
    if ((vsize_t) prog.get_graph().get_curr_level() < dist) {
        //vertex_id_t id = prog.get_vertex_id(*this);
        //printf("Vertex: %u dropping to from: %u to: %u\n", id, dist,
                //prog.get_graph().get_curr_level());

        dist = prog.get_graph().get_curr_level();
        directed_vertex_request req(prog.get_vertex_id(*this), g_edge_type);
        request_partial_vertices(&req, 1);

        // Update per-thread cuml dist
        closeness_vertex_program& vprog = (closeness_vertex_program&) prog;
        vprog.cuml_dist_peq(dist);
    } // else do nothing
}

//// Helpers
//template <typename T>
//void p(std::vector<T> v) {

    //std::cout << "[ ";
    //for (const auto& _ : v)
        //std::cout << _ << " ";
    //std::cout << "]\n";
//}
}

class dist_reset: public vertex_initializer
{
    public:
        void init(compute_vertex &v) {
            closeness_vertex &dv = (closeness_vertex &) v;
            dv.reset();
        }
};

namespace fg
{

std::vector<double> compute_closeness_centrality(FG_graph::ptr fg,
        std::vector<vsize_t>& ids, edge_type traverse_e)
{

	if (!fg->get_graph_header().is_directed_graph()) {
        throw std::runtime_error(
                "This algorithm currently works on a directed graph\n");
	} else {
        g_edge_type = traverse_e;
    }

	graph_index::ptr index = NUMA_graph_index<closeness_vertex>::create(
			fg->get_graph_header());
	graph_engine::ptr graph = fg->create_engine(index);

	if (ids.size() == 0)
        for (vertex_id_t vid = 0; vid < fg->get_num_vertices(); vid++)
            ids.push_back(vid);

	struct timeval start, end;
	gettimeofday(&start, NULL);

    // Run the computation
    std::vector<double> res;

    // FIXME: Account for 0 and 1 degree nodes
    // TODO: Run multiple start_vertexs and batch them

    for (const auto& start_vertex : ids) {
        printf("Running vertex %u\n", start_vertex);

        // Reset to max_dist
        graph->init_all_vertices(vertex_initializer::ptr(
                    new dist_reset()));
        graph->start(&start_vertex, 1, vertex_initializer::ptr(),
                vertex_program_creater::ptr(
                    new closeness_vertex_program_creater()));
        graph->wait4complete();

        // Gather the cuml dist
        std::vector<vertex_program::ptr> progs;
        graph->get_vertex_programs(progs);

        // Reduction (+)
        size_t tot_dist = 0;
        for (auto prog : progs) {
            auto cvp = closeness_vertex_program::cast2(prog);
            tot_dist += cvp->get_cuml_dist();
        }

        if (tot_dist > 0)
            res.push_back(1.0/static_cast<double>(tot_dist));
        else
            res.push_back(0);
    }

	gettimeofday(&end, NULL);
    printf("Closeness took %.5f sec to complete\n", time_diff(start, end));

    assert(ids.size() == res.size());
    printf("Closeness values: \n");
    for (vsize_t id = 0; id < ids.size(); id++) {
        printf("idx: %u, val: %.5f\n", id, res[id]);
    }

	return res;
}
}
