/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#ifndef __SEM_KMEANS_H__
#define __SEM_KMEANS_H__

#include <math.h>

#include <vector>
#include <algorithm>
#include <map>
#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include "graph_engine.h"
#include "graph_config.h"
#include "FGlib.h"
#include "save_result.h"

#include "sem_kmeans_util.h"
#include "prune_stats.h"
#include "dist_matrix.h"
#include "kmeans_types.h"

#define PAGE_ROW
#define KM_TEST 1
#define VERBOSE 0

using namespace fg;
using namespace km;

namespace {
    typedef safs::page_byte_array::seq_const_iterator<double> data_seq_iter;
    enum kmspp_stage_t { ADDMEAN, DIST }; // Either adding a mean / computing dist

    static const unsigned INVALID_CLUST_ID = -1;
    static unsigned NUM_COLS;
    static unsigned K;
    static std::map<vertex_id_t, unsigned> g_init_hash; // Used for forgy init
    static std::vector<double> g_kmspp_distance; // Used for kmeans++ init

    class base_kmeans_vertex: public compute_vertex
    {
        unsigned cluster_id;

        public:
        base_kmeans_vertex(vertex_id_t id):
            compute_vertex(id) {
                cluster_id = INVALID_CLUST_ID;
            }

        unsigned get_result() const {
            return cluster_id;
        }

        const unsigned get_cluster_id() const {
            return cluster_id;
        }

        void set_cluster_id(const unsigned id ) { cluster_id = id; }
        void run(vertex_program &prog) {
            BOOST_ASSERT_MSG(0, "base_kmeans_vertex run(vertex_program&) must not be called!");
        }
        void run(vertex_program& prog, const page_vertex &vertex) {
            BOOST_ASSERT_MSG(0, "base_kmeans_vertex run(vertex_program&,"
                "const page_vertex&) must not be called!");
        }
        void run_on_message(vertex_program& prog, const vertex_message& msg) { }
    };

    /* Used in per thread cluster formation */
    template <typename T, typename ClusterType>
    class base_kmeans_vertex_program: public vertex_program_impl<T>
    {
        unsigned pt_changed;
        typename ClusterType::ptr pt_clusters;

        public:
        typedef std::shared_ptr<base_kmeans_vertex_program<T, ClusterType> > ptr;

        //TODO: Opt only add cluster when a vertex joins it
        base_kmeans_vertex_program() {
            this->pt_changed = 0;

            pt_clusters = ClusterType::create(K, NUM_COLS);
        }

        static ptr cast2(vertex_program::ptr prog) {
            return std::static_pointer_cast<base_kmeans_vertex_program<T, ClusterType>,
                   vertex_program>(prog);
        }

        typename ClusterType::ptr get_pt_clusters() {
            return pt_clusters;
        }

        void add_member(const unsigned id, data_seq_iter& count_it) {
            pt_clusters->add_member(count_it, id);
        }

        void add_member(const unsigned id, const double* row) {
            pt_clusters->add_member(row, id);
        }

        const unsigned get_pt_changed() { return pt_changed; }

        void pt_changed_pp() {
            pt_changed++;
        }
    };

    // Begin helpers
    void print_sample(vertex_id_t my_id, data_seq_iter& count_it) {
        std::vector<std::string> v;
        while (count_it.has_next()) {
            char buffer [1024];
            double e = count_it.next();
            assert(sprintf(buffer, "%e", e));
            v.push_back(std::string(buffer));
        }
        printf("V%u's vector: \n", my_id); print_vector<std::string>(v);
    }

    std::string s (const double d) {
        if (d == std::numeric_limits<double>::max())
            return "max";
        else
            return std::to_string(d);
    }
    // End helpers //
}

namespace fg
{
    // A class used a return object for R bindings
    class sem_kmeans_ret
    {
        private:
            FG_vector<unsigned>::ptr cluster_assignments;
            std::vector<std::vector<double>> centers;
            std::vector<unsigned> size;
            unsigned iters;

            sem_kmeans_ret(const FG_vector<unsigned>::ptr cluster_assignments,
                    const std::vector<std::vector<double>> centers,
                    const std::vector<unsigned>& size, const unsigned iters) {
                this->cluster_assignments = cluster_assignments;
                this->centers = centers;
                this->size = size;
                this->iters = iters;
            }

        public:
            typedef typename std::shared_ptr<sem_kmeans_ret> ptr;

            static ptr create(const FG_vector<unsigned>::ptr cluster_assignments,
                    const std::vector<std::vector<double>> centers,
                    const std::vector<unsigned>& size, const unsigned iters) {
                return ptr(new sem_kmeans_ret(cluster_assignments, centers, size, iters));
            }

            const FG_vector<unsigned>::ptr get_cluster_assignments() const {
                return this->cluster_assignments;
            }

            const unsigned get_iters() const {
                return this->iters;
            }

            const std::vector<unsigned>& get_size() const {
                return this->size;
            }

            const std::vector<std::vector<double>>& get_centers() const {
                return this->centers;
            }
    };

    /**
     * \brief Compute Semi-External Memory kmeans
     * \param fg The FlashGraph graph object for which you want to compute.
     * \param k The number of clusters.
     * \param init Initialization type [random, forgy, kmeanspp].
     * \param max_iters The max number of iterations to compute for.
     * \param tolerance The min fraction of changes from 1 iter to next required to converge.
     * \param num_rows The # of rows in the dataset
     * \param num_cols The # of columns in the dataset
     * \param centers To skip any initialization use this to pass initialized centers
     */
    sem_kmeans_ret::ptr compute_sem_kmeans(FG_graph::ptr fg, const unsigned k, const std::string init,
            const unsigned max_iters, const double tolerance, const unsigned num_rows=0,
            const unsigned num_cols=0, std::vector<double>* centers=NULL);

    /**
     * \brief Compute Semi-External Memory Triangle Inequality pruned kmeans
     * see: `compute_sem_kmeans` for description of parameter list
     * *SEE: http://users.cecs.anu.edu.au/~daa/courses/GSAC6017/kmeansicml03.pdf
     */
    sem_kmeans_ret::ptr compute_triangle_sem_kmeans(FG_graph::ptr fg, const unsigned k, const std::string init,
            const unsigned max_iters, const double tolerance, const unsigned num_rows=0,
            const unsigned num_cols=0, std::vector<double>* centers=NULL);

    /**
     * \brief Minimized version of Semi-External Memory Triangle Inequality pruned kmeans
     * see: `compute_sem_kmeans` for description of parameter list
     * *SEE: http://users.cecs.anu.edu.au/~daa/courses/GSAC6017/kmeansicml03.pdf
     * \param cache_size_gb the size of row cache in gb
     */
    sem_kmeans_ret::ptr compute_min_triangle_sem_kmeans(FG_graph::ptr fg, const unsigned k,
            const std::string init, const unsigned max_iters, const double tolerance,
            const unsigned num_rows, const unsigned num_cols,
            std::vector<double>* centers=NULL, const double cache_size_gb=0,
            const unsigned rc_update_start_interval=5);
}
#endif
