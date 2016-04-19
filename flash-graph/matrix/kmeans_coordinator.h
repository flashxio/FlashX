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
#ifndef __KMEANS_COORDINATOR_H__
#define __KMEANS_COORDINATOR_H__

#include <vector>
#include <unordered_map>
#include <memory>

#include "log.h"
#include "base_kmeans_coordinator.h"
#include "libgraph-algs/kmeans_types.h"
#include "thread_state.h"

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

namespace km {
    class clusters;
}

class kmeans_thread;

namespace km {
    typedef std::vector<std::shared_ptr<kmeans_thread> >::iterator thread_iter;

    class kmeans_coordinator : public km::base_kmeans_coordinator {
        private:
            std::vector<std::shared_ptr<kmeans_thread> > threads;
            // Metadata
            // max index stored within each threads partition
            std::vector<unsigned> thd_max_row_idx;
            std::shared_ptr<km::clusters> cltrs;

            kmeans_coordinator(const std::string fn, const size_t nrow,
                    const size_t ncol, const unsigned k, const unsigned max_iters,
                    const unsigned nnodes, const unsigned nthreads,
                    const double* centers, const km::init_type_t it,
                    const double tolerance, const km::dist_type_t dt);

        public:
            typedef std::shared_ptr<kmeans_coordinator> ptr;

            static ptr create(const std::string fn, const size_t nrow,
                    const size_t ncol, const unsigned k, const unsigned max_iters,
                    const unsigned nnodes, const unsigned nthreads,
                    const double* centers=NULL, const std::string init="kmeanspp",
                    const double tolerance=-1, const std::string dist_type="eucl") {

                km::init_type_t _init_t;
                if (init == "random")
                    _init_t = km::init_type_t::RANDOM;
                else if (init == "forgy")
                    _init_t = km::init_type_t::FORGY;
                else if (init == "kmeanspp")
                    _init_t = km::init_type_t::PLUSPLUS;
                else if (init == "none")
                    _init_t = km::init_type_t::NONE;
                else {
                    BOOST_LOG_TRIVIAL(fatal) << "[ERROR]: param init must be one of:"
                       " [random | forgy | kmeanspp]. It is '" << init << "'";
                    exit(-1);
                }

                km::dist_type_t _dist_t;
                if (dist_type == "eucl")
                    _dist_t = km::dist_type_t::EUCL;
                else if (dist_type == "cos")
                    _dist_t = km::dist_type_t::COS;
                else {
                    BOOST_LOG_TRIVIAL(fatal) << "[ERROR]: param dist_type must be one of:"
                       " 'eucl', 'cos'.It is '" << dist_type << "'";
                    exit(-1);
                }
#if KM_TEST
                printf("kmeans coordinator => NUMA nodes: %u, nthreads: %u, "
                        "nrow: %lu, ncol: %lu, init: '%s', dist_t: '%s', fn: '%s'"
                        "\n\n", nnodes, nthreads, nrow, ncol, init.c_str(),
                        dist_type.c_str(), fn.c_str());
#endif
                return ptr(new kmeans_coordinator(fn, nrow, ncol, k, max_iters,
                            nnodes, nthreads, centers, _init_t, tolerance, _dist_t));
            }

            std::pair<unsigned, unsigned> get_rid_len_tup(const unsigned thd_id);
            // Pass file handle to threads to read & numa alloc
            void create_thread_map();
            void run_kmeans();
            void update_clusters();
            void kmeanspp_init();
            void wake4run(km::thread_state_t state);
            void destroy_threads();
            void set_thread_clust_idx(const unsigned clust_idx);
            double reduction_on_cuml_sum();
            void set_thd_dist_v_ptr(double* v);
            void run_init();
            void random_partition_init();
            void forgy_init();
            const double* get_thd_data(const unsigned row_id) const;
            ~kmeans_coordinator();
    };
}
#endif
