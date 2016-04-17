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
#ifndef __BASE_KMEANS_COORDINATOR__
#define __BASE_KMEANS_COORDINATOR__

#include <vector>
#include <unordered_map>
#include <memory>
#include <atomic>

#include "libgraph-algs/kmeans_types.h"
#include "thread_state.h"

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

namespace km {
    class base_kmeans_coordinator {
        protected:
            unsigned nthreads, nnodes;
            size_t nrow, ncol;
            std::string fn; // File on disk
            unsigned* cluster_assignments;
            unsigned* cluster_assignment_counts;
            unsigned k;
            init_type_t _init_t;
            dist_type_t _dist_t;
            double tolerance;
            unsigned max_iters;
            unsigned num_changed; // total # samples changed in an iter
            std::atomic<unsigned> pending_threads; // How many threads have not completed their task

            // Threading
            pthread_mutex_t mutex;
            pthread_cond_t cond;
            pthread_mutexattr_t mutex_attr;

            base_kmeans_coordinator(const std::string fn, const size_t nrow,
                    const size_t ncol, const unsigned k, const unsigned max_iters,
                    const unsigned nnodes, const unsigned nthreads,
                    const double* centers, const init_type_t it,
                    const double tolerance, const dist_type_t dt);

        public:
            // Pass file handle to threads to read & numa alloc
            virtual void run_init() = 0;
            virtual void random_partition_init() = 0;
            virtual void forgy_init() = 0;

            virtual void run_kmeans() = 0;
            virtual void kmeanspp_init() = 0;
            virtual void wake4run(thread_state_t state) = 0;
            virtual const double* get_thd_data(const unsigned row_id) const = 0;

            virtual void set_thread_clust_idx(const unsigned clust_idx) = 0;
            virtual double reduction_on_cuml_sum() = 0;
            virtual void destroy_threads() = 0;
            virtual void set_thd_dist_v_ptr(double* v) = 0;

            void wait4complete();
    };
}
#endif
