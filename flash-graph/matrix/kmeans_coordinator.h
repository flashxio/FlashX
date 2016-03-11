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
#ifndef __KMEANS_COORDINATOR__
#define __KMEANS_COORDINATOR__

#include <vector>
#include <unordered_map>
#include <memory>

#include "kmeans_thread.h"

namespace {
    typedef std::vector<kmeans_thread::ptr>::iterator thread_iter;

    class kmeans_coordinator {
        private:
            unsigned nthreads, nnodes, nrow, ncol;
            std::vector<kmeans_thread::ptr> threads; 
            std::string fn; // File on disk
            unsigned* cluster_assignments;
            unsigned* cluster_assignment_counts;
            unsigned k;
            clusters::ptr g_clusters;


            kmeans_coordinator(const unsigned nthreads, const unsigned nnodes,
                    const unsigned nrow, const unsigned ncol, const unsigned k,
                    const std::string fn) {
                this->nthreads = nthreads;
                this->nnodes = nnodes;
                this->nrow = nrow;
                this->ncol = ncol;
                this->fn = fn;
                this->k = k;

                BOOST_VERIFY(cluster_assignments = new unsigned [nrow]);
                BOOST_VERIFY(cluster_assignment_counts = new unsigned [k]);

                std::fill(&cluster_assignments[0],
                        (&cluster_assignments[0])+nrow, -1);
                std::fill(&cluster_assignment_counts[0],
                        (&cluster_assignment_counts[0])+k, 0);

                 clusters::ptr g_clusters = clusters::create(k, ncol);

                // NUMA node affinity binding policy is round-robin
                for (unsigned thd_id = 0; thd_id < nthreads; thd_id++) {
                    std::pair<size_t, unsigned> tup = get_offset_len_tup(thd_id);
                    threads.push_back(kmeans_thread::create((thd_id % nnodes),
                                thd_id, tup.first, tup.second, ncol, g_clusters,
                                cluster_assignments, fn));
                }
            }

        public:
            typedef std::shared_ptr<kmeans_coordinator> ptr;
            static ptr create(const unsigned nthreads, const unsigned nnodes,
                    const unsigned nrow, const unsigned ncol, const unsigned k,
                    const std::string fn) {
                return ptr(new kmeans_coordinator(nthreads, nnodes,
                            nrow, ncol, k, fn));
            }

            std::pair<size_t, unsigned> get_offset_len_tup(const unsigned thd_id) {
                unsigned rows_per_thread = nrow / nthreads;
                size_t start_offset = (thd_id*rows_per_thread*ncol);

                if (thd_id == nthreads - 1)
                    rows_per_thread += nrow % nthreads;
                return std::pair<size_t, unsigned>(start_offset, rows_per_thread);
            }

            // Pass file handle to threads to read & numa alloc
            void numa_alloc_data();
            void join_threads();

            ~kmeans_coordinator() {
                std::vector<kmeans_thread::ptr>::iterator it = threads.begin();
                for (; it != threads.end(); ++it)
                    (*it)->destroy_numa_mem();

                delete [] cluster_assignments;
                delete [] cluster_assignment_counts;
            }

            void create_thread_map();
    };
}
#endif
