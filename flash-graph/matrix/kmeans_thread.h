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

#ifndef __KMEANS_THREAD_H__
#define __KMEANS_THREAD_H__

#include <pthread.h>
#include <numa.h>

#include <memory>
#include <utility>
#include <boost/format.hpp>

#include "log.h"
#include "libgraph-algs/clusters.h"
#include "kmeans.h"

namespace {
    enum thread_state_t {
        TEST, /*just for testing*/
        ALLOC_DATA, /*moving data for reduces rma*/
        KMSPP_INIT,
        EM /*EM steps of kmeans*/
    };

    class kmeans_thread {
        private:
            pthread_t hw_thd;
            unsigned node_id; // Which NUMA node are you on?
            unsigned thd_id;
            size_t start_offset; // With respect to the original data
            unsigned nprocrows; // How many rows to process
            unsigned ncol; // How many columns in the data
            double* local_data; // Pointer to where the data begins that the thread works on
            clusters::ptr g_clusters; // Pointer to global cluster data
            clusters::ptr local_clusters;
            std::vector<size_t> local_num_changed;
            FILE* f; // Data file on disk

            unsigned* cluster_assignments;

            thread_state_t state;

            friend void* callback(void* arg);
            void run();
            kmeans_thread(const int node_id, const unsigned thd_id, const size_t start_offset,
                    const unsigned nprocrows, const unsigned ncol,
                    clusters::ptr g_clusters, unsigned* cluster_assignments,
                    const std::string fn) {
                this->node_id = node_id;
                this->thd_id = thd_id;
                this->start_offset = start_offset;
                this->nprocrows = nprocrows;
                this->ncol = ncol;
                this->g_clusters = g_clusters;
                this->cluster_assignments = cluster_assignments;
                BOOST_VERIFY(this->f = fopen(fn.c_str(), "rb"));

                printf("Starting thread. Metadata: thd_id: %u, node_id: %u"
                        ", start_offset: %lu, nprocrows: %u, ncol: %u\n", this->thd_id,
                        this->node_id, this->start_offset, this->nprocrows, this->ncol);
            }
            // Testing
            unsigned val;
        public:
            typedef std::shared_ptr<kmeans_thread> ptr;

            static ptr create(const int node_id, const unsigned thd_id,
                    const size_t start_offset, const unsigned nprocrows,
                    const unsigned ncol, clusters::ptr g_clusters,
                    unsigned* cluster_assignments, const std::string fn) {
                return ptr(new kmeans_thread(node_id, thd_id, start_offset,
                            nprocrows, ncol, g_clusters, cluster_assignments, fn));
            }

            void join() {
                void* join_status;
                printf("Thread %u calling join()\n", thd_id);

                int rc = pthread_join(hw_thd, &join_status);
                if (rc) {
                    fprintf(stderr, "[FATAL]: Return code from pthread_join() "
                            "is %d\n", rc);
                    exit(rc);
                }
                printf("Sucessful join on thread: %u with status %lu\n", thd_id,
                        (long)join_status);
            }

            const unsigned get_thd_id() const {
                return thd_id;
            }

            // Testing
            void set_val(const unsigned val) { this->val = val; }
            const unsigned get_val() const { return val; }
            void test();
            // END Testing

            void start(const thread_state_t state);
            // Allocate and move data using this thread
            void numa_alloc_mem();
            void print_local_data();

            void destroy_numa_mem() {
                numa_free(local_data, get_data_size());
            }

            size_t get_data_size() {
                return sizeof(double)*nprocrows*ncol;
            }

            // Once the algorithm ends we should deallocate the memory we moved
            void close_file_handle() {
                int rc = fclose(f);
                if (rc) {
                    fprintf(stderr, "[FATAL]: fclose() failed with code: %d\n", rc);
                    exit(rc);
                }
                f = NULL;
            }

            ~kmeans_thread() {
                if (f)
                    close_file_handle();
                printf("Thread %u being destroyed\n", thd_id);
            }
    };
}
#endif
