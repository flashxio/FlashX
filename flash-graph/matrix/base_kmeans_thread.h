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

#ifndef __BASE_KMEANS_THREAD_H__
#define __BASE_KMEANS_THREAD_H__

#include <pthread.h>
#include <numa.h>

#include <memory>
#include <utility>
#include <boost/format.hpp>

#include "log.h"
#include "libgraph-algs/clusters.h"
#include "kmeans.h"

#define KM_TEST 1
#define VERBOSE 0

namespace {
    enum thread_state_t {
        TEST, /*just for testing*/
        ALLOC_DATA, /*moving data for reduces rma*/
        KMSPP_INIT,
        EM, /*EM steps of kmeans*/
        IDLE /*When the thread is waiting to be run or killed*/
    };

    union metaunion {
        unsigned num_changed; // Used during kmeans
        unsigned clust_idx; // Used during kms++
    };

    class base_kmeans_thread {
        friend class kmeans_thread; // FIXME: RM
        private:
            pthread_t hw_thd;
            unsigned node_id; // Which NUMA node are you on?
            unsigned thd_id;
            size_t start_offset; // With respect to the original data
            unsigned ncol; // How many columns in the data
            double* local_data; // Pointer to where the data begins that the thread works on
            size_t data_size; // true size of local_data at any point
            clusters::ptr local_clusters;

            metaunion meta;
            //unsigned num_changed;

            FILE* f; // Data file on disk

            unsigned* cluster_assignments;

            thread_state_t state;
            double* dist_v;
            double cuml_dist;

            friend void* callback(void* arg);

            base_kmeans_thread(const int node_id, const unsigned thd_id, const unsigned ncol,
                    const unsigned nclust, unsigned* cluster_assignments, size_t start_offset,
                    const std::string fn) {
                this->node_id = node_id;
                this->thd_id = thd_id;
                this->ncol = ncol;
                this->cluster_assignments = cluster_assignments;
                this->start_offset = start_offset;
                BOOST_VERIFY(this->f = fopen(fn.c_str(), "rb"));

                local_clusters = clusters::create(nclust, ncol);
                meta.num_changed = 0; // Same as meta.clust_idx = 0;
            }

        public:

            virtual void start(const thread_state_t state) = 0;
            // Allocate and move data using this thread
            virtual void EM_step() = 0;
            virtual void kmspp_dist() = 0;
            virtual const unsigned get_global_data_id(const unsigned row_id) const = 0;
            virtual void run() = 0;
            virtual void destroy_numa_mem() = 0;

            void numa_alloc_mem();
            void test();

            void set_dist_v_ptr(double* v) {
                dist_v = v;
            }

            void join() {
                void* join_status;
                int rc = pthread_join(hw_thd, &join_status);
                if (rc) {
                    fprintf(stderr, "[FATAL]: Return code from pthread_join() "
                            "is %d\n", rc);
                    exit(rc);
                }
                this->state = IDLE;
            }

            const thread_state_t get_state() {
                return this->state;
            }

            const unsigned get_thd_id() const {
                return thd_id;
            }

            const double* get_local_data() const {
                return local_data;
            }

            const unsigned get_num_changed() const {
                return meta.num_changed;
            }

            const clusters::ptr get_local_clusters() const {
                return local_clusters;
            }

            void set_clust_idx(const unsigned idx) {
                meta.clust_idx = idx;
            }

            const double get_culm_dist() const {
                return cuml_dist;
            }

            void set_data_size(const size_t data_size) {
                this->data_size = data_size;
            }

            const size_t get_data_size() const {
                return this->data_size;
            }

            //FIXME: bunch of getters and setters
            // Once the algorithm ends we should deallocate the memory we moved
            void close_file_handle() {
                int rc = fclose(f);
                if (rc) {
                    fprintf(stderr, "[FATAL]: fclose() failed with code: %d\n", rc);
                    exit(rc);
                }
#if VERBOSE
                printf("Thread %u closing the file handle.\n",thd_id);
#endif
                f = NULL;
            }

            ~base_kmeans_thread() {
                if (f)
                    close_file_handle();
#if VERBOSE
                printf("Thread %u being destroyed\n", thd_id);
#endif
            }
    };

    void base_kmeans_thread::test() {
        printf("Hello from thread: %u\n", get_thd_id());
    }

    static void bind2node_id(int node_id) {
        struct bitmask *bmp = numa_allocate_nodemask();
        numa_bitmask_setbit(bmp, node_id);
        numa_bind(bmp);
        numa_free_nodemask(bmp);
    }

    // Move data ~equally to all nodes
    void base_kmeans_thread::numa_alloc_mem() {
        BOOST_ASSERT_MSG(f, "File handle invalid, can only alloc once!");
        size_t blob_size = get_data_size();
        local_data = static_cast<double*>(numa_alloc_onnode(blob_size, node_id));
        fseek(f, start_offset*sizeof(double), SEEK_CUR); // start position
        BOOST_VERIFY(1 == fread(local_data, blob_size, 1, f));
        close_file_handle();
    }
}
#endif
