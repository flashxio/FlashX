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
        EM, /*EM steps of kmeans*/
        PAUSED /*When the thread is waiting to be run or killed*/
    };

    class kmeans_thread {
        private:
            pthread_t hw_thd;
            unsigned node_id; // Which NUMA node are you on?
            unsigned thd_id;
            size_t start_offset; // With respect to the original data
            unsigned nprocrows; // How many rows to process
            unsigned thds_row; // How many rows to process
            unsigned ncol; // How many columns in the data
            double* local_data; // Pointer to where the data begins that the thread works on
            clusters::ptr g_clusters; // Pointer to global cluster data
            clusters::ptr local_clusters;
            unsigned num_changed;
            FILE* f; // Data file on disk

            unsigned* cluster_assignments;

            thread_state_t state;

            friend void* callback(void* arg);
            void run();
            kmeans_thread(const int node_id, const unsigned thd_id, const size_t start_offset,
                    const unsigned nprocrows, const unsigned thds_row, const unsigned ncol,
                    clusters::ptr g_clusters, unsigned* cluster_assignments,
                    const std::string fn) {
                this->node_id = node_id;
                this->thd_id = thd_id;
                this->start_offset = start_offset;
                this->nprocrows = nprocrows;
                this->thds_row = thds_row; // Threads per row other than possibly the last one
                this->ncol = ncol;
                this->g_clusters = g_clusters;
                this->cluster_assignments = cluster_assignments;
                BOOST_VERIFY(this->f = fopen(fn.c_str(), "rb"));

                local_clusters = clusters::create(g_clusters->get_nclust(), ncol);
                this->num_changed = 0;

#if KM_TEST
                printf("Initializing thread. Metadata: thd_id: %u, node_id: %u"
                        ", start_offset: %lu, nprocrows: %u, ncol: %u\n", this->thd_id,
                        this->node_id, this->start_offset, this->nprocrows, this->ncol);
#endif
            }
#if KM_TEST
            unsigned val;
#endif
        public:
            typedef std::shared_ptr<kmeans_thread> ptr;

            static ptr create(const int node_id, const unsigned thd_id,
                    const size_t start_offset, const unsigned nprocrows,
                    const unsigned thds_row, const unsigned ncol,
                    clusters::ptr g_clusters, unsigned* cluster_assignments,
                    const std::string fn) {
                return ptr(new kmeans_thread(node_id, thd_id, start_offset,
                            nprocrows, thds_row, ncol, g_clusters,
                            cluster_assignments, fn));
            }

            void join() {
                void* join_status;
                int rc = pthread_join(hw_thd, &join_status);
                if (rc) {
                    fprintf(stderr, "[FATAL]: Return code from pthread_join() "
                            "is %d\n", rc);
                    exit(rc);
                }
                this->state = PAUSED;
            }

            const unsigned get_thd_id() const {
                return thd_id;
            }

#if KM_TEST
            void set_val(const unsigned val) { this->val = val; }
            const unsigned get_val() const { return val; }
            void test();
#endif

            void start(const thread_state_t state);
            // Allocate and move data using this thread
            void numa_alloc_mem();
            void EM_step();

            const void print_local_data() const;
            const unsigned get_global_data_id(const unsigned row_id) const;

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
                printf("Thread %u closing the file handle.\n",thd_id);
                f = NULL;
            }

            const double* get_local_data() const {
                return local_data;
            }

            const unsigned get_num_changed() const {
                return num_changed;
            }

            const clusters::ptr get_local_clusters() const {
                return local_clusters;
            }

            ~kmeans_thread() {
                if (f)
                    close_file_handle();
                printf("Thread %u being destroyed\n", thd_id);
            }
    };

    void kmeans_thread::run() {
        switch(state) {
            case TEST:
                test();
                break;
            case ALLOC_DATA:
                numa_alloc_mem();
                break;
            case KMSPP_INIT:
                assert(0);
                break;
            case EM: /*E step of kmeans*/
                EM_step();
                break;
            case PAUSED:
                fprintf(stderr, "[FATAL]: Thread state is PAUSED but running!\n");
                exit(EXIT_FAILURE);
            default:
                fprintf(stderr, "Unknown thread state\n");
                exit(EXIT_FAILURE);
        }
    }

    void kmeans_thread::test() {
        printf("Hello from thread: %u\n", get_thd_id());
        set_val(get_thd_id());
    }

    static void bind2node_id(int node_id) {
        struct bitmask *bmp = numa_allocate_nodemask();
        numa_bitmask_setbit(bmp, node_id);
        numa_bind(bmp);
        numa_free_nodemask(bmp);
    }

    void* callback(void* arg) {
        kmeans_thread* t = static_cast<kmeans_thread*>(arg);
        bind2node_id(t->node_id);
        t->run();
        pthread_exit(NULL);
    }

    void kmeans_thread::start(const thread_state_t state) {
        this->state = state;
        int rc = pthread_create(&hw_thd, NULL, callback, this);
        if (rc) {
            fprintf(stderr, "[FATAL]: Thread creation failed with code: %d\n", rc);
            exit(rc);
        }
    }

    // Move data ~equally to all nodes
    void kmeans_thread::numa_alloc_mem() {
        BOOST_ASSERT_MSG(f, "File handle invalid, can only alloc once!");
        size_t blob_size = get_data_size();
        local_data = static_cast<double*>(numa_alloc_onnode(blob_size, node_id));
        fseek(f, start_offset*sizeof(double), SEEK_CUR); // start position
        BOOST_VERIFY(1 == fread(local_data, blob_size, 1, f));
        close_file_handle();
    }

    /*TODO: Check if it's cheaper to do per thread `cluster_assignments`*/
    const unsigned kmeans_thread::
        get_global_data_id(const unsigned row_id) const {
        return (thd_id*thds_row)+row_id;
    }

    void kmeans_thread::EM_step() {
        num_changed = 0; // Always reset at the beginning of an EM-step
        local_clusters->clear();

        for (unsigned row = 0; row < nprocrows; row++) {
            size_t asgnd_clust = INVALID_CLUSTER_ID;
            double best, dist;
            dist = best = std::numeric_limits<double>::max();

            for (unsigned clust_idx = 0;
                    clust_idx < g_clusters->get_nclust(); clust_idx++) {
                dist = dist_comp_raw(&local_data[row*ncol],
                        &(g_clusters->get_means()[clust_idx*ncol]), ncol);

                if (dist < best) {
                    best = dist;
                    asgnd_clust = clust_idx;
                }
            }

            BOOST_VERIFY(asgnd_clust != INVALID_CLUSTER_ID);
            unsigned true_row_id = get_global_data_id(row);

            if (asgnd_clust != cluster_assignments[true_row_id])
                num_changed++;

            cluster_assignments[true_row_id] = asgnd_clust;
            local_clusters->add_member(&local_data[row*ncol], asgnd_clust);
        }
    }

    const void kmeans_thread::print_local_data() const {
        print_mat(local_data, nprocrows, ncol);
    }
}
#endif
