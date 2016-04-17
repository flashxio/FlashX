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
#include <atomic>
#include <boost/assert.hpp>

#include "log.h"
#include "common.h"
#include "thread_state.h"

#define VERBOSE 0
#define INVALID_THD_ID -1

namespace km {
    class clusters;
}

namespace km {
    union metaunion {
        unsigned num_changed; // Used during kmeans
        unsigned clust_idx; // Used during kms++
    };

    class base_kmeans_thread {
        protected:
            pthread_t hw_thd;
            unsigned node_id; // Which NUMA node are you on?
            int thd_id;
            unsigned start_rid; // With respect to the original data
            unsigned ncol; // How many columns in the data
            double* local_data; // Pointer to where the data begins that the thread works on
            size_t data_size; // true size of local_data at any point
            std::shared_ptr<km::clusters> local_clusters;

            pthread_mutex_t mutex;
            pthread_cond_t cond;
            pthread_mutexattr_t mutex_attr;

            pthread_cond_t* parent_cond;
            std::atomic<unsigned>* parent_pending_threads;

            metaunion meta;
            //unsigned num_changed;

            FILE* f; // Data file on disk
            unsigned* cluster_assignments;

            thread_state_t state;
            double* dist_v;
            double cuml_dist;

            friend void* callback(void* arg);

            base_kmeans_thread(const int node_id, const unsigned thd_id, const unsigned ncol,
                    const unsigned nclust, unsigned* cluster_assignments, const unsigned start_rid,
                    const std::string fn) {


                pthread_mutexattr_init(&mutex_attr);
                pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
                pthread_mutex_init(&mutex, &mutex_attr);
                pthread_cond_init(&cond, NULL);
                this->node_id = node_id;
                this->thd_id = thd_id;
                this->ncol = ncol;
                this->cluster_assignments = cluster_assignments;
                this->start_rid = start_rid;
                BOOST_VERIFY(this->f = fopen(fn.c_str(), "rb"));

                meta.num_changed = 0; // Same as meta.clust_idx = 0;
                set_thread_state(WAIT);
            }

        protected:
            void set_thread_state(thread_state_t state) {
                this->state = state;
            }

        public:
            typedef std::shared_ptr<base_kmeans_thread> ptr;

            virtual void start(const thread_state_t state) = 0;
            // Allocate and move data using this thread
            virtual void EM_step() = 0;
            virtual void kmspp_dist() = 0;
            virtual const unsigned get_global_data_id(const unsigned row_id) const = 0;
            virtual void run() = 0;
            virtual void sleep() = 0;

            void numa_alloc_mem();

            void test() {
                //printf("%u ", get_thd_id());
            }

            void set_dist_v_ptr(double* v) {
                dist_v = v;
            }

            const thread_state_t get_state() const {
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

            const std::shared_ptr<km::clusters> get_local_clusters() const {
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

            pthread_mutex_t& get_lock() {
                return mutex;
            }

            pthread_cond_t& get_cond() {
                return cond;
            }

            unsigned get_node_id() {
                return node_id;
            }

            void set_parent_cond(pthread_cond_t* cond) {
                parent_cond = cond;
            }

            void set_parent_pending_threads(std::atomic<unsigned>* ppt) {
                parent_pending_threads = ppt;
            }

            void destroy_numa_mem() {
                numa_free(local_data, get_data_size());
            }

            void join();
            void close_file_handle();
            ~base_kmeans_thread();
    };

    void base_kmeans_thread::join() {
        void* join_status;
        int rc = pthread_join(hw_thd, &join_status);
        if (rc) {
            fprintf(stderr, "[FATAL]: Return code from pthread_join() "
                    "is %d\n", rc);
            exit(rc);
        }
        thd_id = INVALID_THD_ID;
    }

    // Once the algorithm ends we should deallocate the memory we moved
    void base_kmeans_thread::close_file_handle() {
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

    base_kmeans_thread::~base_kmeans_thread() {
        pthread_cond_destroy(&cond);
        pthread_mutex_destroy(&mutex);
        pthread_mutexattr_destroy(&mutex_attr);

        if (f)
            close_file_handle();
#if VERBOSE
        printf("Thread %u being destroyed\n", thd_id);
#endif
        if (thd_id != INVALID_THD_ID)
            join();
    }

    void bind2node_id(int node_id) {
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
        fseek(f, start_rid*ncol*sizeof(double), SEEK_CUR); // start position
        BOOST_VERIFY(1 == fread(local_data, blob_size, 1, f));
        close_file_handle();
    }
}
#endif
