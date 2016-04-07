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
#include "base_kmeans_thread.h"

namespace {
    class kmeans_thread : public base_kmeans_thread {
        private:
            clusters::ptr g_clusters; // Pointer to global cluster data
            unsigned nprocrows; // How many rows to process

            kmeans_thread(const int node_id, const unsigned thd_id, const unsigned start_rid,
                    const unsigned nprocrows, const unsigned ncol,
                    clusters::ptr g_clusters, unsigned* cluster_assignments,
                    const std::string fn) : base_kmeans_thread(node_id, thd_id, ncol,
                        g_clusters->get_nclust(), cluster_assignments, start_rid, fn) {

                this->nprocrows = nprocrows;
                this->g_clusters = g_clusters;

                set_data_size(sizeof(double)*nprocrows*ncol);
#if VERBOSE
                printf("Initializing thread. Metadata: thd_id: %u, node_id: %u"
                        ", start_rid: %lu, nprocrows: %u, ncol: %u\n",
                        this->thd_id, this->node_id, this->start_rid,
                        this->nprocrows, this->ncol);
#endif
            }

        public:
            typedef std::shared_ptr<kmeans_thread> ptr;

            static ptr create(const int node_id, const unsigned thd_id,
                    const unsigned start_rid, const unsigned nprocrows,
                    const unsigned ncol, clusters::ptr g_clusters,
                    unsigned* cluster_assignments, const std::string fn) {
                return ptr(new kmeans_thread(node_id, thd_id, start_rid,
                            nprocrows, ncol, g_clusters,
                            cluster_assignments, fn));
            }

            void start(const thread_state_t state);
            // Allocate and move data using this thread
            void EM_step();
            void kmspp_dist();
            const unsigned get_global_data_id(const unsigned row_id) const;
            void run();
            void wait();
            void wake(thread_state_t state);

            const void print_local_data() const;

            void destroy_numa_mem() {
                numa_free(local_data, get_data_size());
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
                kmspp_dist();
                break;
            case EM: /*E step of kmeans*/
                EM_step();
                break;
            case EXIT:
                fprintf(stderr, "[FATAL]: Thread state is EXIT but running!\n");
                exit(EXIT_FAILURE);
            default:
                fprintf(stderr, "[FATAL]: Unknown thread state\n");
                exit(EXIT_FAILURE);
        }

        int rc;

        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");

        (*parent_pending_threads)--;
        set_thread_state(WAIT);

        if (*parent_pending_threads == 0) {
            rc = pthread_cond_signal(parent_cond); // Wake up parent thread
            if (rc) perror("pthread_cond_signal");
        }
        pthread_mutex_unlock(&mutex);
    }

    void kmeans_thread::wait() {
        int rc;
        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");

        while (state == WAIT) {
            //printf("Thread %d begin cond_wait\n", thd_id);
            rc = pthread_cond_wait(&cond, &mutex);
            if (rc) perror("pthread_cond_wait");
        }

        rc = pthread_mutex_unlock(&mutex);
        if (rc) perror("pthread_mutex_unlock");
    }

    void kmeans_thread::wake(thread_state_t state) {
        int rc;
        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");
        set_thread_state(state);
        rc = pthread_mutex_unlock(&mutex);
        if (rc) perror("pthread_mutex_unlock");

        rc = pthread_cond_signal(&cond);
    }

    void* callback(void* arg) {
        kmeans_thread* t = static_cast<kmeans_thread*>(arg);
        bind2node_id(t->get_node_id());

        while (true) { // So we can receive task after task
            if (t->get_state() == WAIT)
                t->wait();

            if (t->get_state() == EXIT) {// No more work to do
                //printf("Thread %d exiting ...\n", t->thd_id);
                break;
            }

            //printf("Thread %d awake and doing a run()\n", t->thd_id);
            t->run(); // else
        }

        // We've stopped running so exit
        pthread_exit(NULL);
    }

    void kmeans_thread::start(const thread_state_t state=WAIT) {
        this->state = state;
        int rc = pthread_create(&hw_thd, NULL, callback, this);
        if (rc) {
            fprintf(stderr, "[FATAL]: Thread creation failed with code: %d\n", rc);
            exit(rc);
        }
    }

    const unsigned kmeans_thread::
        get_global_data_id(const unsigned row_id) const {
        return start_rid+row_id;
    }

    void kmeans_thread::EM_step() {
        meta.num_changed = 0; // Always reset at the beginning of an EM-step
        local_clusters->clear();

        for (unsigned row = 0; row < nprocrows; row++) {
            unsigned asgnd_clust = INVALID_CLUSTER_ID;
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
                meta.num_changed++;

            cluster_assignments[true_row_id] = asgnd_clust;
            local_clusters->add_member(&local_data[row*ncol], asgnd_clust);
        }
    }

    /** Method for a distance computation vs a single cluster.
      * Used in kmeans++ init
      */
    void kmeans_thread::kmspp_dist() {
        cuml_dist = 0;
        unsigned clust_idx = meta.clust_idx;
        for (unsigned row = 0; row < nprocrows; row++) {
            unsigned true_row_id = get_global_data_id(row);

            double dist = dist_comp_raw(&local_data[row*ncol],
                    &((g_clusters->get_means())[clust_idx*ncol]), ncol);

            if (dist < dist_v[true_row_id]) { // Found a closer cluster than before
                dist_v[true_row_id] = dist;
                cluster_assignments[true_row_id] = clust_idx;
            }
            cuml_dist += dist_v[true_row_id];
        }
    }

    const void kmeans_thread::print_local_data() const {
        print_mat(local_data, nprocrows, ncol);
    }
}
#endif
