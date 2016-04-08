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

#ifndef __KMEANS_TASK_THREAD_H__
#define __KMEANS_TASK_THREAD_H__

#include "base_kmeans_thread.h"
#include "task_queue.h"
#include "libgraph-algs/dist_matrix.h"
#include "thd_safe_bool_vector.h"
#include <atomic>

using namespace km;

namespace prune {
    class kmeans_task_thread : public base_kmeans_thread {
        private:
            prune_clusters::ptr g_clusters; // Pointer to global cluster data
            unsigned start_rid; // The row id of the first item in this partition

            std::vector<kmeans_task_thread::ptr>* workers; // all other worker threads
            task_queue tasks;
            task curr_task;

            bool prune_init;
            prune::dist_matrix::ptr dm; // global
            thd_safe_bool_vector::ptr recalculated_v; // global

            kmeans_task_thread(const int node_id, const unsigned thd_id,
                    const unsigned start_rid, const unsigned nlocal_rows,
                    const unsigned ncol, prune_clusters::ptr g_clusters,
                    unsigned* cluster_assignments,
                    const std::string fn) : base_kmeans_thread(node_id, thd_id, ncol,
                        g_clusters->get_nclust(), cluster_assignments, start_rid, fn) {

                this->g_clusters = g_clusters;
                // Init task queue
                tasks.set_start_rid(start_rid);
                tasks.set_nrow(nlocal_rows);
                tasks.set_ncol(ncol);
                prune_init = true;

                set_data_size(sizeof(double)*nlocal_rows*ncol);
#if VERBOSE
                printf("Initializing thread. Metadata: thd_id: %u, node_id: %u"
                        ", start_rid: %lu, nlocal_rows: %u, ncol: %u\n",
                        this->thd_id, this->node_id, this->start_rid,
                        nlocal_rows, this->ncol);
#endif
            }

        public:
            typedef std::shared_ptr<kmeans_task_thread> ptr;

            static ptr create(const int node_id, const unsigned thd_id,
                    const unsigned start_rid, const unsigned nlocal_rows,
                    const unsigned ncol, prune_clusters::ptr g_clusters,
                    unsigned* cluster_assignments, const std::string fn) {
                return ptr(new kmeans_task_thread(node_id, thd_id, start_rid,
                            nlocal_rows, ncol, g_clusters,
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
            void request_task();
            void lock_sleep();
            void sleep();

            const void print_local_data() const;
            double* get_dist_v_ptr() { return &dist_v[0]; }

            void set_parent_cond(pthread_cond_t* cond) {
                parent_cond = cond;
            }

            void set_parent_pending_threads(std::atomic<unsigned>* ppt) {
                parent_pending_threads = ppt;
            }

            void set_prune_init(const bool prune_init) {
                this->prune_init = prune_init;
            }

            const bool is_prune_init() {
                return prune_init;
            }

            void set_recalc_v_ptr(thd_safe_bool_vector::ptr recalculated_v) {
                this->recalculated_v = recalculated_v;
            }

            void set_dist_mat_ptr(prune::dist_matrix::ptr dm) {
                this->dm = dm;
            }
    };

    void kmeans_task_thread::request_task() {
        int rc;
        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");

        if (tasks.has_task()) {
            // Grab another task
            //delete curr_task;
            curr_task = tasks.get_task();
            BOOST_VERIFY(curr_task.get_nrow() <= tasks.get_nrow());

            // FIXME: someone got the last task
            //printf("request_task: Thd: %u, Task ==> ", get_thd_id()); curr_task.print();
            BOOST_ASSERT_MSG(curr_task.get_nrow(), "FIXME: Empty task");
        } else {
            // TODO: Here is where we should steal tasks first
            sleep();
        }
        pthread_mutex_unlock(&mutex);
    }

    void kmeans_task_thread::lock_sleep() {
        int rc;
        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");

        (*parent_pending_threads)--;
        set_thread_state(WAIT);

        if (*parent_pending_threads == 0) {
            rc = pthread_cond_signal(parent_cond); // Wake up parent thread
            if (rc) perror("pthread_cond_signal");
        }
        rc = pthread_mutex_unlock(&mutex);
        if (rc) perror("pthread_mutex_unlock");
    }

    // Assumes caller has lock already ... or else ...
    void kmeans_task_thread::sleep() {
        (*parent_pending_threads)--;
        set_thread_state(WAIT);

        if (*parent_pending_threads == 0) {
            int rc = pthread_cond_signal(parent_cond); // Wake up parent thread
            if (rc) perror("pthread_cond_signal");
        }
    }

    void kmeans_task_thread::run() {
        switch(state) {
            case TEST:
                test();
                lock_sleep();
                break;
            case ALLOC_DATA:
                numa_alloc_mem();
                tasks.set_data_ptr(local_data); // We now have real data
                lock_sleep();
                break;
            case KMSPP_INIT:
                kmspp_dist();
                request_task();
                break;
            case EM: /* Super-E-step */
                EM_step();
                request_task();
                break;
            case EXIT:
                fprintf(stderr, "[FATAL]: Thread state is EXIT but running!\n");
                exit(EXIT_FAILURE);
            default:
                fprintf(stderr, "[FATAL]: Unknown thread state\n");
                exit(EXIT_FAILURE);
        }
    }

    void kmeans_task_thread::wait() {
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

    void kmeans_task_thread::wake(thread_state_t state) {
        int rc;
        rc = pthread_mutex_lock(&mutex);
        if (rc) perror("pthread_mutex_lock");
        set_thread_state(state);

        if (state == thread_state_t::EM ||
                state == thread_state_t::KMSPP_INIT) {
            // Threads only sleep if they AND all other threads have no tasks
            tasks.reset(); // NOTE: Only place this is reset
            curr_task = tasks.get_task();
            BOOST_VERIFY(curr_task.get_nrow() <= tasks.get_nrow());

            meta.num_changed = 0; // Always reset at the beginning of an EM-step
            local_clusters->clear();

            //printf("wake: Thd: %u, Task ==> ", get_thd_id()); curr_task.print();
        }

        rc = pthread_mutex_unlock(&mutex);
        if (rc) perror("pthread_mutex_unlock");

        rc = pthread_cond_signal(&cond);
    }

    void* callback(void* arg) {
        kmeans_task_thread* t = static_cast<kmeans_task_thread*>(arg);
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

    void kmeans_task_thread::start(const thread_state_t state=WAIT) {
        //printf("Thread %d started ...\n", thd_id);
        this->state = state;
        int rc = pthread_create(&hw_thd, NULL, callback, this);
        if (rc) {
            fprintf(stderr, "[FATAL]: Thread creation failed with code: %d\n", rc);
            exit(rc);
        }
    }

    const unsigned kmeans_task_thread::
        get_global_data_id(const unsigned row_id) const {
            return row_id + curr_task.get_start_rid();
        }

    void kmeans_task_thread::EM_step() {
        for (unsigned row = 0; row < curr_task.get_nrow(); row++) {
            unsigned true_row_id = get_global_data_id(row);
            unsigned old_clust = cluster_assignments[true_row_id];

            if (prune_init) {
                double dist = std::numeric_limits<double>::max();

                for (unsigned clust_idx = 0;
                        clust_idx < g_clusters->get_nclust(); clust_idx++) {
                    dist = dist_comp_raw(&curr_task.get_data_ptr()[row*ncol],
                            &(g_clusters->get_means()[clust_idx*ncol]), ncol);

                    if (dist < dist_v[true_row_id]) {
                        dist_v[true_row_id] = dist;
                        cluster_assignments[true_row_id] = clust_idx;
                    }
                }

            } else {
                recalculated_v->set(true_row_id, false);
                dist_v[true_row_id] +=
                    g_clusters->get_prev_dist(cluster_assignments[true_row_id]);

                if (dist_v[true_row_id] <=
                        g_clusters->get_s_val(cluster_assignments[true_row_id])) {
                    // Skip all rows
                } else {
                    for (unsigned clust_idx = 0;
                            clust_idx < g_clusters->get_nclust(); clust_idx++) {

                        if (dist_v[true_row_id] <= dm->get(cluster_assignments
                                    [true_row_id], clust_idx)) {
                            // Skip this cluster
                            continue;
                        }

                        if (!recalculated_v->get(true_row_id)) {
                            dist_v[true_row_id] = dist_comp_raw(
                                    &curr_task.get_data_ptr()[row*ncol],
                                    &(g_clusters->get_means()[cluster_assignments
                                        [true_row_id]*ncol]), ncol);
                            recalculated_v->set(true_row_id, true);
                        }

                        if (dist_v[true_row_id] <=
                            dm->get(cluster_assignments[true_row_id], clust_idx)) {
                            // Skip this cluster
                            continue;
                        }

                        // Track 5
                        double jdist = dist_comp_raw(
                                &curr_task.get_data_ptr()[row*ncol],
                                &(g_clusters->get_means()[clust_idx*ncol]), ncol);

                        if (jdist < dist_v[true_row_id]) {
                            dist_v[true_row_id] = jdist;
                            cluster_assignments[true_row_id] = clust_idx;
                        }
                    } // endfor
                }
            }

            BOOST_VERIFY(cluster_assignments[true_row_id] >= 0 &&
                    cluster_assignments[true_row_id] < g_clusters->get_nclust());

            if (prune_init) {
                meta.num_changed++;
            local_clusters->add_member(&(curr_task.get_data_ptr()[row*ncol]),
                    cluster_assignments[true_row_id]);
            } else if (old_clust != cluster_assignments[true_row_id]) {
                meta.num_changed++;
                local_clusters->swap_membership(
                        &(curr_task.get_data_ptr()[row*ncol]),
                        old_clust, cluster_assignments[true_row_id]);
            }
        }
    }



    /** Method for a distance computation vs a single cluster.
      * Used in kmeans++ init
      */
    void kmeans_task_thread::kmspp_dist() {
        BOOST_ASSERT_MSG(false, "Unimplemented!\n");
#if 0
        cuml_dist = 0;
        unsigned clust_idx = meta.clust_idx;
        for (unsigned row = 0; row < nlocal_rows; row++) {
            unsigned true_row_id = get_global_data_id(row);

            double dist = dist_comp_raw(&local_data[row*ncol],
                    &((g_clusters->get_means())[clust_idx*ncol]), ncol);

            if (dist < dist_v[true_row_id]) { // Found a closer cluster than before
                dist_v[true_row_id] = dist;
                cluster_assignments[true_row_id] = clust_idx;
            }
            cuml_dist += dist_v[true_row_id];
        }
#endif
    }

    const void kmeans_task_thread::print_local_data() const {
        print_mat(local_data, (get_data_size()/(sizeof(double)*ncol)), ncol);
    }
}
#endif
