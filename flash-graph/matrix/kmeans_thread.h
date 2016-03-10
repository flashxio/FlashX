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

#include <utility>
#include <boost/format.hpp>
#include "log.h"

namespace {
    class kmeans_thread {

        private:
            pthread_t hw_thd;
            unsigned node_id; // Which NUMA node are you on?
            unsigned thd_id;
            unsigned start_offset; // With respect to the original data
            unsigned nprocrows; // How many rows to process

            friend void* callback(void* arg);
            // Testing
            unsigned val;
        public:
            kmeans_thread(const int node_id, const unsigned thd_id,
                    const size_t start_offset, const unsigned nprocrows) {
                this->node_id = node_id;
                this->thd_id = thd_id;
                this->start_offset = start_offset;
                this->nprocrows = nprocrows;

                printf("Starting thread. Metadata: thd_id: %u, node_id: %u"
                        ", start_offset: %lu, nprocrows: %u\n", thd_id, node_id,
                        start_offset, nprocrows);
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

            // TODO: Figure out why destructor is called so often ...
            ~kmeans_thread() {
                printf("Thread %u being destroyed\n", thd_id);
            }

            const unsigned get_thd_id() {
                return thd_id;
            }

            // Testing
            void set_val(const unsigned val) { this->val = val; }
            const unsigned get_val() const { return val; }

            void run();
            void start();
    };

    void bind2node_id(int node_id) /*TODO: Add static*/
    {
        struct bitmask *bmp = numa_allocate_nodemask();
        numa_bitmask_setbit(bmp, node_id);
        numa_bind(bmp);
        numa_free_nodemask(bmp);
    }
}
#endif
