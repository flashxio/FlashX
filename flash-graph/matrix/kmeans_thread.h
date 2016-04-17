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
#include "thread_state.h"

namespace km {
    class clusters;
}

class kmeans_thread : public km::base_kmeans_thread {
    private:
        std::shared_ptr<km::clusters> g_clusters; // Pointer to global cluster data
        unsigned nprocrows; // How many rows to process

        kmeans_thread(const int node_id, const unsigned thd_id, const unsigned start_rid,
                const unsigned nprocrows, const unsigned ncol,
                std::shared_ptr<km::clusters> g_clusters, unsigned* cluster_assignments,
                const std::string fn);
    public:
        typedef std::shared_ptr<kmeans_thread> ptr;

        static ptr create(const int node_id, const unsigned thd_id,
                const unsigned start_rid, const unsigned nprocrows,
                const unsigned ncol, std::shared_ptr<km::clusters> g_clusters,
                unsigned* cluster_assignments, const std::string fn) {
            return ptr(new kmeans_thread(node_id, thd_id, start_rid,
                        nprocrows, ncol, g_clusters,
                        cluster_assignments, fn));
        }

        void start(const km::thread_state_t state);
        // Allocate and move data using this thread
        void EM_step();
        void kmspp_dist();
        const unsigned get_global_data_id(const unsigned row_id) const;
        void run();
        void wait();
        void sleep();
        void wake(km::thread_state_t state);
        const void print_local_data() const;
};
#endif
