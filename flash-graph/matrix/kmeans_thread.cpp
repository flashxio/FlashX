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

#include <iostream>
#include <boost/assert.hpp>

#include "kmeans_thread.h"

namespace {
    void kmeans_thread::run() {
        std::cout << "Hello from thread: " <<
            get_thd_id() << "\n";
        set_val(get_thd_id());
    }

    void* callback(void* arg) {
        kmeans_thread* t = static_cast<kmeans_thread*>(arg);
        t->run();
        pthread_exit(NULL);
    }

    void kmeans_thread::start() {
        //BOOST_VERIFY(hw_thd == 0);
        int rc = pthread_create(&hw_thd, NULL, callback, this);
        printf("Thread completed run method\n");
        if (rc) {
            fprintf(stderr, "[FATAL]: Thread creation failed with code: %d\n", rc);
                exit(rc);
        }
    }

#if 0
    // Move data equally to all nodes
    std::vector<double*> numa_data_move(double* data, const unsigned nnodes,
            const unsigned nrow, const unsigned ncol) {
        nrow_tuple rt = get_rows_per_node(nrow, ncol, nnodes);

        std::vector<double*> ret_addrs;
        for (unsigned node_id = 0; node_id < nnodes; node_id++) {
            size_t start_offset = node_id*rt.first*ncol;
            if (node_id + 1 < nnodes) {
                size_t blob_size = sizeof(double)*rt.first*ncol;
                ret_addrs.push_back(
                        (double*)numa_alloc_onnode(blob_size, node_id));
                std::move(&(data[start_offset]),
                        &(data[start_offset+(ncol*rt.first)]), ret_addrs.back());
            } else {
                size_t blob_size = sizeof(double)*rt.second*ncol;
                ret_addrs.push_back(
                        (double*)numa_alloc_onnode(blob_size, node_id));
                std::move(&(data[start_offset]),
                        &(data[start_offset+(ncol*rt.second)]), ret_addrs.back());
            }
        }
        return ret_addrs;
    }

    // <First -> rows_per_node, Last -> rows_last_node>
    void numa_destroy_mem(std::vector<double*>& data_addrs, const unsigned nnodes,
            const unsigned nrow, const unsigned ncol) {
        nrow_tuple rt = compute_rows_per_node(nrow, ncol, nnodes);
        for (unsigned node_id = 0; node_id < nnodes; node_id++) {
            size_t blob_size;
            if (node_id + 1 < nnodes) {
                blob_size = sizeof(double)*rt.first*ncol;
            } else {
                blob_size = sizeof(double)*rt.second*ncol;
            }
            BOOST_VERIFY(blob_size > 0);
            //printf("numa_free: %lu bytes\n", blob_size);
            numa_free(data_addrs[node_id], blob_size);
        }
    }

#endif
}
