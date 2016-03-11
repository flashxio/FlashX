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
        switch(state) {
            case TEST:
                test();
                break;
            case ALLOC_DATA:
                numa_alloc_mem();
                assert(0);
                break;
            case KMSPP_INIT:
                assert(0);
                break;
            case EM: /*EM steps of kmeans*/
                assert(0);
                break;
            default:
                fprintf(stderr, "Unknown thread state\n");
                exit(EXIT_FAILURE);
        }
    }

    void kmeans_thread::test() {
        std::cout << "Hello from thread: " << get_thd_id() << "\n";
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

            //std::copy(&(full_data[start_offset]),
            //        &(full_data[start_offset+(ncol*nprocrows)]), local_data);
    }

    void kmeans_thread::print_local_data() {
        print_mat(local_data, nprocrows, ncol);
    }
}
