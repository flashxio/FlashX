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

#include "kmeans_coordinator.h"

namespace {
    void kmeans_coordinator::numa_alloc_data() {
        //TODO: Pass file handle to threads to read & numa alloc
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->start(ALLOC_DATA);
        join_threads();
    }

    void kmeans_coordinator::join_threads() {
        for (thread_iter it = threads.begin(); it != threads.end(); ++it)
            (*it)->join();
    }
}
