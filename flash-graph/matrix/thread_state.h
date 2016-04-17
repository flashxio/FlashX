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

#ifndef __KM_THREAD_STATE_H__
#define __KM_THREAD_STATE_H__

namespace km {
    enum thread_state_t {
        TEST, /*just for testing*/
        ALLOC_DATA, /*moving data for reduces rma*/
        KMSPP_INIT,
        EM, /*EM steps of kmeans*/
        WAIT, /*When the thread is waiting for a new task*/
        EXIT /* Say goodnight */
    };
}

#endif
