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

#ifndef __ROW_CACHE_H__
#define __ROW_CACHE_H__

#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <atomic>
#include <numeric>
#include <boost/assert.hpp>
#include <pthread.h>
#include "log.h"
#include "sem_kmeans_util.h"

namespace {
    static std::vector<unsigned> g_cache_hits;

    template <typename T>
    class partition_cache {
        private:
            std::vector<double> data;
            std::vector<unsigned> ids;
            std::unordered_map<unsigned, T*> data_map;

            std::vector<unsigned> elem_added; // elems added since we updated the global count
            std::vector<unsigned> tot_elem_added; // total elems added
            std::vector<size_t> end_index;

            std::atomic<unsigned> numel;
            unsigned max_numel, numel_sync, elem_len;
            unsigned pt_elem, nthread;
            static constexpr unsigned MAX_SYNC_ELEM = 200;

            typedef typename std::unordered_map<unsigned, T*>::iterator cache_map_iter;

            /**
              * \param numel_sync how many items to insert into a single threads
              *     data elements before synchronizing the cache
              */
            partition_cache(const unsigned nthread, const unsigned elem_len,
                    const unsigned numel_sync, const unsigned max_numel) {
                // Preallocate the mem for each of these
                this->nthread = nthread;
                this->elem_len = elem_len;
                this->pt_elem = ceil(max_numel/(float)nthread);

                data.resize(nthread*pt_elem*elem_len);
                ids.resize(nthread*pt_elem);
                tot_elem_added.resize(nthread);

                for (unsigned thd = 0; thd < nthread; thd++)
                    end_index.push_back(thd*pt_elem*elem_len);

                elem_added.resize(nthread); // Default == 0
                if (g_cache_hits.empty()) {
                    printf("Resizing g_cache_hits!\n");
                    g_cache_hits.resize(nthread);
                }

                this->max_numel = max_numel;
                this->numel = 0; // TODO: See if ok non-volatile in practice

                this->numel_sync = (0 == numel_sync) ? 1 :
                    numel_sync > MAX_SYNC_ELEM ? MAX_SYNC_ELEM:
                    numel_sync;
                BOOST_ASSERT_MSG(numel_sync >= 0, "[ERROR]: param numel_sync <= 0");

                BOOST_LOG_TRIVIAL(info) << "\nParams ==> nthread: " << nthread <<
                    ", elem_len: " << this->elem_len << ", max_numel:"
                    << this->max_numel << ", numel_sync: " << this->numel_sync;
            }

            static void print_arr(T* arr, unsigned len) {
                printf("[ ");
                for (unsigned i = 0; i < len; i++) {
                    std::cout << arr[i] << " ";
                }
                printf("]\n");
            }

            const void print_data_mat() {
                BOOST_LOG_TRIVIAL(info) << "\nEchoing the per-thread data";
                for (unsigned thd = 0; thd < this->elem_added.size(); thd++) {
                    for(unsigned idx = 0; idx < tot_elem_added[thd]; idx++) {
                        if (idx == 0) printf("thd: %u:\n", thd);

                        printf("row: %u ==> ", ids[(thd*pt_elem)+idx]);
                        print_arr(&(data[(thd*pt_elem*elem_len)+idx]), elem_len);
                    }
                }
            }

        public:
            typedef std::shared_ptr<partition_cache> ptr;
            static ptr create(const unsigned nthread, const unsigned elem_len,
                    const unsigned numel_sync, const unsigned max_numel) {
                return ptr(new partition_cache(nthread,
                            elem_len, numel_sync, max_numel));
            }

            // Get ids associated with a thread
            const unsigned* get_ids(const unsigned thd) {
                return &(ids[thd*pt_elem]);
            }

            // Each id is added by thread
            // If I cannot add any more ids then this returns false
            bool add_id(const unsigned thd, const unsigned id) {
                if (is_full(thd)) /* Ok because || is short-circuiting */
                    return false;

                unsigned tmp = thd*pt_elem;
                std::vector<unsigned>::iterator begin = ids.begin()+tmp;
                std::vector<unsigned>::iterator end = begin + tot_elem_added[thd];
                if (std::find(begin, end, id) != end)
                    return false;

                ids[tmp+tot_elem_added[thd]++] = id;
                elem_added[thd]++;
                return true;
            }

            // Each vector is added one elem at a time
            void add(const unsigned thd, const T elem, const bool is_end) {
                data[end_index[thd]++] = elem;
                if (is_end) {
                    if (elem_added[thd] == numel_sync) {
                        numel = numel + elem_added[thd];
                        elem_added[thd] = 0; // reset
                    }
                }
            }

            const bool is_full(const unsigned thd) const {
                return tot_elem_added[thd] == pt_elem;
            }

            const bool index_empty() {
                return data_map.empty();
            }

            void build_index() {
                //BOOST_LOG_TRIVIAL(info) << "Printing data matrix";
                //print_data_mat();
                //BOOST_LOG_TRIVIAL(info) << "Building hash index";

                for (unsigned thd = 0; thd < nthread; thd++) {
                    T* start_addr = &(data[thd*pt_elem*elem_len]);
                    unsigned tmp = thd*pt_elem;
                    for (unsigned idx = 0; idx < tot_elem_added[thd]; idx++) {
                        data_map[ids[tmp+idx]] = start_addr + (idx*elem_len);
                    }
                }
                //printf("Printing the cache:\n"); print();
                //verify();
            }

            T* get(const unsigned id, const unsigned thd){
                cache_map_iter it = data_map.find(id);
                if (it == data_map.end())
                    return NULL;
                else {
                    g_cache_hits[thd]++;
                    return it->second;
                }
            }

            const size_t get_cache_hits() {
                return std::accumulate(g_cache_hits.begin(), g_cache_hits.end(), 0);
            }

            void print() {
                printf("Printing hashed data with %lu #elems & "
                        "numel = %u\n", data_map.size(), (unsigned)numel);
                for (cache_map_iter it = data_map.begin(); it != data_map.end();
                        ++it) {
                    printf("r:%u ==> ", it->first);
                    print_arr(it->second, elem_len);
                }
            }

            bool none_nan(const T* row, const unsigned len, const unsigned id) {
                bool nonenan = true;
                for (unsigned i = 0; i < len; i++) {
                    if (isnan((double)row[i])) {
                        printf("row: %u isnan at index: %u\n", id, i);
                        nonenan = false;
                    }
                }
                return nonenan;
            }

            void verify() {
                for (cache_map_iter it=data_map.begin(); it != data_map.end(); ++it) {
                    none_nan(it->second, elem_len, it->first);
                    /*BOOST_ASSERT_MSG(none_nan(it->second, elem_len, it->first),
                            "[Error] in cache row ID");*/
                }
            }
    };
}
#endif
