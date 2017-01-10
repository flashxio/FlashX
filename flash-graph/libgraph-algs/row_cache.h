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
    template <typename T>
    class partition_cache {
        private:
            // Each thread has it's own
            typedef typename std::vector<std::unordered_map<unsigned,
                    std::vector<T> > > local_caches;
            typedef typename std::unordered_map<unsigned,
                    std::vector<T> >::iterator cache_iter;

             local_caches pt_data;

            std::unordered_map<unsigned, T*> data_map;

            // elems added since we updated the global count
            std::vector<unsigned> elem_added;
            std::vector<unsigned> tot_elem_added; // total elems added

            std::atomic<unsigned> numel;
            unsigned max_numel, numel_sync, elem_len;
            unsigned pt_elem, nthread;
            static constexpr unsigned MAX_SYNC_ELEM = 200;

            typedef typename std::unordered_map<unsigned, T*>::iterator
                cache_map_iter;

            std::vector<size_t> cache_hits; // Per thread so no locks

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

                for (unsigned thd = 0; thd < nthread; thd++) {
                    // All on the stack!
                    std::unordered_map<unsigned, std::vector<T> > cache;
                    pt_data.push_back(cache);
                }

                tot_elem_added.assign(nthread, 0);

                elem_added.assign(nthread, 0); // Default == 0
                cache_hits.assign(nthread, 0);

                this->max_numel = max_numel;
                this->numel = 0;
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
                    printf("thd: %u:\n", thd);
                    cache_iter it = pt_data[thd].begin();

                    for (; it != pt_data[thd].end(); ++it) {
                        printf("row: %u ==> ", it->first);
                        print_vector<T>(it->second);
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

            // Each id is added by thread
            // If I cannot add any more vids OR I already have this vid,
                // then return false
            bool add_id(const unsigned thd, const unsigned id) {
                if (is_full(thd)) /* Ok because || is short-circuiting */
                    return false;

                cache_iter it = pt_data[thd].find(id);
                if (it != pt_data[thd].end()) {
                    // vid already exists in cache
                    return false;
                }

                elem_added[thd]++;
                tot_elem_added[thd]++;
                return true;
            }

            void add(const unsigned thd, const unsigned id,
                    const std::vector<T>& row) {
                pt_data[thd][id] = row; // TODO: Verify - Calls copy constructor

                if (elem_added[thd] == numel_sync) {
                    numel = numel + elem_added[thd];
                    elem_added[thd] = 0; // reset
                }

                if (pt_data[thd].size() > pt_elem)
                    BOOST_LOG_TRIVIAL(info) << "thd: " << thd << " has " <<
                        pt_data[thd].size() << " items";
                BOOST_VERIFY(pt_data[thd].size() < (pt_elem + 1));
            }

            const bool is_full(const unsigned thd) const {
                return tot_elem_added[thd] == pt_elem;
            }

            const size_t size() {
                return std::accumulate(tot_elem_added.begin(),
                        tot_elem_added.end(), 0);
            }

            const bool index_empty() {
                return data_map.empty();
            }

            /** \brief Build an index of vids => addresses.
              *     The hope is to reduce the cost of merging per thread maps
              *     with data into a single queriably item. This reduces the
              *     cost of copying actual data.
              */
            void build_index() {
#ifdef KM_TEST
                BOOST_LOG_TRIVIAL(info) << "Printing data matrix";
                print_data_mat();
                BOOST_LOG_TRIVIAL(info) << "Building hash index";
#endif
                for (unsigned thd = 0; thd < nthread; thd++) {
                    for (cache_iter it = pt_data[thd].begin();
                            it != pt_data[thd].end(); ++it) {
                        data_map[it->first] = &(it->second[0]);
                    }
                }
#ifdef KM_TEST
                print();
                verify();
#endif
            }

            T* get(const unsigned id, const unsigned thd){
                cache_map_iter it = data_map.find(id);
                if (it == data_map.end())
                    return NULL;
                else {
                    cache_hits[thd]++;
                    return it->second;
                }
            }

            const size_t get_cache_hits() {
                return std::accumulate(cache_hits.begin(),
                        cache_hits.end(), 0);
            }

            void print() {
                printf("Printing cache data with %lu #elems & "
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
                for (cache_map_iter it = data_map.begin();
                        it != data_map.end(); ++it) {
                    none_nan(it->second, elem_len, it->first);
                    BOOST_ASSERT_MSG(none_nan(it->second, elem_len, it->first),
                            "[Error] in cache row ID");
                }
            }
    };
}
#endif
