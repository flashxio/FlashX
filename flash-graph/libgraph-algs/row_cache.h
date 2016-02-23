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
#include <map>
#include <boost/assert.hpp>
#include <pthread.h>
#include "sem_kmeans_util.h"
#include "log.h"

namespace {
    template <typename T>
        class row_queue {
            private:
                class node {
                    private:
                        node* prev;
                        node* next;
                        unsigned id; // TODO: change to unsigned for experiments
                        std::vector<T> data;

                        node() {
                            this->prev = NULL;
                            this->next = NULL;
                        }

                        node(std::vector<T>& data, node* next=NULL,
                                node* prev=NULL) {
                            this->prev = next;
                            this->next = prev;
                            this->data = data;
                        }

                    public:
                        static node* create() {
                            return new node();
                        }

                        static node* create(std::vector<T>& data,
                                node* next, node* prev) {
                            return new node(data, next, prev);
                        }

                        void set_next(node* next) {
                            this->next = next;
                        }

                        void set_prev(node* prev) {
                            this->prev = prev;
                        }

                        node* get_next() {
                            return this->next;
                        }

                        node* get_prev() {
                            return this->prev;
                        }
                };

                node* sentinel;
            public:
                row_queue(const unsigned size) {
                    sentinel = node::create(0);
                    sentinel->set_next(sentinel);
                    sentinel->set_prev(sentinel);
                }

                // Always add to the front of the queue
                void push_front(std::vector<T>& data) {
                    node* front_next = sentinel->next;
                    // front_prev == sentinel always
                    node* new_node = node::create(data, front_next, sentinel);
                    front_next->set_prev(new_node);
                    sentinel->next = new_node;
                }

                void drop_rear() {
                    BOOST_VERIFY(0);
                    return NULL;
                }

                bool is_empty() {
                    return sentinel->next == sentinel;
                }

                // Iterate & get a the node
                node* find_node(const unsigned id) {
                    BOOST_VERIFY(0);
                    return NULL;
                }

                // Move a node to the front of the list
                void promote_node(const unsigned id) {
                    node* n = find_node(id);
                    // Disconnect n from its neighbors
                    n->get_next()->set_prev(n->get_prev());
                    n->get_prev()->set_next(n->get_next());
                    // Add to front
                    n->set_prev(sentinel);
                    n->set_next(sentinel->get_next());
                    sentinel->get_next()->set_prev(n);
                    // Update sentinel
                    sentinel->set_next(n);
                }

                ~row_queue() {
                    node* current = sentinel->next;
                    while(current != sentinel) {
                        node* tmp = current->get_next();
                        delete current;
                        current = tmp;
                    }
                }
        };

    template <typename T>
        class lazy_cache
        {
            private:
                T* data;
                std::map<unsigned, T*> dmap;
                unsigned max_size; // Max size
                unsigned size; // Actual Size
                unsigned elem_len;
                bool printed;
                pthread_spinlock_t lock;

                lazy_cache(const unsigned max_numel, unsigned elem_len) {
                    size = 0;
                    this->elem_len = elem_len;
                    max_size = max_numel*elem_len;
                    BOOST_VERIFY(max_size > 0);
                    data = new T[max_size];
                    printed = false;
                    pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
                }

                typedef typename std::map<unsigned, T*>::iterator cache_map_iter;

            public:
                typedef std::shared_ptr<lazy_cache> ptr;

                static ptr create(const unsigned max_numel, unsigned elem_len) {
                    return ptr(new lazy_cache(max_numel, elem_len));
                }

                void add(const T* row, const unsigned id, const unsigned len) {
                    if (size < max_size) {
                        //acquire lock
                        pthread_spin_lock(&lock);

                        dmap[id] = &data[size];
                        std::copy(row, row+len, &data[size]);
                        size += len;

                        //drop lock
                        pthread_spin_unlock(&lock);
                    } else {
                        if (!printed) {
                            printf("Printing full cache @ r:%u\n", id);
                            print();
                            printed = true;
                        }
                    }
                    // Unconventional in that we don't care to keep the most recently
                    // found but instead just care that the cache is full
                }

                void print() {
                    for (cache_map_iter it = dmap.begin(); it != dmap.end(); ++it) {
                        printf("#row id:%u ==> ", it->first);
                        printf("[ ");
                        for (unsigned i = 0; i < elem_len; i++) {
                            std::cout << it->second[i] << " ";
                        }
                        printf("]\n");
                    }
                }

                bool is_full() {
                    return size >= max_size;
                }

                T* get(const unsigned id) {
                    cache_map_iter it = dmap.find(id);
                    if (it == dmap.end())
                        return NULL; // Cache miss
                    else
                        return it->second; // Cache hit
                }

                ~lazy_cache() {
                    delete [] data;
                }
        };

    // A paritioned paritioned cache that will be lazily updated
    template <typename T>
    class partition_cache {
        private:
            std::vector<std::vector<T>> data;
            std::vector<std::vector<unsigned>> ids;
            std::map<unsigned, T*> data_map;

            // How many elems have been added since we updated the global count
            std::vector<unsigned> elem_added;
            unsigned numel;
            unsigned max_numel;
            unsigned numel_sync;
            unsigned elem_len;
            unsigned cache_hits;
            pthread_spinlock_t lock;
            static constexpr unsigned MAX_SYNC_ELEM = 100;

            typedef typename std::map<unsigned, T*>::iterator cache_map_iter;

            /**
              * \param numel_sync how many items to insert into a single threads
              *     data elements before synchronizing the cache
              */
            partition_cache(const unsigned nthread, const unsigned elem_len,
                    const unsigned numel_sync, const unsigned max_numel) {
                data.resize(nthread);
                ids.resize(nthread);
                elem_added.resize(nthread); // Default == 0
                this->elem_len = elem_len;
                this->max_numel = max_numel;
                this->numel = 0; // TODO: See if ok non-volatile in practice
                this->cache_hits = 0;

                this->numel_sync = (0 == numel_sync) ? 1 :
                    numel_sync > MAX_SYNC_ELEM ? MAX_SYNC_ELEM:
                    numel_sync;
                BOOST_ASSERT_MSG(numel_sync >= 0, "[ERROR]: param numel_sync <= 0");
                pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);

                BOOST_LOG_TRIVIAL(info) << "\nParams ==> nthread: " << nthread <<
                    ", elem_len: " << this->elem_len << ", max_numel:" << this->max_numel <<
                    ", numel_sync: " << this->numel_sync;
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
                for (unsigned thd = 0; thd < data.size(); thd++) {
                    BOOST_VERIFY(ids[thd].size() == data[thd].size()/elem_len);
                    printf("thd: %u:, size: %lu\n", thd, ids[thd].size());

                    for(unsigned idx = 0; idx < ids[thd].size(); idx++) {
                        printf("row: %u ==> ", ids[thd][idx]);
                        print_arr(&(data[thd][idx*elem_len]), elem_len);
                    }
                }
            }

        public:
            typedef std::shared_ptr<partition_cache> ptr;
            static ptr create(const unsigned nthread, const unsigned elem_len,
                    const unsigned numel_sync, const unsigned max_numel) {
                return ptr(new partition_cache(nthread, elem_len, numel_sync, max_numel));
            }

            // Get ids associated with a thread
            std::vector<unsigned>& get_ids(const unsigned thd) {
                return ids[thd];
            }

            // Each id is added by thread
            // If I cannot add any more ids then this returns false
            bool add_id(const unsigned thd, const unsigned id) {
                if (is_full() || /* Ok because || is short-circuiting */
                        (std::find(ids[thd].begin(), ids[thd].end(), id) != ids[thd].end()))
                    return false;
                 else {
                    ids[thd].push_back(id);
                    elem_added[thd]++;
                    return true;
                }
            }

            // Each vector is added one elem at a time
            void add(const unsigned thd, const T elem, const bool is_end) {
                data[thd].push_back(elem);
                if (is_end) {
                    BOOST_VERIFY(data[thd].size() % elem_len == 0);
                    if (elem_added[thd] == numel_sync) {
                        /*BOOST_LOG_TRIVIAL(info) <<
                            "Syncing numel in cache for thd:" << thd
                            << " with numel: " << numel;*/
                        pthread_spin_lock(&lock); // lock it
                        numel += elem_added[thd];
                        pthread_spin_unlock(&lock); // unlock it
                        elem_added[thd] = 0; // reset
                    }
                }
            }

            const unsigned get_numel(const unsigned thd) const {
                return data[thd].size();
            }

            const bool is_full() const {
                return numel >= max_numel;
            }

            const bool is_empty() {
                return numel == 0;
            }

            const bool index_empty() {
                return data_map.empty();
            }

            void build_index() {
                //print_data_mat();

                BOOST_LOG_TRIVIAL(info) << "Building hash index";
                std::vector<unsigned> zero(0);
                for (unsigned thd = 0; thd < ids.size(); thd++) {
                    T* start_addr = &(data[thd][0]);
                    for (unsigned idx = 0; idx < ids[thd].size(); idx++) {
                        data_map[ids[thd][idx]] = start_addr + (idx * elem_len);
                    }
                    ids[thd].swap(zero); // TODO: Determine efficiency
                }
                //printf("Printing the cache:\n"); print();
            }

            T* get(const unsigned id){
                typename std::map<unsigned, T*>::iterator it = data_map.find(id);
                if (it == data_map.end())
                    return NULL;
                else {
                    // cache_hits++; // TODO:
                    return it->second;
                }
            }

            void print() {
                printf("Printing hashed data with %lu #elems & "
                        "numel = %u\n", data_map.size(), numel);
                for (cache_map_iter it = data_map.begin(); it != data_map.end();
                        ++it) {
                    printf("r:%u ==> ", it->first);
                    print_arr(it->second, elem_len);
                }
            }
    };
}
#endif
