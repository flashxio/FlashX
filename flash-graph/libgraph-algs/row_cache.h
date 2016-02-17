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
#include <memory>
#include <map>


namespace {
    // TODO: Replace queue 
    template <typename T>
        class row_queue {
            private:
                class node {
                    private:
                        node* prev;
                        node* next;
                        size_t id; // TODO: change to unsigned for experiments
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
                row_queue(const size_t size) {
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
                    assert(0);
                    return NULL; 
                }

                bool is_empty() {
                    return sentinel->next == sentinel; 
                }

                // Iterate & get a the node
                node* find_node(const size_t id) {
                    assert(0);
                    return NULL; 
                }

                // Move a node to the front of the list
                void promote_node(const size_t id) {
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
                std::map<size_t, T*> dmap; 
                size_t max_size; // Max size
                size_t size; // Actual Size

                lazy_cache(const size_t max_numel, const size_t elem_len) {
                    size = 0;
                    max_size = max_numel*elem_len;
                    assert(max_size > 0);
                }

            public:
                typedef std::shared_ptr<lazy_cache> ptr;

                static ptr create(const size_t max_numel,
                        const size_t elem_len) {
                    return ptr(new lazy_cache(max_numel, elem_len));
                }

                void add(T* row, const size_t id, const size_t len) {
                    if (size <= max_size) {
                        dmap[id] = &data[size];
                        std::copy(row, row+len , data[size]);
                        size += len;
                    }

                    // Unconventional in that we don't care to keep the most recently
                    // found but instead just care that the cache is full
                }

                bool is_full() {
                    return size > max_size;
                }

                T* get(const size_t id) {
                    typename std::map<size_t, T*>::iterator it = dmap.find(id);
                    if (it == dmap.end())
                        return NULL; // Cache miss
                    else
                        return it->second; // Cache hit
                }
        };
}
#endif
