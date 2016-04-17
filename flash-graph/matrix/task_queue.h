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

#ifndef __KMEANS_TASK_QUEUE_H__
#define __KMEANS_TASK_QUEUE_H__

#include <memory>
#include <boost/assert.hpp>
#include "kmeans.h"

#define MIN_TASK_ROWS 8192 // TODO: Change
namespace km {
    template <typename T>
        class data_container {
            private:
                T* data;
                unsigned start_rid; // row id of the first elem
                unsigned nrow;
            public:
                data_container() { }

                data_container(T* data, const unsigned start_rid) {
                    set_data_ptr(data);
                    set_start_rid(start_rid);
                }

                data_container(T* data, const unsigned start_rid,
                        const unsigned nrow) {
                    set_data_ptr(data);
                    set_start_rid(start_rid);
                    set_nrow(nrow);
                }

                void set_data_ptr(T* data) {
                    this->data = data;
                }

                void set_start_rid(const unsigned start_rid) {
                    this->start_rid = start_rid;
                }
                void set_nrow(const unsigned nrow) {
                    this->nrow = nrow;
                }

                T* get_data_ptr() const {
                    return data;
                }

                const unsigned get_start_rid() const {
                    return start_rid;
                }

                const unsigned get_nrow() const {
                    return nrow;
                }

                const void print(const unsigned ncol) const {
                    printf("start_rid: %u, nrow: %u\n",
                            get_start_rid(), get_nrow());
                    print_mat<T>(get_data_ptr(), get_nrow(), ncol);
                }
        };

    // Task sent to a thread to process
    class task : public data_container<double> {
        public:
            task():data_container() { }
            task(double* data, const unsigned start_rid):data_container(data, start_rid) { }
            task(double* data, const unsigned start_rid,
                    const unsigned nrow):data_container(data, start_rid, nrow){ }
    };

    template<typename T>
        class task_queue_interface {
            protected:
                bool _has_task;

            public:
                virtual task* get_task() = 0;
                virtual const bool has_task() const = 0;
                virtual ~task_queue_interface() {};
        };


    // Repr of mem alloc'd generally by a thread
    //  bound to numa node
    class task_queue: public data_container<double>, task_queue_interface<double> {
        private:
            unsigned curr_rid; // Last index (local to the task) processed in the Q
            unsigned ncol;
        public:
            task_queue() {}

            task_queue(double* data, const unsigned start_rid, const unsigned nrow,
                    const unsigned ncol): data_container(data, start_rid, nrow) {
                if (nrow > 0)
                    _has_task = true;

                curr_rid = 0;
                this->ncol = ncol;
            }

            // NOTE: This must be called with a lock taken
            task* get_task () {
                if (!has_task()) {
                    printf("[ERROR]: In get_task() with no task left!!\n");
                    return new task(NULL, -1, 0);
                }
                BOOST_VERIFY(curr_rid < get_nrow());

                // TODO: Make better for when there are only
                //  a few left rows if we give away a task
                task* t = new task(&(get_data_ptr()[curr_rid*ncol]), get_start_rid()+curr_rid);
                if ((curr_rid + MIN_TASK_ROWS) < (get_nrow()-1)) {
                    t->set_nrow(MIN_TASK_ROWS);
                    curr_rid += MIN_TASK_ROWS;
                    return t;
                } else {
                    t->set_nrow(get_nrow()-curr_rid);
                    curr_rid = get_nrow()-1;
                    _has_task = false;
                }
                BOOST_VERIFY(t->get_nrow() > 0);
                return t;
            }

            const bool has_task() const {
                return _has_task;
            }

            const unsigned get_curr_rid() const {
                return curr_rid;
            }

            void set_ncol(const unsigned ncol) {
                this->ncol = ncol;
            }

            const unsigned get_ncol() {
                return ncol;
            }

            void reset() {
                curr_rid = 0;
                if (get_nrow() > 0)
                    _has_task = true;
            }
    };
}
#endif
