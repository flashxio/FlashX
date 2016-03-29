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

#define TASK_SIZE 4096 // TODO: Vary me
/**
  * These are given to threads to run
  */
class kmeans_task {
    private:
        double* task_data;
        unsigned n_task_rows;
    public:
        kmeans_task(double* task_data, const unsigned n_task_rows) {
            this->task_data = task_data;
            this->n_task_rows = n_task_rows;
        }
        double* get_task_data() { return task_data; }
        unsigned get_n_task_rows () { return n_task_rows; }
};

class task_queue_interface {
    private:
        bool _has_task;

    public:
        virtual kmeans_task* get_task() = 0;
        virtual bool has_task() = 0;
};

// FIXME: Let's assume no NUMA memory allocation first
// NOTE: This should be locked when accessed
class simple_task_queue : public task_queue_interface {
    private:
        unsigned curr_pos; // Where we are in the queue
        double *data; // pointer to the start of the 
        unsigned NROW_PER_TASK;
        unsigned ncol;
        unsigned data_size;

    public:
        typedef std::shared_ptr<simple_task_queue> ptr;

        simple_task_queue(double* data, const unsigned data_size,
                const unsigned ncol) {
            this->data = data;
            this->ncol = ncol;
            this->data_size = data_size;
            NROW_PER_TASK = TASK_SIZE/ncol; // NOTE: integer div
        }

        /**
         * \param data_size The number of rows in the dataset
         **/
        static ptr create(double* data, const unsigned data_size,
                const unsigned ncol) {
            return ptr(new simple_task_queue(data, data_size, ncol));
        }

        // TODO: Check if we really need to alloc this ...
        kmeans_task* get_task() {
            unsigned task_size; 
            if ((curr_pos + NROW_PER_TASK) >= data_size) {
                task_size = data_size - curr_pos; // TODO: Check me
                curr_pos = data_size - 1;
            } else {
                task_size = NROW_PER_TASK;
                curr_pos += NROW_PER_TASK;
            }

            kmeans_task* t = new kmeans_task(&data[curr_pos*ncol], task_size); 
            return t;
        }

        bool has_task() {
            return curr_pos < data_size - 1;
        }
};
#endif
