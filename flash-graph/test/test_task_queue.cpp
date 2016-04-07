/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa)
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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdlib.h>
#include <pthread.h>
#include <atomic>

#include "matrix/task_queue.h"
#include "matrix/kmeans.h"

void test_queue_get() {
    constexpr unsigned NTHREADS = 1;
    printf("\n\nRunning test_queue_get with"
            " constexpr NTHREADS = %u...\n", NTHREADS);
    constexpr unsigned nrow = 50;
    constexpr unsigned ncol = 5;
    const std::string fn = "/mnt/nfs/disa/data/tiny/matrix_r50_c5_rrw.bin";

    bin_reader<double> br(fn, nrow, ncol);
    double* data = new double [nrow*ncol];
    printf("Bin read data\n");
    br.read(data);
#if 0
    printf("Raw data print out ...\n");
    print_mat<double>(data, nrow, ncol);
    printf("END Raw data ...\n\n");
#endif
    km::task_queue q(data, 0, nrow, ncol);
    printf("Task queue ==> nrow: %u, ncol: %u\n",
            q.get_nrow(), nrow);

    printf("Printing task queue data ...\n");
    for (unsigned i = 0; i < 2; i++) {
        // Test reset
        q.reset();
        while(q.has_task()) {
            km::task t = q.get_task();
#if 0
            print_mat<double>(t.get_data_ptr(), t.get_nrow(), ncol);
#endif
            BOOST_VERIFY(t.get_start_rid() == (q.get_curr_rid()-MIN_TASK_ROWS));
            BOOST_VERIFY(eq_all(t.get_data_ptr(), &(data[t.get_start_rid()*ncol]),
                        MIN_TASK_ROWS*ncol));
        }
    }

    printf("Task queue test SUCCESSful! ...\n");
    delete [] data;
}

int main(int argc, char* argv[]) {
    test_queue_get();
    return (EXIT_SUCCESS);
}
