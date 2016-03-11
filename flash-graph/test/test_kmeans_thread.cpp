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
#include "matrix/kmeans_thread.h"
#include "matrix/kmeans_thread.cpp"
#include "libgraph-algs/sem_kmeans_util.h"

void test_thread_creation(const unsigned NTHREADS, const unsigned nnodes) {
    std::vector<kmeans_thread::ptr> thds;

    // Always: Build state alone
    for (unsigned i = 0; i < NTHREADS; i++) {
        clusters::ptr cl = clusters::create(2,2);
        thds.push_back(kmeans_thread::create
                (i%nnodes, i, 69, 200, 2, cl, NULL, "/dev/null"));
    }

    // Always: Start threads alone
    for (unsigned i = 0; i < NTHREADS; i++)
        thds[i]->start(TEST);

    unsigned sum = 0;
    unsigned verification = 0;
    for (unsigned i = 0; i < thds.size(); i++) {
        thds[i]->join();
        sum += thds[i]->get_val(); // OK to do together
        verification += i;
    }
    BOOST_VERIFY(sum == verification);
    std::cout << "SUCCESS: The sum is: " << sum << "\n";
}

void test_numa_populate_data() {
    constexpr unsigned NTHREADS = 10;
    constexpr unsigned nnodes = 4;
    constexpr unsigned nrow = 50;
    const unsigned nprocrows = nrow/NTHREADS;
    constexpr unsigned ncol = 5;
    const std::string fn = "/mnt/nfs/disa/data/tiny/matrix_r50_c5_rrw.bin";

    std::vector<kmeans_thread::ptr> thds;

    // Always: Build state alone
    for (unsigned i = 0; i < NTHREADS; i++) {
        clusters::ptr cl = clusters::create(2,2);
        thds.push_back(kmeans_thread::create
                (i%nnodes, i, i*nprocrows*ncol, nprocrows, ncol,
                 cl, NULL, fn));
    }

    bin_reader<double> br(fn, nrow, ncol);
    double* data = new double [nrow*ncol];
    printf("Bin read data\n");
    br.read(data);

    // Allocate & move
    std::vector<kmeans_thread::ptr>::iterator it = thds.begin();
    for (; it != thds.end(); ++it)
        (*it)->start(ALLOC_DATA);

    for (it = thds.begin(); it != thds.end(); ++it)
        (*it)->join();

    // Print it back
    for (it = thds.begin(); it != thds.end(); ++it) {
        double *dp = &data[(*it)->get_thd_id()*ncol*nprocrows];
        BOOST_VERIFY(eq_all(dp, (*it)->get_local_data(), nprocrows*ncol));
        printf("Thread %u PASSED numa_mem_alloc()\n", (*it)->get_thd_id());
    }
    delete [] data;
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        fprintf(stderr, "usage: ./test_kmeans_thread nthreads nnodes\n");
        exit(EXIT_FAILURE);
    }

#if KM_TEST
    test_thread_creation(atoi(argv[1]), atoi(argv[2]));
#else
    printf("[FATAL]: Set KM_TEST 1 in kmeans.h\n");
    exit(EXIT_FAILURE);
#endif
    test_numa_populate_data();

    return (EXIT_SUCCESS);
}
