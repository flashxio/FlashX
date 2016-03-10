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
#include <vector>
#include "matrix/kmeans_thread.h"
#include "matrix/kmeans_thread.cpp"

void test_thread_creation(const unsigned NTHREADS, const unsigned nnodes) {
    std::vector<kmeans_thread> thds;

    // Always: Build state alone
    for (unsigned i = 0; i < NTHREADS; i++)
        thds.push_back(kmeans_thread(i%nnodes, i, 69, 200));

    // Always: Start threads alone
    for (unsigned i = 0; i < NTHREADS; i++)
        thds[i].start();

    unsigned sum = 0;
    unsigned verification = 0;
    for (unsigned i = 0; i < thds.size(); i++) {
        thds[i].join();
        sum += thds[i].get_val(); // OK to do together
        verification += i;
    }
    BOOST_VERIFY(sum == verification);

    std::cout << "In main the sum is: " << sum << "\n";
    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        fprintf(stderr, "usage: ./test_kmeans_thread nthreads nnodes\n");
        exit(EXIT_FAILURE);
    }

    test_thread_creation(atoi(argv[1]), atoi(argv[2]));
    return (EXIT_SUCCESS);
}
