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

#include <math.h>
#include <iostream>
#include <boost/assert.hpp>
#include "libgraph-algs/row_cache.h"
#include "matrix/kmeans.h"

unsigned nrow, ncol, cache_size, nthread, numel_sync;

void build_state(double* data, const unsigned nrow, const unsigned ncol) {
    unsigned min = 1;
    unsigned max = 5;

    for (unsigned i = 0; i < nrow*ncol; i++)
        data[i] = min + ((double)random() / (double)RAND_MAX * (max - min));
}

void populate_cache(partition_cache<double>::ptr cache, const double* data) {
    cache = partition_cache<double>::create(nthread, ncol, 
            cache_size/(nthread*2), cache_size);
    for (unsigned row = 0; row < nrow; row++) {
        cache->add_id(row % nthread, row); // Fake threading
        for (unsigned col = 0; col < ncol; col++) {
            cache->add(row % nthread, data[(row*ncol)+col], col == ncol); // Fake threading
        }
    }
    cache->build_index();
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        fprintf(stderr, "usage: ./test_partition_cache nrow ncol cache_size nthread\n");
        exit(EXIT_FAILURE);
    }
    nrow = atoi(argv[1]);
    ncol = atoi(argv[2]);
    cache_size = atoi(argv[3]);
    nthread = atoi(argv[4]);
    numel_sync = cache_size/(nthread);

    printf("Input params: nrow:%u, ncol:%u, cache_size:%u,"
            "nthread:%u, numel_sync:%u\n\n", nrow, ncol,
            cache_size, nthread, numel_sync);

    double* data = new double [nrow*ncol];
    build_state(data, nrow, ncol);
    printf("The data:\n");
    print_mat<double>(data, nrow, ncol);

    partition_cache<double>::ptr cache;

    for (unsigned i = 0; i <= 30; i++) {
        if (i % 5 == 0) {
            printf("building new cache\n");
            cache = partition_cache<double>::create(nthread, ncol, 
                    numel_sync, cache_size);
        } else if (cache->index_empty()){
            cache->build_index();
        }

        for (unsigned row = 0; row < nrow; row++) {
            // Do a get or populate
            if (!cache->get(row)) {
                if (cache->add_id(row % nthread, row)) { // Fake threading
                    for (unsigned col = 0; col < ncol; col++) {
                        // Fake threading
                        cache->add(row % nthread, data[(row*ncol)+col], col+1 == ncol);
                    }
                }
            } else {
                printf("Cache hit for row:%u\n", row); // Cache hit
            }
        }

        printf("Printing the cache:\n");
        cache->print();
        printf("End print!\n");

    }

    delete [] data;
    return EXIT_SUCCESS;
}
