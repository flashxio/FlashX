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

#include <iostream>

#include <boost/assert.hpp>

#include "libgraph-algs/row_cache.h"
#include "matrix/kmeans.h"

unsigned nrow, ncol, cache_size;

void build_state(double* data, const unsigned nrow, const unsigned ncol) {
    unsigned min = 1;
    unsigned max = 5;

    for (unsigned i = 0; i < nrow*ncol; i++)
        data[i] = min + ((double)random() / (double)RAND_MAX * (max - min));
}

void populate_cache(lazy_cache<double>::ptr cache, const double* data) {

    for (unsigned row = 0; row < nrow; row++)
        cache->add(&data[row*ncol], row, ncol);
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        fprintf(stderr, "usage: ./test_lazy_cache nrow ncol cache_size\n");
        exit(EXIT_FAILURE);
    }
    nrow = atoi(argv[1]);
    ncol = atoi(argv[2]);
    cache_size = atoi(argv[3]);

    double* data = new double [nrow*ncol];
    build_state(data, nrow, ncol);
    printf("The data:\n");
    print_mat<double>(data, nrow, ncol);

    lazy_cache<double>::ptr cache;

    for (unsigned i = 0; i <= 30; i++) {
        if (i % 5 == 0) {
            printf("building new cache\n");
            cache = lazy_cache<double>::create(cache_size, ncol);
        }

        for (unsigned row = 0; row < nrow; row++) {
            // Do a get or populate
            if (!cache->get(row)) {
                cache->add(&data[row*ncol], row, ncol);
            } else {
                printf(".");
            }
        }
    }

    delete [] data;
    return EXIT_SUCCESS;
}
