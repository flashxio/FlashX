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
#include <omp.h>
#include <iostream>
#include <boost/assert.hpp>
#include "libgraph-algs/row_cache.h"
//#include "matrix/kmeans.h"

unsigned nrow, ncol, cache_size, nthread, numel_sync;

void build_state(double* data, const unsigned nrow, const unsigned ncol) {
    unsigned min = 1;
    unsigned max = 5;

    for (unsigned i = 0; i < nrow*ncol; i++)
        data[i] = min + ((double)random() / (double)RAND_MAX * (max - min));
}

void populate_cache(partition_cache<double>::ptr cache, const double* data) {
    printf("Populating cache\n");

#pragma omp parallel for schedule(static)
    for (unsigned rid = 0; rid < nrow; rid++) {
        unsigned thd_id = omp_get_thread_num();

        if (cache->add_id(thd_id, rid)) {
            std::vector<double> row_data;

            for (unsigned col = 0; col < ncol; col++) {
                row_data.push_back(data[(rid*ncol)+col]);
            }
            cache->add(thd_id, rid, row_data);
        }
    }
    cache->build_index();
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        fprintf(stderr,
                "usage: ./test_partition_cache nrow ncol cache_size nthread\n");
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

    partition_cache<double>::ptr cache;
    unsigned NITER_REBUILD = 5;

    for (unsigned i = 0; i <= 30; i++) {
        if (i % NITER_REBUILD == 0) {
            printf("Creating new cache\n");
            cache = partition_cache<double>::create(nthread, ncol,
                    numel_sync, cache_size);
            // Populate the cache
            populate_cache(cache, data);
        }

#pragma omp parallel for schedule(static)
        for (unsigned rid = 0; rid < nrow; rid++) {
            unsigned thd_id = omp_get_thread_num();
            double* lookup = cache->get(rid, thd_id);
            if (lookup) {
                BOOST_VERIFY(kpmbase::eq_all<double>(&data[rid*ncol], lookup, ncol));
            }
        }

        BOOST_ASSERT_MSG(cache->size()
                != 0, "Cache size of 0 tests nothing!");
        BOOST_ASSERT_MSG(cache->get_cache_hits()
                == cache->size()*((i % NITER_REBUILD)+1),
                "Incorrect # of cache hits!");
    }

    delete [] data;
    printf("Test successful!\n");
    return EXIT_SUCCESS;
}
