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
#include <iostream>
#include <boost/assert.hpp>

#include "matrix/thd_safe_bool_vector.h"
#include "libgraph-algs/sem_kmeans_util.h"

using namespace prune;

void build_state(std::vector<short> &verifier,
        thd_safe_bool_vector::ptr data, const unsigned len) {
    for (unsigned i = 0; i < len; i++) {
        short rand_val = rand() % 2;
        if (rand_val == 0) {
            verifier[i] = rand_val;
            data->set(i, false);
        } else if (rand_val == 1) {
            verifier[i] = rand_val;
            data->set(i, true);
        } else {
            BOOST_VERIFY(0);
        }
    }
    printf("Built state successfully ...\n");
}

template <typename T>
void test_correctness(const std::vector<T> &verifier,
        const thd_safe_bool_vector::ptr data) {
    BOOST_VERIFY(verifier.size() == data->size());

    for (unsigned i = 0; i < verifier.size(); i++)
        BOOST_VERIFY((bool)verifier[i] == data->get(i));

    printf("Successfully passed correctness ...\n");
}

void test_init_ctor(const unsigned len) {
    thd_safe_bool_vector::ptr data = thd_safe_bool_vector::create(len, true);
    for (unsigned i = 0; i < len; i++)
        BOOST_VERIFY(data->get(i) == true);

    data = thd_safe_bool_vector::create(len, false);
    for (unsigned i = 0; i < len; i++)
        BOOST_VERIFY(data->get(i) == false);

    printf("Successfully init ctor ...\n");
}

void test_thread_safety(const thd_safe_bool_vector::ptr data,
        const unsigned nthreads) {

    constexpr unsigned NUM_TESTS = 50;
    for (unsigned t = 0; t < NUM_TESTS; t++) {
        // Make the _test vector
        std::vector<bool> _test;
        for (unsigned i = 0; i < data->size(); i++) {
            short rand_val = rand() % 2;
            if (rand_val)
                _test.push_back(true);
            else
                _test.push_back(false);
        }

        // Assign _test values to data
#pragma omp parallel for shared(_test)
        for (unsigned i = 0; i < data->size(); i++) {
            data->set(i, _test[i]);
        }
        test_correctness<bool>(_test, data);
    }
    printf("Successfully passed thread safety ...\n");
}

int main(int argc, char* argv[]) {

    if (argc < 3) {
        fprintf(stderr, "usage: ./test_thd_safe_bool_vector nthreads len\n");
        exit(EXIT_FAILURE);
    }

    unsigned nthreads = atol(argv[1]);
    unsigned len = atol(argv[2]);

    thd_safe_bool_vector::ptr data = thd_safe_bool_vector::create(len);
    std::vector<short> verifier(len);

    build_state(verifier, data, len);
    test_correctness<short>(verifier, data);
    test_init_ctor(len);

    test_thread_safety(data, nthreads);

    return (EXIT_SUCCESS);
}
