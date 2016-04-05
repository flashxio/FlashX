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

#include "matrix/kmeans_coordinator.h"

int main (int argc, char* argv []) {

    if (argc < 6) {
        fprintf(stderr, "usage: ./test_kmeans_coordinator nthreads "
                "nrow ncol k datafile [max_iters[tol[initfile]]]\n");
        exit(911);
    }

    size_t nrow = atoi(argv[2]);
    size_t ncol = atoi(argv[3]);
    unsigned k = atoi(argv[4]);
    std::string fn = argv[5];
    unsigned max_iters = -1;
    unsigned nnodes = numa_num_task_nodes();
    unsigned nthreads = atoi(argv[1]);
    double tolerance = -1;

    //std::string init = "kmeanspp";
    //std::string init = "random";
    std::string init = "forgy";

    double* p_centers = NULL;
    if (argc > 6)
        max_iters = atol(argv[6]);
    if (argc > 7)
        tolerance = atof(argv[7]);

    if (argc > 8) {
        p_centers = new double [k*ncol];
        bin_reader<double> br2(argv[7], k, ncol);
        br2.read(p_centers);
        printf("Read centers!\n");
        init = "none";
    }

    kmeans_coordinator::ptr kc = kmeans_coordinator::create(fn,
            nrow, ncol, k, max_iters, nnodes, nthreads, p_centers,
            init, tolerance);
    kc->run_kmeans();
    delete [] p_centers;

    return(EXIT_SUCCESS);
}
