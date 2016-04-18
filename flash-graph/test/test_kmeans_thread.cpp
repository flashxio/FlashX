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

#include "matrix/kmeans_thread.h"
#include "libgraph-algs/clusters.h"
#include "libgraph-algs/sem_kmeans_util.h"

using namespace km;

static std::atomic<unsigned> pending_threads;
//static unsigned pending_threads;
static pthread_mutex_t mutex;
static pthread_cond_t cond;
static pthread_mutexattr_t mutex_attr;

static void wait4complete() {
    //printf("\nParent entering wait4complete ..\n");
    pthread_mutex_lock(&mutex);
    while (pending_threads != 0) {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    //printf("Exiting wait4complete!!\n\n");
}

static void wake4run(std::vector<kmeans_thread::ptr>& threads,
        const unsigned nthreads, const thread_state_t state) {
    pending_threads = nthreads;
    for (unsigned thd_id = 0; thd_id < threads.size(); thd_id++) {
        threads[thd_id]->wake(state);
    }
}

static void test_thread_creation(const unsigned NTHREADS, const unsigned nnodes) {
    std::vector<kmeans_thread::ptr> threads;

    // Always: Build state alone
    for (unsigned i = 0; i < NTHREADS; i++) {
        clusters::ptr cl = clusters::create(2,2);
        threads.push_back(kmeans_thread::create
                (i%nnodes, i, 69, 200, 1, cl, NULL, "/dev/null"));
        threads[i]->set_parent_cond(&cond);
        threads[i]->set_parent_pending_threads(&pending_threads);
        threads[i]->start(WAIT); // Thread puts itself to sleep
    }

    for (unsigned i = 0; i < 2048; i++) {
        wake4run(threads, NTHREADS, TEST);
        wait4complete();
    }

    wake4run(threads, NTHREADS, EXIT);
    std::cout << "SUCCESS: for creation & join\n";
}

void test_numa_populate_data() {
    constexpr unsigned NTHREADS = 10;
    printf("\n\nRunning test_numa_populate_data with"
            " constexpr NTHREADS = %u...\n", NTHREADS);
    constexpr unsigned nnodes = 4;
    constexpr unsigned nrow = 50;
    const unsigned nprocrows = nrow/NTHREADS;
    constexpr unsigned ncol = 5;
    const std::string fn = "/mnt/nfs/disa/data/tiny/matrix_r50_c5_rrw.bin";

    std::vector<kmeans_thread::ptr> threads;

    // Always: Build state alone
    for (unsigned i = 0; i < NTHREADS; i++) {
        clusters::ptr cl = clusters::create(2,2);
        threads.push_back(kmeans_thread::create
                (i%nnodes, i, i*nprocrows, nprocrows, ncol,
                 cl, NULL, fn));
        threads[i]->set_parent_cond(&cond);
        threads[i]->set_parent_pending_threads(&pending_threads);
        threads[i]->start(WAIT); // Thread puts itself to sleep
    }

    bin_reader<double> br(fn, nrow, ncol);
    double* data = new double [nrow*ncol];
    printf("Bin read data\n");
    br.read(data);

    wake4run(threads, NTHREADS, ALLOC_DATA);
    wait4complete();

    std::vector<kmeans_thread::ptr>::iterator it = threads.begin();
    // Print it back
    for (it = threads.begin(); it != threads.end(); ++it) {
        double *dp = &data[(*it)->get_thd_id()*ncol*nprocrows];
        BOOST_VERIFY(eq_all(dp, (*it)->get_local_data(), nprocrows*ncol));
        printf("Thread %u PASSED numa_mem_alloc()\n", (*it)->get_thd_id());
    }

    wake4run(threads, NTHREADS, EXIT);
    delete [] data;
    printf("SUCCESS test_numa_populate_data ..\n");
}

int main(int argc, char* argv[]) {
    pending_threads = 0; // NOTE: This must be initialized
    pthread_mutexattr_init(&mutex_attr);
    pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
    pthread_mutex_init(&mutex, &mutex_attr);
    pthread_cond_init(&cond, NULL);

    if (argc < 2) {
        fprintf(stderr, "usage: ./test_kmeans_thread nthreads nnodes\n");
        exit(EXIT_FAILURE);
    }

    test_thread_creation(atoi(argv[1]), atoi(argv[2]));
    test_numa_populate_data();

    pthread_cond_destroy(&cond);
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mutex_attr);

    return (EXIT_SUCCESS);
}
