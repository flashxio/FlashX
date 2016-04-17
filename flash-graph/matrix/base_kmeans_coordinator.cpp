#include <boost/assert.hpp>

#include "common.h"
#include "log.h"
#include "base_kmeans_coordinator.h"

using namespace km;

base_kmeans_coordinator::base_kmeans_coordinator(const std::string fn, const size_t nrow,
        const size_t ncol, const unsigned k, const unsigned max_iters,
        const unsigned nnodes, const unsigned nthreads,
        const double* centers, const init_type_t it,
        const double tolerance, const dist_type_t dt) {

    this->fn = fn;
    this->nrow = nrow;
    this->ncol = ncol;
    this->k = k;
    BOOST_ASSERT_MSG(k >= 1, "[FATAL]: 'k' must be >= 1");
    this->max_iters = max_iters;
    this->nnodes = nnodes;
    this->nthreads = nthreads;
    if (nthreads >  (unsigned)get_num_omp_threads()) {
        BOOST_LOG_TRIVIAL(warning) << "[WARNING]: Exceeded system"
            " #virtual cores of: " << get_num_omp_threads();
    }
    this->_init_t = it;
    this->tolerance = tolerance;
    this->_dist_t = dt;
    num_changed = 0;
    pending_threads = 0;

    BOOST_VERIFY(cluster_assignments = new unsigned [nrow]);
    BOOST_VERIFY(cluster_assignment_counts = new unsigned [k]);

    std::fill(&cluster_assignments[0],
            (&cluster_assignments[0])+nrow, -1);
    std::fill(&cluster_assignment_counts[0],
            (&cluster_assignment_counts[0])+k, 0);

    // Threading
    pending_threads = 0; // NOTE: This must be initialized
    pthread_mutexattr_init(&mutex_attr);
    pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
    pthread_mutex_init(&mutex, &mutex_attr);
    pthread_cond_init(&cond, NULL);
}

void base_kmeans_coordinator::wait4complete() {
    //printf("Coordinator entering wait4complete ..\n");
    pthread_mutex_lock(&mutex);
    while (pending_threads != 0) {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
    //printf("Coordinator exiting wait4complete!!\n");
}
