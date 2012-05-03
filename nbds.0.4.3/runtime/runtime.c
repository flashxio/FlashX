/*
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#define _POSIX_C_SOURCE 1 // for rand_r()
#include <stdlib.h>
#include <pthread.h>
#include "common.h"
#include "runtime.h"
#include "rlocal.h"
#include "mem.h"
#include "tls.h"

DECLARE_THREAD_LOCAL(tid_, int);
DECLARE_THREAD_LOCAL(rx_, uint32_t);
DECLARE_THREAD_LOCAL(ry_, uint32_t);
DECLARE_THREAD_LOCAL(rz_, uint32_t);
DECLARE_THREAD_LOCAL(rc_, uint32_t);

typedef struct thread_info {
    int thread_id;
    void *(*start_routine)(void *);
    void *restrict arg;
} thread_info_t;

__attribute__ ((constructor)) void nbd_init (void) {
    INIT_THREAD_LOCAL(r);
    INIT_THREAD_LOCAL(tid_);
    SET_THREAD_LOCAL(tid_, 0);
    mem_init();
    lwt_thread_init(0);
    rcu_thread_init(0);
    srand((uint32_t)rdtsc());
}

static void *worker (void *arg) {
    thread_info_t *ti = (thread_info_t *)arg;
    SET_THREAD_LOCAL(tid_, ti->thread_id);
    LOCALIZE_THREAD_LOCAL(tid_, int);

    SET_THREAD_LOCAL(rx_, rand());
    SET_THREAD_LOCAL(ry_, rand());
    SET_THREAD_LOCAL(rz_, rand());
    SET_THREAD_LOCAL(rc_, rand());

    lwt_thread_init(ti->thread_id);
    rcu_thread_init(ti->thread_id);

    void *ret = ti->start_routine(ti->arg);
    nbd_free(ti);
    return ret;
}

int nbd_thread_create (pthread_t *restrict thread, int thread_id, void *(*start_routine)(void *), void *restrict arg) {
    thread_info_t *ti = (thread_info_t *)nbd_malloc(sizeof(thread_info_t));
    ti->thread_id = thread_id;
    ti->start_routine = start_routine;
    ti->arg = arg;
    return pthread_create(thread, NULL, worker, ti);
}

// George Marsaglia's KISS generator
uint64_t nbd_rand (void) {
    LOCALIZE_THREAD_LOCAL(rx_, unsigned);
    LOCALIZE_THREAD_LOCAL(ry_, unsigned);
    LOCALIZE_THREAD_LOCAL(rz_, unsigned);
    LOCALIZE_THREAD_LOCAL(rc_, unsigned);

    uint32_t rx = 69069 * rx_ + 12345;
    uint32_t ry = ry_;
    uint32_t rz = rz_;
    ry ^= (ry << 13);
    ry ^= (ry >> 17);
    ry ^= (ry <<  5);
    uint64_t t = rz * 698769069LL + rc_;
    uint64_t r = rx + ry + (rz = t);

    SET_THREAD_LOCAL(rx_, rx);
    SET_THREAD_LOCAL(ry_, ry);
    SET_THREAD_LOCAL(rz_, rz);
    SET_THREAD_LOCAL(rc_, t >> 32);

    return r;
}

// Fairly fast random numbers
uint64_t nbd_rand_seed (int i) {
    return rdtsc() + -715159705 + i * 129;
}

int nbd_next_rand (uint64_t *r) {
    *r = (*r * 0x5DEECE66DLL + 0xBLL) & MASK(48);
    return (*r >> 17) & 0x7FFFFFFF;
}
