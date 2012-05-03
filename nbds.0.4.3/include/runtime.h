/*
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#ifndef RUNTIME_H
#define RUNTIME_H

#include <pthread.h>
#include "tls.h"

extern DECLARE_THREAD_LOCAL(tid_, int);

int nbd_thread_create (pthread_t *restrict thread, int thread_id, void *(*start_routine)(void *), void *restrict arg);
uint64_t nbd_rand (void);
uint64_t nbd_rand_seed (int i);
int nbd_next_rand (uint64_t *r);

#endif//RUNTIME_H
