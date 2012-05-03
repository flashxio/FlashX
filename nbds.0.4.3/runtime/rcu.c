/*
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * safe memory reclamation using a simple technique from rcu
 *
 * WARNING: not robust enough for real-world use
 */
#include <string.h>
#include "common.h"
#include "rlocal.h"
#include "lwt.h"
#include "mem.h"
#include "tls.h"
#include "rcu.h"

#define RCU_POST_THRESHOLD 10
#define RCU_QUEUE_SCALE 20

typedef struct fifo {
    uint32_t head;
    uint32_t tail;
    uint32_t scale;
    void *x[0];
} fifo_t;

#define MOD_SCALE(x, b) ((x) & MASK(b))
static uint64_t rcu_[MAX_NUM_THREADS][MAX_NUM_THREADS] = {};
static uint64_t rcu_last_posted_[MAX_NUM_THREADS][MAX_NUM_THREADS] = {};
static fifo_t *pending_[MAX_NUM_THREADS] = {};
static int num_threads_ = 0;

static fifo_t *fifo_alloc(int scale) {
    fifo_t *q = (fifo_t *)nbd_malloc(sizeof(fifo_t) + (1ULL << scale) * sizeof(void *));
    memset(q, 0, sizeof(fifo_t));
    q->scale = scale;
    q->head = 0;
    q->tail = 0;
    return q;
}

void rcu_thread_init (int id) {
    assert(id < MAX_NUM_THREADS);
    if (pending_[id] == NULL) {
        pending_[id] = fifo_alloc(RCU_QUEUE_SCALE);
        (void)SYNC_ADD(&num_threads_, 1);
    }
}

void rcu_update (void) {
    LOCALIZE_THREAD_LOCAL(tid_, int);
    assert(tid_ < num_threads_);
    int next_thread_id = (tid_ + 1) % num_threads_;
    TRACE("r1", "rcu_update: updating thread %llu", next_thread_id, 0);
    int i;
    for (i = 0; i < num_threads_; ++i) {
        if (i == tid_)
            continue;

        // No need to post an update if the value hasn't changed
        if (rcu_[tid_][i] == rcu_last_posted_[tid_][i])
            continue;

        uint64_t x = rcu_[tid_][i];
        rcu_[next_thread_id][i] = rcu_last_posted_[tid_][i] = x;
        TRACE("r2", "rcu_update: posted updated value (%llu) for thread %llu", x, i);
    }

    // free
    fifo_t *q = pending_[tid_];
    while (q->tail != rcu_[tid_][tid_]) {
        uint32_t i = MOD_SCALE(q->tail, q->scale);
        TRACE("r0", "rcu_update: freeing %p from queue at position %llu", q->x[i], q->tail);
        nbd_free(q->x[i]);
        q->tail++;
    }
}

void rcu_defer_free (void *x) {
    assert(x);
    LOCALIZE_THREAD_LOCAL(tid_, int);
    fifo_t *q = pending_[tid_];
    assert(MOD_SCALE(q->head + 1, q->scale) != MOD_SCALE(q->tail, q->scale));
    uint32_t i = MOD_SCALE(q->head, q->scale);
    q->x[i] = x;
    TRACE("r0", "rcu_defer_free: put %p on queue at position %llu", x, q->head);
    q->head++;

    if (pending_[tid_]->head - rcu_last_posted_[tid_][tid_] >= RCU_POST_THRESHOLD) {
        TRACE("r0", "rcu_defer_free: posting %llu", pending_[tid_]->head, 0);
        int next_thread_id = (tid_ + 1) % num_threads_;
        rcu_[next_thread_id][tid_] = rcu_last_posted_[tid_][tid_] = pending_[tid_]->head;
    }
}
