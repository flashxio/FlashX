/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * lightweight tracing 
 */
#include <stdio.h>
#include "common.h"
#include "rlocal.h"
#include "lwt.h"
#include "mem.h"

#define LWT_BUFFER_SCALE 20
#define LWT_BUFFER_SIZE (1ULL << LWT_BUFFER_SCALE)
#define LWT_BUFFER_MASK (LWT_BUFFER_SIZE - 1)

volatile int halt_ = 0;

typedef struct lwt_record {
    uint64_t timestamp;
    uint64_t format;
    size_t value1;
    size_t value2;
} lwt_record_t;

typedef struct lwt_buffer {
    uint32_t head;
    lwt_record_t x[0];
} lwt_buffer_t;

lwt_buffer_t *lwt_buf_[MAX_NUM_THREADS] = {};
char flag_state_[256] = {};
static const char *flags_ = "";

void lwt_thread_init (int thread_id)
{
    assert(thread_id < MAX_NUM_THREADS);
    if (lwt_buf_[thread_id] == NULL) {
        lwt_buf_[thread_id] = (lwt_buffer_t *)nbd_malloc(sizeof(lwt_buffer_t) + sizeof(lwt_record_t) * LWT_BUFFER_SIZE);
        memset(lwt_buf_[thread_id], 0, sizeof(lwt_buffer_t));
    }
}

void lwt_set_trace_level (const char *flags)
{
    assert(strlen(flags) % 2 == 0); // a well formed <flags> should be an even number of characters long
    flags_ = flags;
    memset(flag_state_, 0, sizeof(flag_state_));
    for (int i = 0; flags[i]; i+=2) {
        flag_state_[(unsigned)flags[i]] = flags[i+1];
    }
}

static void dump_record (FILE *file, int thread_id, lwt_record_t *r, uint64_t offset)
{
    // print the record if its trace category is enabled at a high enough level
    int flag  =  r->format >> 56;
    int level = (r->format >> 48) & 0xFF;
    if (flag_state_[(unsigned)flag] >= level) {
        char s[3] = {flag, level, '\0'};
        fprintf(file, "%09llu %d %s ", ((uint64_t)r->timestamp - offset) >> 5, thread_id, s);
        const char *format = (const char *)(size_t)(r->format & MASK(48)); // strip out the embedded flags
        fprintf(file, format, r->value1, r->value2);
        fprintf(file, "\n");
    }
}

static void dump_buffer (FILE *file, int thread_id, uint64_t offset)
{
    lwt_buffer_t *tb = lwt_buf_[thread_id]; 
    assert(tb);
    if (tb->head > LWT_BUFFER_SIZE) {
        for (int i = tb->head & LWT_BUFFER_MASK; i < LWT_BUFFER_SIZE; ++i) {
            dump_record(file, thread_id, tb->x + i, offset);
        }
    }

    for (int i = 0; i < (tb->head & LWT_BUFFER_MASK); ++i) {
        dump_record(file, thread_id, tb->x + i, offset);
    }
}

void lwt_halt (void) {
    halt_ = 1;
}

void lwt_dump (const char *file_name)
{
    halt_ = 1;
    uint64_t offset = (uint64_t)-1;

    for (int i = 0; i < MAX_NUM_THREADS; ++i) {
        if (lwt_buf_[i] != NULL && lwt_buf_[i]->head != 0) {
            uint64_t x = lwt_buf_[i]->x[0].timestamp;
            if (x < offset) {
                offset = x;
            }
            if (lwt_buf_[i]->head > LWT_BUFFER_SIZE)
            {
                x = lwt_buf_[i]->x[lwt_buf_[i]->head & LWT_BUFFER_MASK].timestamp;
                if (x < offset) {
                    offset = x;
                }
            }
        }
    }

    if (offset != (uint64_t)-1) {
        FILE *file = fopen(file_name, "w");
        assert(file);
        for (int i = 0; i < MAX_NUM_THREADS; ++i) {
            if (lwt_buf_[i] != NULL) {
                dump_buffer(file, i, offset);
            }
        }
        fflush(file);
        fclose(file);
    }
}

void lwt_trace_i (uint64_t format, size_t value1, size_t value2) {
    while (halt_) {}
    LOCALIZE_THREAD_LOCAL(tid_, int);
    lwt_buffer_t *tb = lwt_buf_[tid_];
    if (tb != NULL) {
        unsigned int u, l;
        __asm__ __volatile__("rdtsc" : "=a" (l), "=d" (u)); 
        uint64_t timestamp = ((uint64_t)u << 32) | l; 
        lwt_record_t temp = { timestamp, format, value1, value2 };

        tb->x[tb->head++ & LWT_BUFFER_MASK] = temp;
    }
}
