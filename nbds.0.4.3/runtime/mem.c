/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * Extreamly fast multi-threaded malloc.
 */
#ifndef USE_SYSTEM_MALLOC
#define _BSD_SOURCE // so we get MAP_ANON on linux
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>
#include "common.h"
#include "rlocal.h"
#include "lwt.h"

#ifndef NBD32
#define MAX_SCALE        36 // allocate blocks up to 64GB (arbitrary, could be bigger)
#define MIN_SCALE         3 // smallest allocated block is 8 bytes
#define MAX_POINTER_BITS 48
#define PAGE_SCALE       21 // 2MB pages
#else
#define MAX_SCALE        31 
#define MIN_SCALE         2 // smallest allocated block is 4 bytes
#define MAX_POINTER_BITS 32
#define PAGE_SCALE       12 // 4KB pages
#endif
#define PAGE_SIZE        (1ULL << PAGE_SCALE)
#define HEADERS_SIZE     (((size_t)1ULL << (MAX_POINTER_BITS - PAGE_SCALE)) * sizeof(header_t))

typedef struct block {
    struct block *next;
} block_t;

// TODO: Break the page header into two parts. The first part is located in the header region. The 
//       second part is located on the page and is only used when there are free items.
typedef struct header {
#ifdef  RECYCLE_PAGES
    struct header *next;
    struct header *prev;
    block_t *free_list; // list of free blocks
    int num_in_use;
#endif//RECYCLE_PAGES
    uint8_t owner; // thread id of owner
    uint8_t scale; // log2 of the block size
} header_t;

#ifdef RECYCLE_PAGES
typedef struct size_class {
    header_t *active_page;
    header_t *oldest_partial;
    header_t *newest_partial;
} size_class_t;
#endif//RECYCLE_PAGES

typedef struct tl {
#ifndef RECYCLE_PAGES
    block_t *free_list[MAX_SCALE+1];
#else
    header_t *free_pages;
    size_class_t size_class[MAX_SCALE+1];
#endif//RECYCLE_PAGES
    block_t *blocks_from[MAX_NUM_THREADS];
    block_t *blocks_to[MAX_NUM_THREADS];
} __attribute__((aligned(CACHE_LINE_SIZE))) tl_t;

static header_t *headers_ = NULL;

static tl_t tl_[MAX_NUM_THREADS] = {};

static inline header_t *get_header (void *r) {
    ASSERT(((size_t)r >> PAGE_SCALE) < HEADERS_SIZE);
    return headers_ + ((size_t)r >> PAGE_SCALE);
}

static void *get_new_region (int block_scale) {
    LOCALIZE_THREAD_LOCAL(tid_, int);
#ifdef RECYCLE_PAGES
    tl_t *tl = &tl_[tid_]; // thread-local data
    if (block_scale <= PAGE_SCALE && tl->free_pages != NULL) {
        void *region = tl->free_pages;
        tl->free_pages = tl->free_pages->next;
        get_header(region)->scale = block_scale;
        return region;
    }
#endif//RECYCLE_PAGES
    size_t region_size = (1ULL << block_scale);
    if (region_size < PAGE_SIZE) {
        region_size = PAGE_SIZE;
    }
    void *region = mmap(NULL, region_size, PROT_READ|PROT_WRITE, MAP_NORESERVE|MAP_ANON|MAP_PRIVATE, -1, 0);
    TRACE("m1", "get_new_region: mmapped new region %p (size %p)", region, region_size);
    if (region == (void *)-1) {
        perror("get_new_region: mmap");
        exit(-1);
    }
    if ((size_t)region & (region_size - 1)) {
        TRACE("m0", "get_new_region: region not aligned", 0, 0);
        munmap(region, region_size);
        region = mmap(NULL, region_size * 2, PROT_READ|PROT_WRITE, MAP_NORESERVE|MAP_ANON|MAP_PRIVATE, -1, 0);
        if (region == (void *)-1) {
            perror("get_new_region: mmap");
            exit(-1);
        }
        TRACE("m0", "get_new_region: mmapped new region %p (size %p)", region, region_size * 2);
        void *aligned = (void *)(((size_t)region + region_size) & ~(region_size - 1));
        size_t extra = (char *)aligned - (char *)region;
        if (extra) {
            munmap(region, extra);
            TRACE("m0", "get_new_region: unmapped extra memory %p (size %p)", region, extra);
        }
        extra = ((char *)region + region_size) - (char *)aligned;
        if (extra) {
            munmap((char *)aligned + region_size, extra);
            TRACE("m0", "get_new_region: unmapped extra memory %p (size %p)", (char *)aligned + region_size, extra);
        }
        region = aligned;
    }
    assert(region);

    header_t *h = get_header(region);
    TRACE("m1", "get_new_region: header %p (%p)", h, h - headers_);
    assert(h->scale == 0);
    h->scale = block_scale;
    h->owner = tid_;

    return region;
}

void mem_init (void) {
    assert(headers_ == NULL);
    // Allocate space for the page headers. This could be a big chunk of memory on 64 bit systems,
    // but it just takes up virtual address space. Physical space used by the headers is still 
    // proportional to the amount of memory the user mallocs.
    headers_ = mmap(NULL, HEADERS_SIZE, PROT_READ|PROT_WRITE, MAP_ANON|MAP_PRIVATE, -1, 0);
    TRACE("m1", "mem_init: header page %p", headers_, 0);

    // initialize spsc queues
    for (int i = 0; i < MAX_NUM_THREADS; ++i) {
        for (int j = 0; j < MAX_NUM_THREADS; ++j) {
            if (i != j) {
                tl_[i].blocks_to[j] = (block_t *)&(tl_[j].blocks_from[i]);
            }
        }
    }
}

void nbd_free (void *x) {
    TRACE("m1", "nbd_free: block %p page %p", x, (size_t)x & ~MASK(PAGE_SCALE));
    ASSERT(x);
    LOCALIZE_THREAD_LOCAL(tid_, int);
    block_t  *b = (block_t *)x;
    header_t *h = get_header(x);
    int b_scale = h->scale;
    TRACE("m1", "nbd_free: header %p scale %llu", h, b_scale);
    ASSERT(b_scale && b_scale <= MAX_SCALE);
#ifdef RECYCLE_PAGES
    if (b_scale > PAGE_SCALE) {
        int rc = munmap(x, 1ULL << b_scale);
        ASSERT(rc == 0);
        rc = rc;
    }
#endif
#ifndef NDEBUG
    memset(b, 0xcd, (1ULL << b_scale)); // bear trap
#endif
    tl_t *tl = &tl_[tid_]; // thread-local data
    if (h->owner == tid_) {
        TRACE("m1", "nbd_free: private block, old free list head %p", tl->free_list[b_scale], 0);

#ifndef RECYCLE_PAGES
        b->next = tl->free_list[b_scale];
        tl->free_list[b_scale] = b;
#else //RECYCLE_PAGES
        b->next = h->free_list;
        h->free_list = b;
        h->num_in_use--;
        size_class_t *sc = &tl->size_class[b_scale];
        if (sc->active_page != h) {
            if (h->num_in_use == 0) {
                // remove <h> from the partial-page list
                if (h->next != NULL) { h->next->prev = h->prev; }
                if (h->prev != NULL) { h->prev->next = h->next; }
                // put <h> on the free-page list
                h->next = tl->free_pages;
                tl->free_pages = h;
            } else {
                // move <h> to the top of the partial-page list
                if (h->next != NULL) {
                    h->next->prev = h->prev;
                    if (h->prev != NULL) { h->prev->next = h->next; }
                    h->prev = sc->newest_partial;
                    h->next = NULL;
                    sc->newest_partial = h;
                }
            }
        }
#endif//RECYCLE_PAGES
    } else {
        // push <b> onto it's owner's queue
        int b_owner = h->owner;
        TRACE("m1", "nbd_free: owner %llu", b_owner, 0);
        
        // The assignment statements are volatile to prevent the compiler from reordering them.
        VOLATILE_DEREF(b).next = NULL; 
        VOLATILE_DEREF(tl->blocks_to[b_owner]).next = b;

        tl->blocks_to[b_owner] = b;
    }
}

static inline void process_incoming_blocks (tl_t *tl) {
    for (int p = 0; p < MAX_NUM_THREADS; ++p) {
        block_t *b = tl->blocks_from[p];
        if (EXPECT_FALSE(b == NULL)) continue; // the queue is completely empty

        // Leave the last block on the queue. Removing the last block on the queue would create a
        // race with the producer thread putting a new block on the queue.
        for (block_t *next = b->next; next != NULL; b = next, next = b->next) {
            // push <b> onto the appropriate free list
#ifndef RECYCLE_PAGES
            int b_scale = get_header(b)->scale;
            b->next = tl->free_list[b_scale];
            tl->free_list[b_scale] = b;
#else //RECYCLE_PAGES
            header_t *h = get_header(b);
            b->next = h->free_list;
            h->free_list = b;
#endif//RECYCLE_PAGES
        }
        tl->blocks_from[p] = b;
    }
}

static inline block_t *pop_free_list (tl_t *tl, int scale) {
#ifndef RECYCLE_PAGES
    block_t **free_list = &tl->free_list[scale];
#else //RECYCLE_PAGES
    size_class_t *sc = &tl->size_class[scale];
    if (EXPECT_FALSE(sc->active_page == NULL))
        return NULL;
    block_t **free_list = &sc->active_page->free_list;
#endif//RECYCLE_PAGES
    block_t *b = *free_list;
    if (EXPECT_FALSE(b == NULL))
        return NULL;
    ASSERT(get_header(b)->scale == scale);
    *free_list = b->next;
    return b;
}

// Allocate a block of memory at least size <n>. Blocks are binned in powers-of-two. Round up <n> to
// the nearest power of two. 
//
// First check the current thread's free list for an available block. If there are no blocks on the
// free list, pull items off of the current thread's incoming block queues and push them onto the 
// free list. If we didn't get an appropriate size block off of the block queues then allocate a new
// page, break it up into blocks and push them onto the free list. 
void *nbd_malloc (size_t n) {
    // the scale is the log base 2 of <n>, rounded up
    int b_scale = (sizeof(void *) * __CHAR_BIT__) - __builtin_clzl((n) - 1);
    TRACE("m1", "nbd_malloc: size %llu (scale %llu)", n, b_scale);

    if (EXPECT_FALSE(b_scale < MIN_SCALE)) { b_scale = MIN_SCALE; }
    if (EXPECT_FALSE(b_scale > MAX_SCALE)) { return NULL; }

    LOCALIZE_THREAD_LOCAL(tid_, int);
    tl_t *tl = &tl_[tid_]; // thread-local data

    block_t *b = pop_free_list(tl, b_scale);
    if (b != NULL) {
        TRACE("m1", "nbd_malloc: returning block %p", b, 0);
        return b;
    assert(b);
    }

    // The free list is empty so process blocks freed from other threads and then check again.
    process_incoming_blocks(tl);
    b = pop_free_list(tl, b_scale);
    if (b != NULL) {
        TRACE("m1", "nbd_malloc: returning block %p", b, 0);
        return b;
    assert(b);
    }

#ifdef  RECYCLE_PAGES
    // The current active page is completely allocated. Make the oldest partially allocated page 
    // the new active page.
    size_class_t *sc = &tl->size_class[b_scale];
    if (sc->oldest_partial != NULL) {
        sc->active_page = sc->oldest_partial;
        sc->oldest_partial = sc->oldest_partial->next;
        sc->oldest_partial->prev = NULL;
        b = pop_free_list(tl, b_scale);
        ASSERT(b != NULL);
        TRACE("m1", "nbd_malloc: returning block %p", b, 0);
        return b;
    assert(b);
    }
    // There are no partially allocated pages so get a new page.

#endif//RECYCLE_PAGES

    // Get a new page.
    char *page = get_new_region(b_scale);
    b = (block_t *)page; // grab the first block on the page

    // Break up the remainder of the page into blocks and put them on the free list. Start at the
    // end of the page so that the free list ends up in increasing order, for ease of debugging.
    if (b_scale < PAGE_SCALE) {
        size_t block_size = (1ULL << b_scale);
        block_t *head = NULL;
        for (int offset = PAGE_SIZE - block_size; offset > 0; offset -= block_size) {
            block_t *x = (block_t *)(page + offset);
            x->next = head; head = x;
        }
#ifndef RECYCLE_PAGES
        tl->free_list[b_scale] = head;
#else //RECYCLE_PAGES
        sc->active_page = get_header(page);
        sc->active_page->free_list = head;
#endif//RECYCLE_PAGES
    }

    TRACE("m1", "nbd_malloc: returning block %p from new region %p", b, (size_t)b & ~MASK(PAGE_SCALE));
    assert(b);
    return b;
}
#else//USE_SYSTEM_MALLOC
#include <stdlib.h>
#include "common.h"
#include "rlocal.h"
#include "lwt.h"

void mem_init (void) {
    return;
}

void nbd_free (void *x) {
    TRACE("m1", "nbd_free: %p", x, 0);
#ifndef NDEBUG
    memset(x, 0xcd, sizeof(void *)); // bear trap
#endif//NDEBUG
    free(x);
    return;
}

void *nbd_malloc (size_t n) {
    TRACE("m1", "nbd_malloc: request size %llu", n, 0);
    void *x = malloc(n);
    TRACE("m1", "nbd_malloc: returning %p", x, 0);
    return x;
}
#endif//USE_SYSTEM_MALLOC
