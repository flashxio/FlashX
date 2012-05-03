/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * fast multi-threaded malloc.
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

#define CHUNK_SCALE 12 // 4k chunks
#define PAGE_SCALE  21 // 2MB pages
#define PAGE_SIZE (1ULL << PAGE_SCALE)

// On both linux and Mac OS X the size of the mmap-able virtual address space is between 2^46 and 2^47. Linux has 
// no problem when you grab the whole thing. Mac OS X apparently does some O(n) thing on the first page fault
// that takes over 2 seconds if you mmap 2^46 bytes. So on Mac OS X we only take 2^38 bytes of virtual space. Which
// is OK though, since you can only buy a Mac with up to 32GB of RAM (as of 2/09).
#ifndef NBD32
#ifdef  __MACOSX__
#define TOTAL_SCALE 38
#else //__MACOSX__
#define TOTAL_SCALE 46
#endif//__MACOSX__
#else// NBD32
#define TOTAL_SCALE 32
#endif//NBD32
#define TOTAL_SIZE (1ULL << TOTAL_SCALE)

#define INVALID_SLAB_CLASS          255
#define METASLAB_CLASS_MAX          2
#define NESTED_4K_SLAB_CLASS_MAX    16
#define NESTED_32K_SLAB_CLASS_MAX   39
#define NESTED_256K_SLAB_CLASS_MAX  63
#define NESTED_SLAB_CLASS_MAX       NESTED_256K_SLAB_CLASS_MAX
#define LARGE_SLAB_CLASS_MAX        93
#define HUGE_SLAB_CLASS_MAX         (sizeof(BlockSize) / sizeof(*BlockSize))
#define SLAB_CLASS_MAX              HUGE_SLAB_CLASS_MAX

#define NESTED_SLAB_CASES NESTED_4K_SLAB_CASES: case NESTED_32K_SLAB_CASES: case NESTED_256K_SLAB_CASES
#define NESTED_4K_SLAB_CASES           METASLAB_CLASS_MAX+1 ... NESTED_4K_SLAB_CLASS_MAX
#define NESTED_32K_SLAB_CASES    NESTED_4K_SLAB_CLASS_MAX+1 ... NESTED_32K_SLAB_CLASS_MAX: case 0
#define NESTED_256K_SLAB_CASES  NESTED_32K_SLAB_CLASS_MAX+1 ... NESTED_SLAB_CLASS_MAX: case 1
#define LARGE_SLAB_CASES            NESTED_SLAB_CLASS_MAX+1 ... LARGE_SLAB_CLASS_MAX: case 2
#define HUGE_SLAB_CASES              LARGE_SLAB_CLASS_MAX+1 ... HUGE_SLAB_CLASS_MAX

#define SLAB_CLASS_SCALE(class) ({                       \
        int _scale = 0;                                  \
        switch (class) {                                 \
        case NESTED_4K_SLAB_CASES:   _scale = 12; break; \
        case NESTED_32K_SLAB_CASES:  _scale = 15; break; \
        case NESTED_256K_SLAB_CASES: _scale = 18; break; \
        case LARGE_SLAB_CASES:       _scale = 21; break; \
        }                                                \
        _scale;                                          \
})

// indexed by class
static const uint32_t BlockSize[] = { 
    // meta slab classes (for the nested slabs)
    1 << 12, 1 << 15, 1 << 18
    
    // nested slab classes (4kB, 32kB, and 256kB)
    8,     16,    24,    32,    40,    48,    56,    64,    72,    80,
    88,    96,    112,   120,   128,   144,   160,   176,   192,   224,   
    256,   288,   320,   352,   384,   416,   448,   480,   512,   576,   
    640,   704,   768,   832,   896,   960,   1024,  1152,  1280,  1408,  
    1536,  1664,  1856,  2048,  2240,  2432,  2688,  2944,  3200,  3520,  
    3840,  4160,  4544,  4928,  5312,  5696,  6144,  6592,  7040,  7488,  
    7936,  

    // large slab classes (full page, 2MB)
    8896,  9984,  11200, 12544, 14016, 15616, 17408, 19328, 21440, 23744, 
    26176, 28800, 31616, 34624, 37760, 41024, 44416, 47936, 51584, 55296, 
    59008, 62784, 66496, 70208, 73856, 77376, 80832, 84160, 87360, 90368, 
    93248, 95936, 98496, 100864, 

    // huge slabs (slabs on huge blocks, 2MB-4MB)
    110912,  121984,  134144,  147520,  162240,  178432,  196224,  215808,  237376,  261056,
    287104,  315776,  347328,  382016,  420160,  462144,  508352,  559168,  615040,  676544,
    744192,  818560,  900416,  990400,  1089408, 1198336, 1318144, 1449920, 1594880, 1754368,
    1929792
};

typedef uint8_t class_t;

typedef struct block {
    struct block *next;
} block_t;

typedef struct slab {
    unsigned valid:1;
    unsigned free_list:15;
    unsigned num_in_use:9;
    unsigned class:6;
} __attribute__((packed)) slab_t;

typedef struct metaslab {
    slab_t slab;
    char * data;
    slab_t slab[1 << (PAGE_SCALE - CHUNK_SCALE)];
    struct {
        struct metaslab *older; 
        struct metaslab *newer; 
    } q[NESTED_SLAB_CLASS_MAX+1];
    uint64_t partial_slab_bitmap2[NESTED_32K_SLAB_CLASS_MAX+1];
    uint8_t  partial_slab_bitmap1[NESTED_SLAB_CLASS_MAX+1]; 
} metaslab_t;

char    *MemBase   = NULL;
char    *MemEnd    = NULL;
char    *PageBreak = NULL;
size_t  *PageMap   = NULL;
block_t *FreePages = NULL;
struct { slab_t *slab; char *slab_base; } ActiveSlab[SLAB_CLASS_MAX + 1] = {};

struct {
    size_t slabs_in_use;
    size_t bytes_requested;
    size_t bytes_allocated;
    size_t total_bytes_allocated;
} ClassStats[METASLAB_CLASS_MAX+1];

struct { 
    slab_t *oldest;
    slab_t *newest; 
} PartialSlabQueue[SLAB_CLASS_MAX+1];

struct {
    slab_t *oldest;
} FreshPartialSlabQueue[SLAB_CLASS_MAX+1];

static block_t *get_block (class_t slab_class);

void mem_init (void) {
    ASSERT(INVALID_SLAB_CLASS > SLAB_CLASS_MAX);

    void *buf = mmap(NULL, TOTAL_SIZE, PROT_NONE, MAP_NORESERVE|MAP_ANON|MAP_PRIVATE, -1, 0);
    if (buf == (void *)-1) {
        perror("mmap");
        exit(-1);
    }
    MemEnd  = buf + TOTAL_SIZE;
    MemBase = (char *)( ((size_t)buf + PAGE_SIZE-1) & ~(PAGE_SIZE-1) ); // align to a page boundry

    size_t page_map_size = sizeof(void *) >> (TOTAL_SCALE - PAGE_SCALE);
    mprotect(MemBase, chunk_map_size, PROT_READ|PROT_WRITE);
    PageBreak = MemBase + chunk_map_size;
    PageMap  = (size_t *)MemBase;
}

static class_t get_slab_class (size_t size) {
    for (int i = METASLAB_CLASS_MAX + 1; i <= SLAB_CLASS_MAX; ++i) {
        if (size <= BlockSize[i])
            return i;
    }
    return INVALID_SLAB_CLASS;
}

static class_t get_meta_class (class_t class) {
    int scale = SLAB_CLASS_SCALE(class);
    if (scale == PAGE_SCALE || scale == 0)
        return INVALID_SLAB_CLASS;
    return (scale - 12) / 3;
}

static void *get_page (void) {
    block_t *p = FreePages;
    if (p == NULL) {
        p = (block_t *)PageBreak;
        PageBreak += PAGE_SIZE;
        return p;
    }
    FreePages = p->next;
    return p;
}

static void free_page (void *p) {
    ASSERT(p < (void *)PageBreak);
    block_t *b = (block_t *)p;
    b->next = FreePages;
    FreePages = b;
}

static void init_slab (void *b, class_t slab_class) {
}

static slab_t *new_large_slab (class_t slab_class) {
    return NULL;
}

static int find_partial_slab(metaslab_t *metaslab, class_t target_class, int target_index) {
    switch (target_class) {
        case NESTED_4K_SLAB_CASSES:
            {
                // search nearby the target first
                int base_index = (target_index & ~0x7);
                for (int i = 0; i < 8; ++i) {
                    if (base_index + i == target_index)
                        continue;
                    if (metaslab->slab[base_index + i].class == target_class)
                        return base_index + i;
                }
                do {
                    metaslab->partial_slab_bitmap2[target_class] &= ~(1ULL << (base_index >> 3));
                    uint64_t bitmap = metaslab->partial_slab_bitmap2[target_class];
                    if (bitmap == 0)
                        return NULL;
                    int n = base_index >> 3;
                    if (bitmap & (0xFF << (n & ~0x7))) {
                        bitmap &= 0xFF << (n & ~0x7); // search nearby the target first
                    }
                    base_index = COUNT_TRAILING_ZEROS(bitmap) << 3;
                    for (int i = 0; i < 8; ++i) {
                        if (metaslab->slab[base_index + i].class == target_class)
                            return base_index + i;
                    }
                } while (1);
            }
        case NESTED_32K_SLAB_CASSES:
            {
                uint64_t bitmap = metaslab->partial_slab_bitmap2[target_class];
                if (bitmap == 0)
                    return NULL;
                int n = target_index >> 3;
                if (bitmap & (0xFF << (n & ~0x7))) {
                    bitmap &= 0xFF << (n & ~0x7); // search nearby the target first
                }
                return COUNT_TRAILING_ZEROS(bitmap) << 3;
            }
        case NESTED_256K_SLAB_CASSES:
            {
                uint8_t bitmap = metaslab->partial_slab_bitmap1[target_class];
                if (bitmap == 0)
                    return NULL;
                return COUNT_TRAILING_ZEROS(bitmap) << 6;
            }
        default:
            ASSERT(FALSE);
            return -1;
    }
}

static void activate_new_slab (class_t slab_class) {
    slab_t *new_slab;
    switch (slab_class) {
        case NESTED_SLAB_CASES:
            int slab_index       = ActiveSlab[slab_class].slab_index;
            metaslab_t *metaslab = ActiveSlab[slab_class].metaslab;

            // First look for a partial slab on the same metaslab as the old active slab.
            new_slab = find_partial_slab(metaslab, slab_class);
            if (new_slab == NULL) {
                // No partial slab on the same metaslab. Remove a metaslab from the front of the queue.
                metaslab_t *metaslab = (metaslab_t *)PartialSlabQueue[slab_class].oldest;
                if (metaslab != NULL) {
                    ASSERT(metaslab->q[slab_class].older == NULL);
                    PartialSlabQueue[slab_class].newest = (slab_t *)metaslab->q[slab_class].newer;
                    metaslab->q[slab_class].newer->q[slab_class].older = NULL;
                    new_slab = find_partial_slab(metaslab, slab_class);
                } else {
                    // Can't find a partial slab; create a new slab.
                    new_slab = (slab_t *)get_block(get_meta_class(slab_class));
                    init_slab(new_slab, slab_class);
                }
            }
            break;

        case LARGE_SLAB_CASES:
        case HUGE_SLAB_CASES:
            // large or huge slab class
            new_slab = PartialSlabQueue[slab_class].oldest;
            if (new_slab == NULL) {
                ASSERT(new_slab->older == NULL);
                PartialSlabQueue[slab_class].newest = new_slab->newer;
                new_slab->newer->older = NULL;
            }
            if (new_slab == NULL) {
                if (IS_HUGE_SLAB_CLASS(slab_class)) {
                    new_slab = new_large_slab(slab_class);
                } else {
                    ASSERT(IS_LARGE_SLAB_CLASS(slab_class));
                    new_slab = (slab_t *)get_page();
                }
                init_slab(new_slab, slab_class);
            }
            break;

        default:
            ASSERT(FALSE);
    }

    ActiveSlab[slab_class] = new_slab;
}

static void *get_block(class_t slab_class) {

    // Look for a free block on the active slab.
    switch (slab_class) {
        case NESTED_SLAB_CASES:
            int slab_index       = ActiveSlab[slab_class].slab_index;
            metaslab_t *metaslab = ActiveSlab[slab_class].metaslab;
            if (metaslab != NULL) {
                slab_t slab = metaslab->slab[slab_index];
                if (slab.free_list) {
                    char *slab_base = metaslab->data + ( ( slab_index - 1 ) << SLAB_CLASS_SCALE(slab_class) );
                    void *b = (void *)( slab_base + ( ( slab.free_list - 1 ) << 3 ) );
                    metaslab->slab[slab_index].free_list = *(uint16_t *)b;
                    return b;
                }
            }
            break;

        case LARGE_SLAB_CASES:
            //TODO
            break;

        case HUGE_SLAB_CASES:
            //TODO
            break;

        default:
            ASSERT(FALSE);
    }

    // Find another slab, activate it, and try again.
    activate_new_slab(slab_class);
    return get_block(slab_class); // recursive tail-call
}

void *nbd_malloc (size_t n) {
    TRACE("m1", "nbd_malloc: size %llu", n, 0);
    if (n == 0)
        return NULL;

    block_t *b = get_block( get_slab_class(n) );

    TRACE("m1", "nbd_malloc: returning block %p", b, 0);
    return b;
}

void nbd_free (void *x) {
    TRACE("m1", "nbd_free: block %p", x, 0);
    ASSERT(x);
    ASSERT(x >= (void *)MemBase && x < (void *)MemEnd);

    block_t *b = (block_t *)x;
    size_t page_index = (size_t)b >> PAGE_SCALE;
    metaslab_t *metaslab = PageMap[page_index];
    ASSERT(metaslab);
    size_t slab_index = ((size_t)b & MASK(PAGE_SCALE)) >> 12;
    slab_t slab = metaslab->slab[slab_index];

    // if <slab> is not valid <b> is on a larger slab.
    if (slab.valid) {
        b->next = slab.free_list;
        // the <offset> of the block is offset by 1 so 0 can represent NULL.
        slab.free_list = ( ((size_t)b & MASK(12)) >> 3 ) + 1;
    } else {
        // <b> is not on a 4kB slab. 
        slab_index &= 0x7; // Try the 32kB slab.
        slab = metaslab->slab[slab_index];
        if (slab.valid) {
            b->next = slab.free_list;
            slab.free_list = ( ((size_t)b & MASK(15)) >> 3 ) + 1;
        } else {
            // <b> is not on a 32kB slab. 
            slab_index &= 0x3F; // <b> must be on the 256kB slab.
            slab = metaslab->slab[slab_index];
            ASSERT(slab.valid);
            b->next = slab.free_list;
            slab.free_list = ( ((size_t)b & MASK(18)) >> 3 ) + 1;
        }
    }
    --slab.num_in_use;
    metaslab->slab[slab_index] = slab;
    if (slab.num_in_use == 0) {
        free_slab(metaslab, slab_index);
    }
}

#else//USE_SYSTEM_MALLOC
#include <stdlib.h>

void mem_init (void) {
    return;
}

void ndb_free (void *x) {
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
