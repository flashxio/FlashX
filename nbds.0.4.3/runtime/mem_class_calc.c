#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int     uint32_t;

#define CACHE_LINE_SCALE 6

// Return the expected fraction of bytes wasted per slab.
//
// The internal fragmentation due to using size classes is biased by including the space required,
// for a pointer to each block.
double calc_frag(int slab_size, int block_size, int delta)
{
    double quant = (double)delta / 2 / block_size;
    assert(quant >= 0.0);
    int blocks_per_slab = (int)(slab_size / block_size);
 
    // internal fragmentation that comes from tiling non-power-of-2 sized blocks in slabs
    int extra_space = slab_size - blocks_per_slab * block_size; 
    assert(extra_space < block_size);

    // number of different cache line colors needed to evenly distribute cache line accesses
    int num_colors = block_size >> CACHE_LINE_SCALE;
    if (num_colors <= 1)
        return (double)extra_space/slab_size + quant;

    int num_overflow = num_colors - 1 - (extra_space >> CACHE_LINE_SCALE);
    if (num_overflow <= 0)
        return (double)extra_space/slab_size + quant;

    double coloring = (double)num_overflow * block_size / num_colors;
    return ((double)extra_space + coloring)/slab_size + quant;
}

// size classes for various alignments, max 6% expected internal fragmentation

// 2B-128B blocks, 4k slab
static uint8_t  A1_4kB[] = { 2, 3, 5, 7, 9, 11, 14, 17, 20, 24, 28, 33, 39, 46, 53, 62, 70, 80, 91, 105, 120, 128 };
static uint8_t  A2_4kB[] = { 2,    4, 6, 8, 10, 14, 18, 22,     28, 34, 40, 48, 56, 66, 74, 84, 94, 104, 120, 128 };
static uint8_t  A4_4kB[] = {       4,    8, 12,     16, 20, 24,     32, 40, 48, 56, 68,     80, 92, 104, 120, 128 };
static uint8_t  A8_4kB[] = {             8,         16,     24,     32, 40, 48,     64,     80, 96, 112, 120, 128 };
static uint8_t A16_4kB[] = {                        16,             32,     48,     64,     80, 96, 112,      128 };

// 128B-1kB blocks, 32k slab
static uint16_t  A1_32kB[] = { 137, 156, 178, 201, 227, 256, 288, 323, 361, 402, 447, 494, 545, 598, 654, 712, 771, 832, 895, 958, 1022 };
static uint16_t  A8_32kB[] = { 144,      168, 192, 224, 256, 296, 336, 376, 424, 472,      528, 584, 640, 704, 768, 832, 896, 960, 1024 };
static uint16_t A16_32kB[] = { 144,      176, 208,      240, 272, 320, 368, 416, 464,      512, 576, 640, 704, 768, 832, 896, 960, 1024 };

// 1kB-8kB blocks, 256k slab
static uint16_t  A1_256kB[] = { 1152, 1297, 1458,       1636, 1832, 2048, 2284,       2541, 2820, 3124, 3550, 3904, 4280, 4676, 5092,       5525, 5974, 6435, 6906, 7380, 7856 };
static uint16_t  A8_256kB[] = { 1152, 1288, 1440,       1608, 1792, 2000, 2224, 2472, 2744, 3032, 3344, 3680, 4040,       4416, 4816, 5232, 5664, 6112, 6568, 7032, 7504, 7976 };
static uint16_t A64_256kB[] = { 1152, 1280, 1408, 1536, 1664, 1856, 2048, 2240, 2432, 2688, 2944, 3200, 3520, 3840, 4160, 4544, 4928, 5312, 5696, 6144, 6592, 7040, 7488, 7936 };

// 8kB-100kB blocks, 2MB slab
static uint32_t A64_2MB[] = { 
    8896,  9984,  11200, 12544, 14016, 15616, 17408, 19328, 21440, 23744, 26176, 28800, 31616, 34624, 37760, 41024, 
    44416, 47936, 51584, 55296, 59008, 62784, 66496, 70208, 73856, 77376, 80832, 84160, 87360, 90368, 93248, 95936,
    98496, 100864
};

int main (void) {

    double x = 100864;
    int n;
    for (n = 0; n < 40 && x < (1 << 21); ++n) {
        x *= 1.1;
        x = (uint32_t)x & ~63;
        printf("%u, ", (uint32_t)x);
    }
    printf("\n%d\n", n);
    return 0;
    const int start1 = 120832;
    const int start2 = 1408;
    const int alignment = 64;
#define ischosen(x) \
    (x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || \
     x == 0 || x == 0 || x == 0 || x == 0 || x == 0 || x == 0)

    const int slab_size = 1 << 21;
    const double thresh = .06;
    int block_size;
    int i = 0;
    for (block_size = start1; i < 87 && block_size < (slab_size >> 3); ++i, block_size += alignment) {
        printf("%5d ", block_size);

        int d;
        double min = 1;
        int ch = block_size + alignment;
        for (d = block_size; d >= alignment; d-=alignment) { 
            int x = block_size - d;
            if (ischosen(x)) {
                double f = calc_frag(slab_size, block_size, d);
                if (f < thresh && f < min) { min = f; ch = d; }
            }
        }

        for (d = start2; d > start2 - 1024; d-=alignment) {
            if (d <= block_size && d <= ch) {
                double f = calc_frag(slab_size, block_size, d);
                if (f < thresh) {
                    if (d == ch) {
                        printf(" *%3.1f%% ", f*100);
                    } else {
                        printf(" %4.1f%% ", f*100);
                    }
                    continue;
                } 
            }
            if (d-1 <= block_size && d-alignment <= ch && calc_frag(slab_size, block_size, d - alignment) < thresh) {
                printf("%6d ", block_size);
                continue;
            }
            printf("       ");
        }
            
        if (ischosen(block_size)) {
            printf("%5d*", block_size);
        } else {
            printf("%5d", block_size);
        }
        printf("\n");
    }
    return 0;
}
