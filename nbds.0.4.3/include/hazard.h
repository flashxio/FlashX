/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * hazard pointers
 *
 * www.research.ibm.com/people/m/michael/ieeetpds-2004.pdf
 *
 */
#ifndef HAZARD_H
#define HAZARD_H

#define STATIC_HAZ_PER_THREAD 2

typedef void (*free_t) (void *);
typedef void *haz_t;

//static inline void haz_set (volatile haz_t *haz, void *x) { *haz = x; haz_t y = *haz; y = y; }

static inline void haz_set (volatile haz_t *haz, void *x) { *haz = x; __asm__ __volatile__("mfence"); }

haz_t *haz_get_static         (int n);
void   haz_register_dynamic   (haz_t *haz);
void   haz_unregister_dynamic (haz_t *haz);
void   haz_defer_free         (void *p, free_t f);

#endif
