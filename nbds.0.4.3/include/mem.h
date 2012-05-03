/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#ifndef MEM_H
#define MEM_H
void *nbd_malloc (size_t n) __attribute__((malloc, alloc_size(1)));
void nbd_free (void *x) __attribute__((nonnull));
#endif//MEM_H
