/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * A platform independant wrapper around thread-local storage. On platforms that don't support 
 * __thread variables (e.g. Mac OS X), we have to use the pthreads library for thread-local storage 
 */
#ifndef TLS_H
#define TLS_H

#ifdef __ELF__ // use gcc thread-local storage (i.e. __thread variables)
#define DECLARE_THREAD_LOCAL(name, type) __thread type name
#define INIT_THREAD_LOCAL(name) 
#define SET_THREAD_LOCAL(name, value) name = value
#define LOCALIZE_THREAD_LOCAL(name, type)

#else//!__ELF__

#include <pthread.h>

#define DECLARE_THREAD_LOCAL(name, type) pthread_key_t name##_KEY

#define INIT_THREAD_LOCAL(name) \
    do { \
        if (pthread_key_create(&name##_KEY, NULL) != 0) { assert(FALSE); } \
    } while (0)

#define SET_THREAD_LOCAL(name, value) pthread_setspecific(name##_KEY, (void *)(size_t)value);

#define LOCALIZE_THREAD_LOCAL(name, type) type name = (type)(size_t)pthread_getspecific(name##_KEY)

#endif//__ELF__
#endif//TLS_H
