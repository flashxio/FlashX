#ifndef DATATYPE_H
#define DATATYPE_H

#include "common.h"

typedef int      (*cmp_fun_t)   (void *, void *);
typedef void *   (*clone_fun_t) (void *);
typedef uint32_t (*hash_fun_t)  (void *);

typedef struct datatype {
    cmp_fun_t   cmp;
    hash_fun_t  hash;
    clone_fun_t clone;
} datatype_t;

#endif//DATATYPE_H
