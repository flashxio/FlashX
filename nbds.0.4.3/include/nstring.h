#ifndef NSTRING_H
#define NSTRING_H

#include "common.h"
#include "datatype.h"

typedef struct nstring {
    uint32_t len;
    char data[];
} nstring_t;

nstring_t * ns_alloc (uint32_t len);
int         ns_cmp   (const nstring_t *ns1, const nstring_t *ns2);
uint32_t    ns_hash  (const nstring_t *ns);
nstring_t * ns_dup   (const nstring_t *ns);

extern const datatype_t DATATYPE_NSTRING;

#endif//NSTRING_H 
