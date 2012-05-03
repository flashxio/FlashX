#include "common.h"
#include "nstring.h"
#include "murmur.h"
#include "mem.h"

const datatype_t DATATYPE_NSTRING = { (cmp_fun_t)ns_cmp, (hash_fun_t)ns_hash, (clone_fun_t)ns_dup };

nstring_t *ns_alloc (uint32_t len) {
    nstring_t *ns = nbd_malloc(sizeof(nstring_t) + len);
    ns->len = len;
    return ns;
}

int ns_cmp (const nstring_t *ns1, const nstring_t *ns2) {
    int d = memcmp(ns1->data, ns2->data, (ns1->len < ns2->len) ? ns1->len : ns1->len);
    return (d == 0) ? ns1->len - ns2->len : d;
}

uint32_t ns_hash (const nstring_t *ns) {
    return murmur32(ns->data, ns->len);
}

nstring_t *ns_dup (const nstring_t *ns1) {
    nstring_t *ns2 = ns_alloc(ns1->len);
    memcpy(ns2->data, ns1->data, ns1->len);
    return ns2;
}
