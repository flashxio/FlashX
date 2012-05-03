#ifndef LIST_H
#define LIST_H

#include "map.h"

typedef struct ll list_t;
typedef struct ll_iter ll_iter_t;

list_t *   ll_alloc   (const datatype_t *key_type);
map_val_t  ll_cas     (list_t *ll, map_key_t key, map_val_t expected_val, map_val_t new_val);
map_val_t  ll_lookup  (list_t *ll, map_key_t key);
map_val_t  ll_remove  (list_t *ll, map_key_t key);
size_t     ll_count   (list_t *ll);
void       ll_print   (list_t *ll, int verbose);
void       ll_free    (list_t *ll);
map_key_t  ll_min_key (list_t *sl);

ll_iter_t * ll_iter_begin (list_t *ll, map_key_t key);
map_val_t   ll_iter_next  (ll_iter_t *iter, map_key_t *key_ptr);
void        ll_iter_free  (ll_iter_t *iter);

static const map_impl_t MAP_IMPL_LL = { 
    (map_alloc_t)ll_alloc, (map_cas_t)ll_cas, (map_get_t)ll_lookup, (map_remove_t)ll_remove, 
    (map_count_t)ll_count, (map_print_t)ll_print, (map_free_t)ll_free, (map_iter_begin_t)ll_iter_begin,
    (map_iter_next_t)ll_iter_next, (map_iter_free_t)ll_iter_free
};

#endif//LIST_H
