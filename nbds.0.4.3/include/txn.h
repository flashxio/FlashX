/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#ifndef TXN_H
#define TXN_H

#include "map.h"

typedef enum { TXN_RUNNING, TXN_VALIDATING, TXN_VALIDATED, TXN_ABORTED } txn_state_e;

typedef struct txn txn_t;

txn_t *     txn_begin  (map_t *map);
void        txn_abort  (txn_t *txn);
txn_state_e txn_commit (txn_t *txn);

map_val_t   txn_map_get (txn_t *txn, map_key_t key);
void        txn_map_set (txn_t *txn, map_key_t key, map_val_t value);

#endif//TXN_H
