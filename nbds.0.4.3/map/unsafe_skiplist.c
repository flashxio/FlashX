/* 
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 *
 * non thread safe skiplist
 */

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "skiplist.h"
#include "runtime.h"
#include "mem.h"

#define MAX_LEVELS 24

typedef struct node {
    map_key_t key;
    map_val_t val;
    int num_levels;
    struct node *next[1];
} node_t;

struct sl_iter {
    node_t *next;
};

struct sl {
    node_t *head;
    const datatype_t *key_type;
    int high_water; // max level of any item in the list
};

static int random_levels (skiplist_t *sl) {
    uint64_t r = nbd_rand();
    int z = __builtin_ctz(r);
    int levels = (int)(z / 1.5);
    if (levels == 0)
        return 1;
    if (levels > sl->high_water) {
        levels = SYNC_ADD(&sl->high_water, 1);
        TRACE("s2", "random_levels: increased high water mark to %lld", sl->high_water, 0);
    }
    if (levels > MAX_LEVELS) { levels = MAX_LEVELS; }
    return levels;
}

static node_t *node_alloc (int num_levels, map_key_t key, map_val_t val) {
    assert(num_levels > 0 && num_levels <= MAX_LEVELS);
    size_t sz = sizeof(node_t) + (num_levels - 1) * sizeof(node_t *);
    node_t *item = (node_t *)nbd_malloc(sz);
    memset(item, 0, sz);
    item->key = key;
    item->val = val;
    item->num_levels = num_levels;
    TRACE("s2", "node_alloc: new node %p (%llu levels)", item, num_levels);
    return item;
}

skiplist_t *sl_alloc (const datatype_t *key_type) {
    skiplist_t *sl = (skiplist_t *)nbd_malloc(sizeof(skiplist_t));
    sl->key_type = key_type;
    sl->high_water = 1;
    sl->head = node_alloc(MAX_LEVELS, 0, 0);
    memset(sl->head->next, 0, MAX_LEVELS * sizeof(skiplist_t *));
    return sl;
}

void sl_free (skiplist_t *sl) {
    node_t *item = sl->head->next[0];
    while (item) {
        node_t *next = item->next[0];
        if (sl->key_type != NULL) {
            nbd_free((void *)item->key);
        }
        nbd_free(item);
        item = next;
    }
}

size_t sl_count (skiplist_t *sl) {
    size_t count = 0;
    node_t *item = sl->head->next[0];
    while (item) {
        count++;
        item = item->next[0];
    }
    return count;
}

static node_t *find_preds (node_t **preds, node_t **succs, int n, skiplist_t *sl, map_key_t key, int unlink) {
    node_t *pred = sl->head;
    node_t *item = NULL;
    TRACE("s2", "find_preds: searching for key %p in skiplist (head is %p)", key, pred);
    int d = 0;

    // Traverse the levels of <sl> from the top level to the bottom
    for (int level = sl->high_water - 1; level >= 0; --level) {
        node_t *next = pred->next[level];
        if (next == DOES_NOT_EXIST && level >= n)
            continue;
        TRACE("s3", "find_preds: traversing level %p starting at %p", level, pred);
        item = next;
        while (item != NULL) {
            next = item->next[level];

            if (EXPECT_TRUE(sl->key_type == NULL)) {
                d = item->key - key;
            } else {
                d = sl->key_type->cmp((void *)item->key, (void *)key);
            }

            if (d >= 0) {
                if (d == 0 && unlink) {
                    pred->next[level] = next;
                    TRACE("s3", "find_preds: unlinked item from pred %p", pred, 0);
                    item = next;
                    next = (item != NULL) ? item->next[level] : DOES_NOT_EXIST;
                }
                break;
            }

            pred = item;
            item = next;
        }

        TRACE("s3", "find_preds: found pred %p next %p", pred, item);

        if (level < n) { 
            if (preds != NULL) {
                preds[level] = pred;
            }
            if (succs != NULL) {
                succs[level] = item;
            }
        }
    }

    if (d == 0) {
        TRACE("s2", "find_preds: found matching item %p in skiplist, pred is %p", item, pred);
        return item;
    }
    TRACE("s2", "find_preds: found proper place for key %p in skiplist, pred is %p. returning null", key, pred);
    return NULL;
}

// Fast find that does not return the node's predecessors.
map_val_t sl_lookup (skiplist_t *sl, map_key_t key) {
    TRACE("s1", "sl_lookup: searching for key %p in skiplist %p", key, sl);
    node_t *item = find_preds(NULL, NULL, 0, sl, key, FALSE);

    // If we found an <item> matching the <key> return its value.
    if (item != NULL) {
        map_val_t val = item->val;
        return val;
    }

    TRACE("s1", "sl_lookup: no item in the skiplist matched the key", 0, 0);
    return DOES_NOT_EXIST;
}

map_key_t sl_min_key (skiplist_t *sl) {
    node_t *item = sl->head->next[0];
    while (item != NULL)
        return item->key;
    return DOES_NOT_EXIST;
}

map_val_t sl_cas (skiplist_t *sl, map_key_t key, map_val_t expectation, map_val_t new_val) {
    TRACE("s1", "sl_cas: key %p skiplist %p", key, sl);
    TRACE("s1", "sl_cas: expectation %p new value %p", expectation, new_val);
    ASSERT((int64_t)new_val > 0);

    node_t *preds[MAX_LEVELS];
    node_t *nexts[MAX_LEVELS];
    node_t *new_item = NULL;
    int n = random_levels(sl);
    node_t *old_item = find_preds(preds, nexts, n, sl, key, FALSE);

    // If there is already an item in the skiplist that matches the key just update its value.
    if (old_item != NULL) {
        map_val_t old_val = old_item->val;
        if (expectation == CAS_EXPECT_DOES_NOT_EXIST || 
           (expectation != CAS_EXPECT_WHATEVER && expectation != CAS_EXPECT_EXISTS && expectation != old_val)) {
            TRACE("s1", "sl_cas: the expectation was not met; the skiplist was not changed", 0, 0);
            return old_val;
        } 
        old_item->val = new_val;
        return old_val;
    }

    if (EXPECT_FALSE(expectation != CAS_EXPECT_DOES_NOT_EXIST && expectation != CAS_EXPECT_WHATEVER)) {
        TRACE("s1", "sl_cas: the expectation was not met, the skiplist was not changed", 0, 0);
        return DOES_NOT_EXIST; // failure, the caller expected an item for the <key> to already exist 
    }

    TRACE("s3", "sl_cas: inserting a new item between %p and %p", preds[0], nexts[0]);

    // Create a new node and insert it into the skiplist.
    map_key_t new_key = sl->key_type == NULL ? key : (map_key_t)sl->key_type->clone((void *)key);
    new_item = node_alloc(n, new_key, new_val);

    // Set <new_item>'s next pointers to their proper values
    for (int level = 0; level < new_item->num_levels; ++level) {
        new_item->next[level] = nexts[level];
    }

    // Link <new_item> into <sl> 
    for (int level = 0; level < new_item->num_levels; ++level) {
        preds[level]->next[level] = new_item;
    }

    return DOES_NOT_EXIST; // success, inserted a new item
}

map_val_t sl_remove (skiplist_t *sl, map_key_t key) {
    TRACE("s1", "sl_remove: removing item with key %p from skiplist %p", key, sl);
    node_t *preds[MAX_LEVELS];
    node_t *item = find_preds(preds, NULL, sl->high_water, sl, key, FALSE);
    if (item == NULL) {
        TRACE("s3", "sl_remove: remove failed, an item with a matching key does not exist in the skiplist", 0, 0);
        return DOES_NOT_EXIST;
    }
    map_val_t val = item->val; 

    // unlink the item
    find_preds(NULL, NULL, 0, sl, key, TRUE);

    // free the node
    if (sl->key_type != NULL) {
        nbd_free((void *)item->key);
    }
    nbd_free(item);

    return val;
}

void sl_print (skiplist_t *sl) {

    printf("high water: %d levels\n", sl->high_water);
    for (int level = MAX_LEVELS - 1; level >= 0; --level) {
        node_t *item = sl->head;
        if (item->next[level] == DOES_NOT_EXIST)
            continue;
        printf("(%d) ", level);
        int i = 0;
        while (item) {
            node_t *next = item->next[level];
            printf("%p ", item);
            item = next;
            if (i++ > 30) {
                printf("...");
                break;
            }
        }
        printf("\n");
        fflush(stdout);
    }
    node_t *item = sl->head;
    int i = 0;
    while (item) {
        printf("%p:0x%llx ", item, (uint64_t)item->key);
        if (item != sl->head) {
            printf("[%d]", item->num_levels);
        } else {
            printf("[HEAD]");
        }
        for (int level = 1; level < item->num_levels; ++level) {
            node_t *next = item->next[level];
            printf(" %p", next);
            if (item == sl->head && item->next[level] == DOES_NOT_EXIST)
                break;
        }
        printf("\n");
        fflush(stdout);
        item = item->next[0];
        if (i++ > 30) {
            printf("...\n");
            break;
        }
    }
}

sl_iter_t *sl_iter_begin (skiplist_t *sl, map_key_t key) {
    sl_iter_t *iter = (sl_iter_t *)nbd_malloc(sizeof(sl_iter_t));
    if (key != DOES_NOT_EXIST) {
        find_preds(NULL, &iter->next, 1, sl, key, FALSE);
    } else {
        iter->next = sl->head->next[0];
    }
    return iter;
}

map_val_t sl_iter_next (sl_iter_t *iter, map_key_t *key_ptr) {
    assert(iter);
    node_t *item = iter->next;
    if (item == NULL) {
        iter->next = NULL;
        return DOES_NOT_EXIST;
    }
    iter->next = item->next[0];
    if (key_ptr != NULL) {
        *key_ptr = item->key;
    }
    return item->val;
}

void sl_iter_free (sl_iter_t *iter) {
    nbd_free(iter);
}
