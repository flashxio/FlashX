/*
 * Written by Josh Dybnis and released to the public domain, as explained at
 * http://creativecommons.org/licenses/publicdomain
 */
#include "common.h"
#include "txn.h"
#include "mem.h"
#include "rcu.h"
#include "lwt.h"
#include "skiplist.h"

#define UNDETERMINED_VERSION 0
#define ABORTED_VERSION      TAG_VALUE(0, TAG1)
#define INITIAL_WRITES_SIZE  4
#define PTR_TO_VAL(x) ((size_t)(x) >> 2)
#define VAL_TO_PTR(x) ((update_t *)((x) << 2))

typedef struct update_rec update_t;
typedef map_key_t version_t;

struct update_rec {
    version_t version; // tagged versions are txn_t pointers, untagged are actual version numbers
    map_val_t value;
    map_val_t next; // an earlier update
};

typedef struct write_rec {
    map_key_t key;
    update_t *rec;
} write_rec_t;

struct txn {
    version_t rv;
    version_t wv;
    map_t *map;
    write_rec_t *writes;
    size_t writes_size;
    size_t writes_count;
    size_t validate_scan;
    txn_state_e state;
};

static txn_state_e txn_validate (txn_t *txn);

static skiplist_t *active_ = NULL;

static version_t version_ = 1;

// Validate the updates for <key>. Validation fails if there is a write-write conflict. That is if after our
// read version another transaction committed a change to an entry we are also trying to change.
//
// If we encounter a potential conflict with a transaction that is in the process of validating, we help it
// complete validating. It must be finished before we can decide to rollback or commit.
//
static txn_state_e validate_key (txn_t *txn, map_key_t key) {
    assert(txn->state != TXN_RUNNING);

    map_val_t val = map_get(txn->map, key);
    update_t *update = NULL;
    for (; val != DOES_NOT_EXIST; val = update->next) {

        // If the update or its version is not tagged it means the update is committed.
        //
        // We can stop at the first committed record we find that is at least as old as our read version. All
        // the other committed records following it will be older. And all the uncommitted records following it
        // will eventually conflict with it and abort.
        if (!IS_TAGGED(val, TAG2))
            return TXN_VALIDATED;
        update = VAL_TO_PTR(val);
        if (!IS_TAGGED(update->version, TAG1))
            return (update->version <= txn->rv) ? TXN_VALIDATED : TXN_ABORTED;

        // If the update's version is tagged then either the update was aborted or the the version number is
        // actually a pointer to a running transaction's txn_t.

        // Skip aborted transactions.
        if (EXPECT_FALSE(update->version == ABORTED_VERSION))
            continue;

        // The update's transaction is still in progress. Access its txn_t.
        txn_t *writer = (txn_t *)VAL_TO_PTR(update->version);
        if (writer == txn)
            continue; // Skip our own updates.
        txn_state_e writer_state = writer->state;

        // Any running transaction will only be able to acquire a wv greater than ours. A transaction changes its
        // state to validating before aquiring a wv. We can ignore an unvalidated transaction if its version is
        // greater than ours. See the next comment below for the explination why.
        if (writer_state == TXN_RUNNING)
            continue;

        // If <writer> has a later version than us we can safely ignore its updates. It will not commit until
        // we have completed validation (in order to remain non-blocking it will help us validate if necessary).
        // This protocol ensures a deterministic resolution to every conflict and avoids infinite ping-ponging
        // between validating two conflicting transactions.
        if (writer_state == TXN_VALIDATING) {
            if (writer->wv > txn->wv)
                continue;
            // Help <writer> commit. We need to know if <writer> aborts or commits before we can decide what to
            // do. But we don't want to block, so we assist.
            writer_state = txn_validate(writer);
        }

        // Skip updates from aborted transactions.
        if (writer_state == TXN_ABORTED)
            continue;

        assert(writer_state == TXN_VALIDATED);
        return (writer->wv <= txn->rv) ? TXN_VALIDATED : TXN_ABORTED;
    }

    return TXN_VALIDATED;
}

static txn_state_e txn_validate (txn_t *txn) {
    assert(txn->state != TXN_RUNNING);
    switch (txn->state) {

        case TXN_VALIDATING:
            if (txn->wv == UNDETERMINED_VERSION) {
                version_t wv = SYNC_ADD(&version_, 1);
                (void)SYNC_CAS(&txn->wv, UNDETERMINED_VERSION, wv);
            }

            for (int i = 0; i < txn->writes_count; ++i) {
                txn_state_e s = validate_key(txn, txn->writes[i].key);
                if (s == TXN_ABORTED) {
                    txn->state = TXN_ABORTED;
                    break;
                }
                assert(s == TXN_VALIDATED);
            }
            if (txn->state == TXN_VALIDATING) {
                txn->state =  TXN_VALIDATED;
            }
            break;

        case TXN_VALIDATED:
        case TXN_ABORTED:
            break;

        default:
            assert(FALSE);
    }

    return txn->state;
}

static update_t *alloc_update_rec (version_t ver, map_val_t val) {
    update_t *u = (update_t *)nbd_malloc(sizeof(update_t));
    u->version = ver;
    u->value = val;
    u->next = DOES_NOT_EXIST;
    return u;
}

txn_t *txn_begin (map_t *map) {
    TRACE("x1", "txn_begin: map %p", map, 0);
    txn_t *txn = (txn_t *)nbd_malloc(sizeof(txn_t));
    memset(txn, 0, sizeof(txn_t));
    txn->wv = UNDETERMINED_VERSION;
    txn->state = TXN_RUNNING;
    txn->map = map;
    txn->writes = nbd_malloc(sizeof(*txn->writes) * INITIAL_WRITES_SIZE);
    txn->writes_size = INITIAL_WRITES_SIZE;
    if (EXPECT_FALSE(active_ == NULL)) {
        skiplist_t *a = sl_alloc(NULL);
        if (SYNC_CAS(&active_, NULL, a) != NULL) {
            sl_free(a);
        }
    }

    // acquire the read version for txn. must be careful to avoid a race
    do {
        txn->rv = version_;

        unsigned old_count;
        unsigned temp = 0;
        do {
            old_count = temp;
            temp = sl_cas(active_, txn->rv, old_count, old_count + 1);
        } while (temp != old_count);

        if (txn->rv == version_)
            break;

        temp = 1;
        do {
            old_count = temp;
            temp = sl_cas(active_, (map_key_t)txn->rv, old_count, old_count - 1);
        } while (temp != old_count);
    } while (1);

    TRACE("x1", "txn_begin: returning new transaction %p (read version %p)", txn, txn->rv);
    return txn;
}

void txn_abort (txn_t *txn) {
    if (txn->state != TXN_RUNNING)
        return;

    int i;
    for (i = 0; i < txn->writes_count; ++i) {
        update_t *update = (update_t *)txn->writes[i].rec;
        update->version = ABORTED_VERSION;
    }

    rcu_defer_free(txn->writes);
    rcu_defer_free(txn);
}

txn_state_e txn_commit (txn_t *txn) {
    if (txn->state != TXN_RUNNING)
        return txn->state;

    assert(txn->state == TXN_RUNNING);
    txn->state = TXN_VALIDATING;
    txn_state_e state = txn_validate(txn);

    // Detach <txn> from its updates.
    version_t wv = (txn->state == TXN_ABORTED) ? ABORTED_VERSION : txn->wv;
    int i;
    for (i = 0; i < txn->writes_count; ++i) {
        update_t *update = txn->writes[i].rec;
        update->version = wv;
    }

    // Lower the reference count for <txn>'s read version
    unsigned temp = 2;
    unsigned old_count;
    do {
        old_count = temp;
        temp = sl_cas(active_, (map_key_t)txn->rv, old_count, old_count - 1);
        if (temp == 1 && txn->rv != version_) {
            sl_remove(active_, (map_key_t)txn->rv);
            break;
        }
    } while (old_count != temp);

    rcu_defer_free(txn->writes);
    rcu_defer_free(txn);

    return state;
}

// Get most recent committed version prior to our read version.
map_val_t txn_map_get (txn_t *txn, map_key_t key) {
    TRACE("x1", "txn_map_get: txn %p map %p", txn, txn->map);
    TRACE("x1", "txn_map_get: key %p", key, 0);

    if (txn->state != TXN_RUNNING) {
        TRACE("x1", "txn_map_get: error txn not running (state %p)", txn->state, 0);
        return ERROR_TXN_NOT_RUNNING;
    }

    // Iterate through the update records to find the latest committed version prior to our read version.
    map_val_t newest_val = map_get(txn->map, key);
    map_val_t val = newest_val;
    update_t *update;
    for ( ; (update = VAL_TO_PTR(val)) != NULL ; val = update->next) {

        // If TAG2 is set in <val> it indicates that <val> is an update record. Otherwise all the following are
        // true: <val> is a literal value, it is older than any currently active transaction, and it is the most
        // recently set value for its key. Therefore it is visible to <txn>.
        if (!IS_TAGGED(val, TAG2)) {
            TRACE("x1", "txn_map_get: found untagged value; returning %p", val, 0);
            return val;
        }

        // If the update's version is not tagged it means the update is committed.
        if (!IS_TAGGED(update->version, TAG1)) {
            if (update->version <= txn->rv) {
                TRACE("x2", "txn_map_get: found committed update %p (version %p)", update, update->version);
                break; // success
            }
            TRACE("x2", "txn_map_get: skipping update %p (version %p)", update, update->version);
            continue;
        }

        // If the update's version is tagged then either the update was aborted or the the version number is
        // actually a pointer to a running transaction's txn_t.

        // Skip updates from aborted transactions.
        if (EXPECT_FALSE(update->version == ABORTED_VERSION)) {
            TRACE("x2", "txn_map_get: skipping aborted update %p", update, 0);
            continue;
        }

        // The update's transaction is still in progress. Access its txn_t.
        txn_t *writer = (txn_t *)VAL_TO_PTR(update->version);
        if (writer == txn) { 
            TRACE("x2", "txn_map_get: found txn's own update %p", update, 0);
            break; // success
        }

        txn_state_e writer_state = writer->state;
        if (writer_state == TXN_RUNNING) {
            TRACE("x2", "txn_map_get: skipping update %p of in-progress transaction %p", update, writer);
            continue;
        }

        if (writer_state == TXN_VALIDATING) {
            TRACE("x2", "txn_map_get: update %p transaction %p validating", update, writer);
            if (writer->wv > txn->rv)
                continue;
            writer_state = txn_validate(writer);
        }

        // Skip updates from aborted transactions.
        if (writer_state == TXN_ABORTED) {
            TRACE("x2", "txn_map_get: skipping aborted update %p", update, 0);
            continue;
        }

        assert(writer_state == TXN_VALIDATED);
        if (writer->wv > txn->rv) {
            TRACE("x2", "txn_map_get: skipping update %p (version %p)", update, update->version);
            continue;
        }
        break; // success
    }

    if (update == NULL) {
        TRACE("x1", "txn_map_get: key does not exist in map", key, 0);
        return DOES_NOT_EXIST;
    }

    map_val_t value = update->value;
    TRACE("x1", "txn_map_get: key found returning value %p", value, 0);

    // collect some garbage
    version_t min_active_version = UNDETERMINED_VERSION;
    update_t *next_update = NULL;
    if (IS_TAGGED(update->next, TAG2)) {
        next_update = VAL_TO_PTR(update->next);

        // If <next_update> (and all update records following it [execpt if it is aborted]) is old enough
        // that it is not visible to any active transaction we can safely free it.
        min_active_version = (version_t)sl_min_key(active_);
        if (next_update->version < min_active_version) {

            // If the <next_update> is aborted, skip over it to look for more recent ones that may follow
            update_t *temp = next_update;
            while (temp->version == ABORTED_VERSION) {
                assert(!IS_TAGGED(temp->version, TAG1));
                map_val_t next = temp->next;
                if (!IS_TAGGED(next, TAG2))
                    break;

                // Bail out of garbage collection if we find a record that might still be accessed by an
                // ongoing transaction.
                if (VAL_TO_PTR(next)->version >= min_active_version)
                    return value;

                temp = VAL_TO_PTR(next);
            }

            // free the next update record and all the ones following it
            temp = next_update;
            map_val_t next;
            do {
                next = SYNC_SWAP(&temp->next, DOES_NOT_EXIST);

                // if we find ourself in a race just back off and let the other thread take care of it
                if (next == DOES_NOT_EXIST)
                    return value;

                nbd_free(temp);

                temp = VAL_TO_PTR(next);

            } while (IS_TAGGED(next, TAG2));
        }
    }

    // If there is one item left and it is visible by all active transactions we can merge it into the map itself.
    // There is no need for an update record.
    if (next_update == NULL && val == newest_val) {
        if (min_active_version == UNDETERMINED_VERSION) {
            min_active_version = (version_t)sl_min_key(active_);
        }
        if (update->version <= min_active_version) {
            if (map_cas(txn->map, key, TAG_VALUE(val, TAG2), value) == TAG_VALUE(val, TAG2)) {
                rcu_defer_free(update);
            }
        }
    }

    return value;
}

void txn_map_set (txn_t *txn, map_key_t key, map_val_t value) {
    TRACE("x1", "txn_map_set: txn %p map %p", txn, txn->map);
    TRACE("x1", "txn_map_set: key %p value %p", key, value);
    assert(!IS_TAGGED(value, TAG1) && !IS_TAGGED(value, TAG2));

    if (txn->state != TXN_RUNNING) {
        TRACE("x1", "txn_map_set: error txn not running (state %p)", txn->state, 0);
        return;
    }

    // create a new update record
    version_t ver = TAG_VALUE(PTR_TO_VAL(txn), TAG1); // tagged versions are txn_t pointers
    update_t *update = alloc_update_rec(ver, value);

    // push the new update record onto <key>'s update list
    map_val_t old_update = map_get(txn->map, key);
    TRACE("x2", "txn_map_set: old update %p new update record %p", old_update, update);
    do {
        update->next = old_update;
        map_val_t temp = map_cas(txn->map, key, old_update, TAG_VALUE(PTR_TO_VAL(update), TAG2));
        if (temp == old_update) 
            break;

        TRACE("x1", "txn_map_set: cas failed; found %p expected %p", temp, old_update);
        old_update = temp;
    } while (1);

    // add <key> to the write set for commit-time validation
    if (txn->writes_count == txn->writes_size) {
        write_rec_t *w = nbd_malloc(sizeof(write_rec_t) * txn->writes_size * 2);
        memcpy(w, txn->writes, txn->writes_size * sizeof(write_rec_t));
        txn->writes_size *= 2;
        nbd_free(txn->writes);
        txn->writes = w;
    }
    int i = txn->writes_count++;
    txn->writes[i].key = key;
    txn->writes[i].rec = update;
}
