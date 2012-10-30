#ifndef __SA_HASH_TABLE_H__
#define __SA_HASH_TABLE_H__

#include <math.h>
#include <pthread.h>

#include <vector>

#include "hashtable.h"
#include "concurrency.h"

template<class KeyT, class ValueT>
class SA_hashtable;

class no_space_exception
{
};

template<class KeyT, class ValueT, int SIZE>
class entry_set
{
	struct entry {
		KeyT key;
		ValueT value;
	};
	struct entry entries[SIZE];
	int num_entries;
	pthread_spinlock_t lock;
	SA_hashtable<KeyT, ValueT> *table;
	long hash;		// the hash value this set is responsible for.
	bool overflow;

	void init() {
		overflow = false;
		num_entries = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		for (int i = 0; i < SIZE; i++) {
			entries[i].key = KeyT();
			entries[i].value = ValueT();
		}
	}

public:
	entry_set() {
		init();
	}

	entry_set(SA_hashtable<KeyT, ValueT> *table, long hash) {
		this->table = table;
		this->hash = hash;
		init();
	}

	/* 
	 * this is to rehash the pages in the current cell
	 * to the cell in the parameter.
	 */
	void rehash(entry_set *set);

	bool is_overflow() {
		return overflow;
	}

	ValueT get(KeyT key);
	bool remove(KeyT key, ValueT value);
	ValueT putIfAbsent(KeyT key, ValueT value);
	bool replace(KeyT key, ValueT expect, ValueT new_value);
};

static const int SET_SIZE = 8;
template<class KeyT, class ValueT>
class SA_hashtable: public hashtable_interface<KeyT, ValueT>
{
	const int init_nsets;
	const int max_nsets;
	enum {
		TABLE_EXPANDING,
	};
	typedef entry_set<KeyT, ValueT, SET_SIZE> entry_set_t;

	/* 
	 * this lock is to protect the table structure.
	 * if a thread wants to change the table structure,
	 * it needs to grab the lock with write.
	 * otherwise, 
	 */
	seq_lock table_lock;

	std::vector<entry_set_t *> sets_table;
	atomic_integer nsets;
	atomic_flags flags;

	/* used for linear hashing */
	int level;
	int split;

	/********* private members *********/

	entry_set_t *get_set(unsigned int global_idx);
	entry_set_t *get_set_key(KeyT key);
	bool expand(entry_set_t *trigger_set);

public:
	SA_hashtable(int init_nsets);
	~SA_hashtable();

	/* the hash function used for the current level. */
	int hash(KeyT key) {
		return key % (init_nsets * (long) pow(2, level));
	}

	/* the hash function used for the next level. */
	int hash1(KeyT key) {
		return key % (init_nsets * (long) pow(2, level + 1));
	}

	/* interface for the hashtable. */

	ValueT get(KeyT key);
	bool remove(KeyT key, ValueT value);
	/**
	 * If the key doesn't exist, store the key-value pair
	 * and return null. Otherwise, store nothing in the 
	 * table and return the value associated to the key.
	 */
	ValueT putIfAbsent(KeyT key, ValueT value);
	bool replace(KeyT key, ValueT expect, ValueT new_value);
};

#endif
