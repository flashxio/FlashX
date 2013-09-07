#ifndef __SA_HASH_TABLE_H__
#define __SA_HASH_TABLE_H__

#include <math.h>
#include <pthread.h>

#include <iostream>
#include <vector>

#include "my_hashtable.h"
#include "concurrency.h"
#include "exception.h"

template<class KeyT, class ValueT>
class SA_hashtable;

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
	atomic_flags<int> flags;

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
		return key % (init_nsets * (1 << level));
	}

	/* the hash function used for the next level. */
	int hash1(KeyT key) {
		return key % (init_nsets * (1 << (level + 1)));
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

template<class KeyT, class ValueT, int SIZE>
ValueT entry_set<KeyT, ValueT, SIZE>::get(KeyT key) {
	pthread_spin_lock(&lock);
	for (int i = 0; i < num_entries; i++) {
		if (entries[i].key == key) {
			ValueT ret = entries[i].value;
			pthread_spin_unlock(&lock);
			return ret;
		}
	}
	pthread_spin_unlock(&lock);
	return ValueT();
}

template<class KeyT, class ValueT, int SIZE>
bool entry_set<KeyT, ValueT, SIZE>::remove(KeyT key, ValueT value) {
	bool ret = false;
	pthread_spin_lock(&lock);
	for (int i = 0; i < num_entries; i++) {
		if (entries[i].key == key && entries[i].value == value) {
			/* if it's not the last entry */
			if (i < num_entries - 1)
				memmove(&entries[i], &entries[i + 1],
						sizeof(entries[0]) * (num_entries - i - 1));
			num_entries--;
			ret = true;
			break;
		}
	}
	pthread_spin_unlock(&lock);
	return ret;
}

template<class KeyT, class ValueT, int SIZE>
ValueT entry_set<KeyT, ValueT, SIZE>::putIfAbsent(KeyT key, ValueT value) {
	pthread_spin_lock(&lock);
	int i;
	for (i = 0; i < num_entries; i++) {
		if (entries[i].key == key) {
			break;
		}
	}
	/* if the key isn't in the set and there is space for the key-value */
	if (i == num_entries && num_entries < SIZE) {
		num_entries++;
		entries[i].key = key;
		entries[i].value = value;
		pthread_spin_unlock(&lock);
		return ValueT();
	}
	/* if the key is already in the set. */
	else if (i < num_entries) {
		pthread_spin_unlock(&lock);
		return entries[i].value;
	}
	pthread_spin_unlock(&lock);

	/* now we don't have enough space */
	overflow = true;
	throw no_space_exception();
}

template<class KeyT, class ValueT, int SIZE>
bool entry_set<KeyT, ValueT, SIZE>::replace(KeyT key,
		ValueT expect, ValueT new_value) {
	bool ret = false;
	pthread_spin_lock(&lock);
	for (int i = 0; i < num_entries; i++) {
		if (entries[i].key == key && entries[i].value == expect) {
			entries[i].value = new_value;
			ret = true;
			break;
		}
	}
	pthread_spin_unlock(&lock);
	return ret;
}

/**
 * rehash the pages in the current cell to the expanded cell.
 *
 * This function is called when a thread is expanding the hashtable,
 * and only one thread can expand the hashtable.
 */
template<class KeyT, class ValueT, int SIZE>
void entry_set<KeyT, ValueT, SIZE>::rehash(entry_set *expanded) {
	pthread_spin_lock(&lock);
	pthread_spin_lock(&expanded->lock);
	for (int i = 0, j = 0; i < num_entries; ) {
		/* 
		 * It's guaranteed the size of the hashtable can't be changed
		 * by other threads. 
		 */
		struct entry *e = &entries[i];
		int hash1 = table->hash1(e->key);
		/* 
		 * if the two hash values don't match,
		 * it means the entry is mapped to the expanded cell,
		 * so we add it to the expanded set.
		 */
		if (this->hash != hash1) {
			struct entry tmp = *e;
			entries[i] = entries[num_entries - 1];
			num_entries--;
			expanded->entries[j] = tmp;
			expanded->num_entries++;
			j++;
		}
		else
			i++;
	}
	pthread_spin_unlock(&expanded->lock);
	pthread_spin_unlock(&lock);
	if (num_entries < SIZE)
		overflow = false;
}

template<class KeyT, class ValueT>
entry_set<KeyT, ValueT, SET_SIZE> *SA_hashtable<KeyT, ValueT>::get_set_key(KeyT key)
{
	int global_idx;
	unsigned long count;
	entry_set_t *set = NULL;
	do {
		table_lock.read_lock(count);
		global_idx = hash(key);
		if (global_idx < split)
			global_idx = hash1(key);
		set = get_set(global_idx);
	} while (!table_lock.read_unlock(count));
	assert(set);
	return set;
}

template<class KeyT, class ValueT>
entry_set<KeyT, ValueT, SET_SIZE> *SA_hashtable<KeyT, ValueT>::get_set(
		unsigned int global_idx)
{
	unsigned int set_idx = global_idx / init_nsets;
	int idx = global_idx % init_nsets;
	assert(set_idx < sets_table.size());
	entry_set_t *sets = sets_table[set_idx];
	if (sets)
		return &sets[idx];
	else
		return NULL;
}

template<class KeyT, class ValueT>
SA_hashtable<KeyT, ValueT>::SA_hashtable(int init_nsets_):
	init_nsets(init_nsets_), max_nsets(init_nsets_ * 1024)
{
	level = 0;
	split = 0;
	entry_set_t *sets = new entry_set_t[init_nsets];
	for (int i = 0; i < init_nsets; i++)
		sets[i] = entry_set_t(this, i);

	sets_table.push_back(sets);
	nsets.inc(1);

	for (int i = 1; i < max_nsets / init_nsets; i++)
		sets_table.push_back(NULL);
}

template<class KeyT, class ValueT>
SA_hashtable<KeyT, ValueT>::~SA_hashtable() {
	for (int i = 0; i < max_nsets / init_nsets; i++)
		delete [] sets_table[i];
}

template<class KeyT, class ValueT>
ValueT SA_hashtable<KeyT, ValueT>::get(KeyT key) {
	ValueT v;
	unsigned long count;
	do {
		table_lock.read_lock(count);
		v = get_set_key(key)->get(key);
	} while (!table_lock.read_unlock(count));
	return v;
}

template<class KeyT, class ValueT>
bool SA_hashtable<KeyT, ValueT>::remove(KeyT key, ValueT value) {
	bool ret = false;
	unsigned long count;
	do {
		table_lock.read_lock(count);
		bool tmp = get_set_key(key)->remove(key, value);
		/* 
		 * it's possible that a key is removed from the hashtable
		 * successfully, but the table structure is changed.
		 * we can still return immediately.
		 */
		if (tmp) {
			table_lock.read_unlock(count);
			return true;
		}
	} while (!table_lock.read_unlock(count));
	return ret;
}

/**
 * expand the hashtable.
 * @trigger_set: the set that triggers the hashtable expansion.
 */
template<class KeyT, class ValueT>
bool SA_hashtable<KeyT, ValueT>::expand(entry_set_t *trigger_set) {
	entry_set_t *sets = NULL;
	unsigned int i;

	/*
	 * if the flag has been set before, it means another thread
	 * is expanding the table, wait for it to finish expanding.
	 */
	while (flags.set_flag(TABLE_EXPANDING)) {
	}

	/* starting from this point, only one thred can be here. */

	/* find the array of sets with the triggered set. */
	for (i = 0; i < sets_table.size(); i++) {
		sets = sets_table[i];
		if (sets == NULL)
			break;
		if (trigger_set >= sets && trigger_set < sets + init_nsets)
			break;
	}
	assert(sets);

	entry_set_t *set = get_set(split);
	long size = (1 << level) * init_nsets;
	while (trigger_set->is_overflow()) {
		unsigned int sets_idx = (split + size) / init_nsets;
		/* 
		 * I'm sure only this thread can change the table,
		 * so it doesn't need to hold a lock when accessing the size.
		 */
		unsigned int orig_size = nsets.get();
		if (sets_idx >= orig_size) {
			bool out_of_memory = false;
			/* create sets and put them in a temporary table. */
			std::vector<entry_set_t *> table;
			for (unsigned int i = orig_size; i <= sets_idx; i++) {
				entry_set_t *sets = new entry_set_t[init_nsets];
				try {
					for (int j = 0; j < init_nsets; j++) {
						sets[j] = entry_set_t(this, i * init_nsets + j);
					}
					table.push_back(sets);
				} catch (oom_exception e) {
					out_of_memory = true;
					delete [] sets;
					break;
				}
			}

			/*
			 * here we need to hold the lock because other threads
			 * might be accessing the table. by using the write lock,
			 * we notify others the table has been changed.
			 */
			table_lock.write_lock();
			for (unsigned int i = 0; i < table.size(); i++) {
				sets_table[orig_size + i] = table[i];
			}
			nsets.inc(table.size());
			table_lock.write_unlock();
			if (out_of_memory)
				return false;
		}

		entry_set_t *expanded_set = get_set(split + size);
		set->rehash(expanded_set);
		table_lock.write_lock();
		split++;
		if (split == size) {
			level++;
			std::cout << "increase level to " << level << std::endl;
			split = 0;
			table_lock.write_unlock();
			break;
		}
		table_lock.write_unlock();
		set = get_set(split);
	}
	flags.clear_flag(TABLE_EXPANDING);
	return true;
}

/**
 * put the key-value pair to the table. if the key has existed in the table,
 * don't put the new value in the table and return the existing value,
 * otherwise, return a null value.
 */
template<class KeyT, class ValueT>
ValueT SA_hashtable<KeyT, ValueT>::putIfAbsent(KeyT key, ValueT value) {
	ValueT ret = ValueT();
	unsigned long count;
	entry_set_t *prev = NULL;
	do {
		table_lock.read_lock(count);
		entry_set_t *set = get_set_key(key);

		/*
		 * If the previous iteration puts the key-value pair
		 * to a wrong set, we need to remove it from the wrong
		 * set, and do it again.
		 */
		if (prev && prev != set && ret == ValueT()) {
			/* 
			 * It's not guaranteed the key is still there.
			 * It might have to removed by another thread
			 * even though the possibility is very small.
			 */
			prev->remove(key, value);
		}

		try {
			ret = set->putIfAbsent(key, value);
			prev = set;
		} catch (no_space_exception e) {
			expand(set);
		}
	} while (!table_lock.read_unlock(count));
	return ret;
	
}

template<class KeyT, class ValueT>
bool SA_hashtable<KeyT, ValueT>::replace(KeyT key, ValueT expect, ValueT new_value) {
	bool ret;
	unsigned long count;
	do {
		table_lock.read_lock(count);
		entry_set_t *set = get_set_key(key);
		ret = set->replace(key, expect, new_value);
		/* if we can replace the old value successfully, we are done. */
		if (ret) {
			table_lock.read_unlock(count);
			return true;
		}
	} while (!table_lock.read_unlock(count));
	return ret;
}

#endif
