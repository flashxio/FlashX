/**
 * This file has the unit test for all functions in SA_hash_table.
 */

#include <stdio.h>
#include <assert.h>

#include "SA_hash_table.h"
#include "SA_hash_table.cpp"

void testEntrySet()
{
	typedef entry_set<int, int, 8> entry_set_t;
	entry_set_t set;
	for(int i = 0; i < 8; i++) {
		int ret = set.putIfAbsent(i, i);
		assert(ret == 0);
		printf("put (%d, %d) return %d\n", i, i, ret);
	}
	try {
		set.putIfAbsent(8, 8);
	} catch (no_space_exception e) {
		printf("the entry is out of space\n");
		assert(set.is_overflow());
	}

	for(int i = 0; i < 8; i++) {
		int ret = set.putIfAbsent(i, i + 8);
		assert(ret == i);
		printf("fail to put (%d, %d), original: %d\n",
				i, i + 8, ret);
	}

	for (int i = 0; i < 8; i++) {
		int value = set.get(i);
		assert(value == i);
		printf("the value for key %d is %d\n", value, i);
	}

	for (int i = 0; i < 8; i++) {
		bool ret = set.replace(i, i, -1);
		assert(ret);
		printf("replace for key %d successful? %d\n", i, ret);
		ret = set.replace(i, i, -1);
		assert(!ret);
		printf("replace for key %d successful? %d\n", i, ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = set.remove(i, -1);
		assert(ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = set.remove(i, -1);
		assert(!ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = set.remove(i + 8, -1);
		assert(!ret);
	}
}

void testHashTableBasic()
{
	SA_hashtable<int, int> table(8);
	for(int i = 0; i < 8; i++) {
		int ret = table.putIfAbsent(i, i);
		assert(ret == 0);
		printf("put (%d, %d) return %d\n", i, i, ret);
	}

	for(int i = 0; i < 8; i++) {
		int ret = table.putIfAbsent(i, i + 8);
		assert(ret == i);
		printf("fail to put (%d, %d), original: %d\n",
				i, i + 8, ret);
	}

	for (int i = 0; i < 8; i++) {
		int value = table.get(i);
		assert(value == i);
		printf("the value for key %d is %d\n", value, i);
	}

	for (int i = 0; i < 8; i++) {
		bool ret = table.replace(i, i, -1);
		assert(ret);
		printf("replace for key %d successful? %d\n", i, ret);
		ret = table.replace(i, i, -1);
		assert(!ret);
		printf("replace for key %d successful? %d\n", i, ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = table.remove(i, -1);
		assert(ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = table.remove(i, -1);
		assert(!ret);
	}

	for (int i = 0; i < 4; i++) {
		bool ret = table.remove(i + 8, -1);
		assert(!ret);
	}
}

void testHashTableExpand()
{
	SA_hashtable<int, int> table(8);

	for (int i = 0; i <= 8; i++) {
		int key = i * 16 + 14;
		int value = table.putIfAbsent(key, key);
		assert(value == 0);
	}

	for (int i = 0; i <= 8; i++) {
		int key = i * 16 + 14;
		int value = table.get(key);
		assert(value == key);
	}
}

int main()
{
//	testEntrySet();
//	testHashTableBasic();
	testHashTableExpand();
}
