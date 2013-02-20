#ifndef __BELADY_ALGO_H__
#define __BELADY_ALGO_H__

#include <map>
#include <deque>
#include <tr1/unordered_map>
#include <tr1/unordered_set>

class indexed_offset_scanner
{
	int *offs;
	int length;
	// The pages that are only accessed once in the entire sequence.
	std::tr1::unordered_set<int> one_accesses;
	// The last accesses of pages in the sequence.
	// The key is the offset of a page.
	// The value is the location of the last access of the page in the sequence.
	std::tr1::unordered_map<int, int> last_access_map;

	int cached_to;
	int shift;
	// We look forward in the sequence and index the accesses we have
	// looked so far.
	// The key is the offset of a page
	// The value is the locations of the accesses of the page in the sequence.
	std::tr1::unordered_map<int, std::deque<int> > cache;

	void cache_offset(int off, int i);
public:
	indexed_offset_scanner(int offs[], int length);
	int find_next_index(int off) const;
	bool has_next() const {
		return shift < length;
	}
	int next();
	int size() const {
		return length;
	}
};

class belady_algo
{
	int cache_size;	// in the number of pages.
	// Index the accesses in the cache.
	// The key is the offset of an access.
	std::tr1::unordered_map<int, struct access *> access_map;
	// Index the accesses with their locations in the sequence.
	// The key is the location of an access in the sequence.
	// All the locations are unique, so we put them in a tree-based map.
	std::map<int, struct access *> loc_map;
	std::deque<struct access *> non_exist_accesses;

public:
	belady_algo(int cache_size) {
		this->cache_size = cache_size;
	}

	void cache_hit(int off, indexed_offset_scanner &scanner);
	void cache_miss(int off, indexed_offset_scanner &scanner);
	void add_to_cache(int off, indexed_offset_scanner &scanner);
	int access(indexed_offset_scanner &scanner, int hit_count_start);
	void print_cache() const;
	void print_cache_by_seq_order() const;
};

#endif
