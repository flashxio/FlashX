#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>

#include "BeladyAlgo.h"
#include "common.h"
#include "workload.h"

void indexed_offset_scanner::cache_offset(int off, int i)
{
	std::tr1::unordered_map<int, std::deque<int> >::iterator it = cache.find(off);
	if (it != cache.end())
		it->second.push_back(i);
	else {
		std::deque<int> queue;
		queue.push_back(i);
		cache.insert(std::pair<int, std::deque<int> >(off, queue));
	}
}

indexed_offset_scanner::indexed_offset_scanner(int offs[], int length)
{
	this->offs = offs;
	this->length = length;
	cached_to = length;
	shift = 0;

	// Find pages that have only one access in the sequence.
	long start = get_curr_ms();
	std::tr1::unordered_map<int, int> access_freq;
	for (int i = 0; i < length; i++) {
		// Can't find the page
		if (access_freq.find(offs[i]) == access_freq.end())
			access_freq.insert(std::pair<int, int>(offs[i], 1));
		else
			access_freq[offs[i]]++;

		// Find the last access of all pages.
		if (last_access_map.find(offs[i]) == last_access_map.end())
			last_access_map.insert(std::pair<int, int>(offs[i], i));
		else
			last_access_map[offs[i]] = i;
	}
	for (std::tr1::unordered_map<int, int>::const_iterator it = access_freq.begin();
			it != access_freq.end(); it++) {
		if (it->second == 1)
			one_accesses.insert(it->first);
	}
	long end = get_curr_ms();
	printf("It takes %ld ms to find one-accesses\n", end - start);
	printf("There are %ld one-accesses, %d accesses and %ld different accesses\n",
			one_accesses.size(), length, access_freq.size());

	// Index the first few accesses.
	start = end;
	for (int i = 0; i < cached_to; i++) {
		int off = offs[i];
		cache_offset(off, i);
	}
	end = get_curr_ms();
	printf("Indexing the first %d accesses takes %ld ms\n", cached_to, end - start);

#ifdef DEBUG
	for (std::tr1::unordered_map<int, std::deque<int> >::const_iterator it = cache.begin();
			it != cache.end(); it++) {
		printf("%d: ", it->first);
		for (unsigned i = 0; i < it->second.size(); i++)
			printf("%d ", it->second[i]);
		printf("\n");
	}
#endif
}

int indexed_offset_scanner::find_next_index(int off) const
{
	if (one_accesses.find(off) != one_accesses.end())
		return length;

	// Check the cache.
	std::tr1::unordered_map<int, std::deque<int> >::const_iterator it1 = cache.find(off);
	if (it1 != cache.end()) {
		return it1->second.front();
	}

	std::tr1::unordered_map<int, int>::const_iterator it = last_access_map.find(off);
	assert(it != last_access_map.end());
	int last_idx = it->second;
	if (last_idx < shift)
		return length;

	for (int i = cached_to; i < length; i++) {
		if (offs[i] == off)
			return i;
	}
	assert(0);
}

int indexed_offset_scanner::next()
{
	int off = offs[shift];
	shift++;
	// Maintain the cache.
	if (one_accesses.find(off) == one_accesses.end()) {
		std::tr1::unordered_map<int, std::deque<int> >::iterator it = cache.find(off);
		it->second.pop_front();
		if (it->second.empty())
			cache.erase(it);
	}
	if (cached_to < length) {
		cache_offset(offs[cached_to], cached_to);
		cached_to++;
	}
	return off;
}

/*
 * It represents the access of a page.
 */
struct access
{
	int seq;
	int off;
	access(int off, int seq) {
		this->seq = seq;
		this->off = off;
	}
};

void belady_algo::cache_hit(int off, indexed_offset_scanner &scanner)
{
	struct access *acc = access_map[off];
	int orig_seq = acc->seq;
	acc->seq = scanner.find_next_index(off);
	loc_map.erase(orig_seq);
	if (acc->seq == scanner.size())
		non_exist_accesses.push_back(acc);
	else
		loc_map.insert(std::pair<int, struct access *>(acc->seq, acc));
}

void belady_algo::cache_miss(int off, indexed_offset_scanner &scanner)
{
	// Get the farthest access
	struct access *acc;
	if (!non_exist_accesses.empty()) {
		acc = non_exist_accesses.front();
		non_exist_accesses.pop_front();
	}
	else {
		std::map<int, struct access *>::reverse_iterator it = loc_map.rbegin();
		acc = it->second;
	}
	print_cache_by_seq_order();
	int next_idx = scanner.find_next_index(off);
	if (next_idx > acc->seq)
		return;

	print_cache_by_seq_order();
	if (acc->seq < scanner.size())
		loc_map.erase(acc->seq);
	access_map.erase(acc->off);
	delete acc;

	acc = new struct access(off, next_idx);
	if (acc->seq == scanner.size())
		non_exist_accesses.push_back(acc);
	else
		loc_map.insert(std::pair<int, struct access *>(acc->seq, acc));
	access_map.insert(std::pair<int, struct access *>(acc->off, acc));
}

void belady_algo::add_to_cache(int off, indexed_offset_scanner &scanner)
{
	int next_idx = scanner.find_next_index(off);
	struct access *acc = new struct access(off, next_idx);
	if (acc->seq == scanner.size())
		non_exist_accesses.push_back(acc);
	else
		loc_map.insert(std::pair<int, struct access *>(acc->seq, acc));
	access_map.insert(std::pair<int, struct access *>(acc->off, acc));
}

int belady_algo::access(indexed_offset_scanner &scanner)
{
	int cache_hits = 0;
	int num_accesses = 0;
	long start = get_curr_ms();
	while (scanner.has_next()) {
		int off = scanner.next();
		print_cache();
		num_accesses++;
		if (num_accesses == 1000000) {
			long end = get_curr_ms();
			printf("%d accesses take %ld ms\n", num_accesses, end - start);
			num_accesses = 0;
			start = end;
		}

		if (access_map.find(off) != access_map.end()) {
			cache_hits++;
			cache_hit(off, scanner);
#ifdef DEBUG
			printf("%d get cache hit\n", off);
#endif
		}
		else if (access_map.size() < cache_size) {
			add_to_cache(off, scanner);
#ifdef DEBUG
			printf("%d gets cache miss, add to cache\n", off);
#endif
		}
		else {
			cache_miss(off, scanner);
#ifdef DEBUG
			printf("%d gets cache miss\n", off);
#endif
		}
	}
	return cache_hits;
}

void belady_algo::print_cache() const
{
#ifdef DEBUG
	for (std::tr1::unordered_map<int, struct access *>::const_iterator it = access_map.begin();
			it != access_map.end(); it++) {
		printf("%d: %d\t", it->second->off, it->second->seq);
	}
	printf("\n");
#endif
}

void belady_algo::print_cache_by_seq_order() const
{
#ifdef DEBUG
	printf("by seq order: ");
	for (std::map<int, struct access *>::const_iterator it = loc_map.begin();
			it != loc_map.end(); it++) {
		printf("%d: %d\t", it->second->off, it->second->seq);
	}
	for (unsigned i = 0; i < non_exist_accesses.size(); i++)
		printf("%d: %d\t", non_exist_accesses[i]->off, non_exist_accesses[i]->seq);
	printf("\n");
#endif
}
