/**
 * This program is to evaluate the cache hit rate given a workload sequence.
 */

#include "workload.h"
#include "cache.h"
#include "associative_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"
#include "common.h"

#include <string>

int nthreads = 1;

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "cache_evaluator workload_file cache_type cache_size [cell_size]\n");
		fprintf(stderr, "supported cache types: associative, hash-index, lru2q\n");
		return 1;
	}
	std::string workload_file = argv[1];
	std::string cache_type = argv[2];
	// Cache size in bytes.
	long cache_size = str2size(argv[3]);
	printf("The cache size is %ld bytes\n", cache_size);
	int cell_size = atoi(argv[4]);
	printf("The cell size of SA-cache is %d\n", cell_size);
	params.init(1, cell_size, 0);

	page_cache *global_cache;
	if (cache_type == "lru2q") {
		global_cache = new LRU2Q_cache(cache_size);
	}
	else if (cache_type == "associative") {
		global_cache = associative_cache::create(cache_size, MAX_CACHE_SIZE, 0, 1, -1);
	}
	else if (cache_type == "hash-index") {
		global_cache = new hash_index_cache(cache_size, 0);
	}
	global_cache->init(NULL);

	int num_hits = 0;
	int num_accesses_in_pages = 0;
	int num_accesses = 0;

	long length = 0;
	workload_t *workloads = load_file_workload(workload_file, length);
	file_workload *gen = new file_workload(workloads, length, 0, length);
	long size = length;
	int num_pages = 0;
	printf("There are %ld requests\n", size);
	while (gen->has_next()) {
		workload_t workload = gen->next();
		off_t off = ROUND_PAGE(workload.off);
		off_t end = workload.off + workload.size;
		num_accesses++;
		assert(off < end);
		while (off < end) {
			off_t old_off = -1;
			thread_safe_page *pg = (thread_safe_page *) global_cache->search(off, old_off);
			// We only count the cache hits for the second half of the workload.
			if (num_accesses >= size / 2) {
				num_accesses_in_pages++;
				if (old_off == -1)
					num_hits++;
			}
			pg->dec_ref();
			off += PAGE_SIZE;
		}
		num_pages += (off - ROUND_PAGE(workload.off)) / PAGE_SIZE;
	}
	printf("There are %d accesses in pages in the second half and %d hits, %d total pages\n",
			num_accesses_in_pages, num_hits, num_pages);
}
