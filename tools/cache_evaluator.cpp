/**
 * This program is to evaluate the cache hit rate given a workload sequence.
 */

#include "workload.h"
#include "cache.h"
#include "associative_cache.h"
#include "hash_index_cache.h"
#include "LRU2Q.h"

#include <string>

int nthreads = 1;

int main(int argc, char *argv[])
{
	if (argc < 4) {
		fprintf(stderr, "cache_evaluator workload_file cache_type cache_size\n");
		fprintf(stderr, "supported cache types: associative, hash-index, lru2q\n");
		return 1;
	}
	std::string workload_file = argv[1];
	std::string cache_type = argv[2];
	// Cache size in bytes.
	int cache_size = atoi(argv[3]) * 4096;

	page_cache *global_cache;
	if (cache_type == "lru2q") {
		global_cache = new LRU2Q_cache(cache_size);
	}
	else if (cache_type == "associative") {
		global_cache = new associative_cache(cache_size, MAX_CACHE_SIZE);
	}
	else if (cache_type == "hash-index") {
		global_cache = new hash_index_cache(cache_size);
	}
	global_cache->init();

	int num_hits = 0;
	int num_accesses = 0;
	workload_gen *gen = new file_workload(workload_file, 1);
	while (gen->has_next()) {
		workload_t workload = gen->next();
		off_t off = ROUND_PAGE(workload.off);
		off_t end = workload.off + workload.size;
		while (off < end) {
			off_t old_off = -1;
			num_accesses++;
			thread_safe_page *pg = (thread_safe_page *) global_cache->search(off, old_off);
			if (old_off == -1)
				num_hits++;
//			pg->set_dirty(false);
			pg->dec_ref();
			off += PAGE_SIZE;
		}
	}
	printf("There are %d accesses and %d hits\n", num_accesses, num_hits);
}
