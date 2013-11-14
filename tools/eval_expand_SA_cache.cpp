/**
 * This program is to evaluate the cache hit rate given a workload sequence.
 */

#include "workload.h"
#include "cache.h"
#include "associative_cache.h"

#include <string>

int nthreads = 1;

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "eval_expand_SA_cache workload_file cache_size [new_cache_size]\n");
		return 1;
	}
	std::string workload_file = argv[1];
	// Cache size in bytes.
	int cache_size = atoi(argv[2]) * 4096;
	int new_cache_size = cache_size;
	if (argc == 4)
		new_cache_size = atoi(argv[3]) * 4096;
	assert(new_cache_size >= cache_size);

	associative_cache *global_cache = associative_cache::create(cache_size,
			MAX_CACHE_SIZE, 0, 1, new_cache_size != cache_size);
	global_cache->init(NULL);

	int num_hits = 0;
	int num_accesses_in_pages = 0;
	int num_accesses = 0;

	long length = 0;
	workload_t *workloads = load_file_workload(workload_file, length);
	file_workload *gen = new file_workload(workloads, length, 0, 1);
	int size = (int) length;
	printf("There are %d requests\n", size);
	// We expand the cache by half of its original size.
	if (new_cache_size - cache_size > 0)
		global_cache->expand((new_cache_size - cache_size) / PAGE_SIZE);
	while (gen->has_next()) {
		workload_t workload = gen->next();
		off_t off = ROUND_PAGE(workload.off);
		off_t end = workload.off + workload.size;
		num_accesses++;
		assert(off < end);
		while (off < end) {
			data_loc_t loc(0, off);
			data_loc_t old_loc;
			thread_safe_page *pg = (thread_safe_page *) global_cache->search(loc, old_loc);
			// We only count the cache hits for the second half of the workload.
			if (num_accesses >= size / 2) {
				num_accesses_in_pages++;
				if (old_loc.get_offset() == -1)
					num_hits++;
			}
//			pg->set_dirty(false);
			pg->dec_ref();
			off += PAGE_SIZE;
		}
	}
	global_cache->print_stat();
	printf("There are %d accesses in pages in the second half and %d hits\n",
			num_accesses_in_pages, num_hits);
}
