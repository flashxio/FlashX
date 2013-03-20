#include <vector>
#include <set>

#include "associative_cache.h"

int main()
{
	associative_cache *cache = new associative_cache(16 * 1024 * 1024, 512 * 1024 * 1024, 0, 1);
	std::vector<off_t> added_offs;
	std::set<off_t> evicted_offs;
	for (int i = 0; i < 4096; i++) {
		off_t off = (random() & 0xfffffff) * 4096;
		added_offs.push_back(off);
		off_t old_off = -1;
		thread_safe_page *pg = (thread_safe_page *) cache->search(off, old_off);
		if (old_off != -1)
			evicted_offs.insert(old_off);
		pg->dec_ref();
	}

	int ret;
	ret = cache->expand(8 * 1024 * 1024 / 4096);
	printf("expand %d pages\n", ret);
	cache->sanity_check();
	for (size_t i = 0; i < added_offs.size(); i++) {
		page *pg = cache->search(added_offs[i]);
		if (pg == NULL)
			assert(evicted_offs.find(added_offs[i]) != evicted_offs.end());
		else
			pg->dec_ref();
	}

	for (int i = 0; i < 4096; i++) {
		off_t off = (random() & 0xfffffff) * 4096;
		added_offs.push_back(off);
		off_t old_off = -1;
		thread_safe_page *pg = (thread_safe_page *) cache->search(off, old_off);
		if (old_off != -1)
			evicted_offs.insert(old_off);
		pg->dec_ref();
	}

	ret = cache->expand(16 * 1024 * 1024 / 4096);
	printf("expand %d pages\n", ret);
	cache->sanity_check();
	for (size_t i = 0; i < added_offs.size(); i++) {
		page *pg = cache->search(added_offs[i]);
		if (pg == NULL) {
			if (evicted_offs.find(added_offs[i]) == evicted_offs.end())
				fprintf(stderr, "can't find off %ld\n", added_offs[i]);
			assert(evicted_offs.find(added_offs[i]) != evicted_offs.end());
		}
		else
			pg->dec_ref();
	}

	int num_shrink = 8 * 1024 * 1024 / 4096;
	char *pages[num_shrink];
	cache->shrink(num_shrink, pages);
	cache->sanity_check();
}
