#ifndef __PART_CACHED_PRIVATE_H__
#define __PART_CACHED_PRIVATE_H__

#include "direct_private.h"

class part_cached_private: public direct_private
{
	page_cache *cache;

public:
	part_cached_private(const char *names[], int num, long size, int idx, long cache_size,
			int entry_size): direct_private(names, num, size, idx, entry_size) {
		/* all local cache has the same size */
		cache = new tree_cache(cache_size, idx * cache_size);
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		ssize_t ret;
		off_t old_off = -1;
		assert(access_method == READ);
		page *p = cache->search(ROUND_PAGE(offset), old_off);

		assert(p->get_offset() == ROUND_PAGE(offset));
		if (!p->data_ready()) {
			ret = read_private::access((char *) p->get_data(),
					ROUND_PAGE(offset), PAGE_SIZE, access_method);
			if (ret < 0) {
				perror("read");
				exit(1);
			}
			p->set_data_ready(true);
		}
		else {
#ifdef STATISTICS
			cache_hits++;
#endif
		}

		offset -= ROUND_PAGE(offset);
		/* I assume the data I read never crosses the page boundary */
		memcpy(buf, (char *) p->get_data() + offset, size);
		ret = size;
		return ret;
	}
};

#endif
