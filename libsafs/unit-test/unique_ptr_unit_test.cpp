#include <memory>

#include "slab_allocator.h"

class slab_deleter
{
	slab_allocator *alloc;
public:
	slab_deleter(slab_allocator *alloc) {
		this->alloc = alloc;
	}

	void operator()(char *ptr) {
		alloc->free(ptr);
	}
};

/**
 * This is to test how unique_ptr works with a slab allocator.
 */
int main()
{

	slab_allocator *alloc = new slab_allocator("test", 8, 1024, 1024, 0);
	std::unique_ptr<char, slab_deleter> req(alloc->alloc(), slab_deleter(alloc));
	printf("size of unique_ptr: %ld\n", sizeof(req));
}
