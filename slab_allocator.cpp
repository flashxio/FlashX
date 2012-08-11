#include <numa.h>

#include "slab_allocator.h"

int slab_allocator::alloc(char **objs, int nobjs) {
	if (list.get_size() < nobjs) {
		if (curr_size < max_size) {
			char *objs = (char *) numa_alloc_local(increase_size);
			for (int i = 0; i < increase_size / obj_size; i++) {
				linked_obj *header = (linked_obj *) (objs
						+ obj_size * i);
				*header = linked_obj();
				list.add(header);
			}
			curr_size += increase_size;
		}
		else
			return 0;
	}
	for (int i = 0; i < nobjs; i++) {
		linked_obj *o = list.pop();
		objs[i] = (char *) o;
	}
	return nobjs;
}

slab_allocator::~slab_allocator() {
}

void slab_allocator::free(char **objs, int nobjs) {
	for (int i = 0; i < nobjs; i++) {
		linked_obj *o = (linked_obj *) objs[i];
		*o = linked_obj();
		list.add(o);
	}
}

