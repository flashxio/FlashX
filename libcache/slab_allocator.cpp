#include <numa.h>

#include "slab_allocator.h"

int slab_allocator::alloc(char **objs, int nobjs) {
#ifdef MEMCHECK
	for (int i = 0; i < nobjs; i++)
		objs[i] = (char *) allocator.alloc(obj_size);
	return nobjs;
#else
	int num = 0;

	while (true) {
		pthread_spin_lock(&lock);
		linked_obj *o = list.pop(nobjs - num);
		pthread_spin_unlock(&lock);
		while (o != NULL) {
			objs[num++] = (char *) o;
			o = o->get_next();
		}
		if (num == nobjs)
			break;

		// This piece of code shouldn't be executed very frequently,
		// otherwise, the performance can be pretty bad.
		pthread_spin_lock(&lock);
		if (curr_size < max_size) {
			// We should increase the current size in advance, so other threads
			// can see what this thread is doing here.
			curr_size += increase_size;
			pthread_spin_unlock(&lock);
			char *objs;
			if (node_id == -1)
				objs = (char *) numa_alloc_local(increase_size);
			else
				objs = (char *) numa_alloc_onnode(increase_size, node_id);
			if (init)
				memset(objs, 0, increase_size);
			linked_obj_list tmp_list;
			for (int i = 0; i < increase_size / obj_size; i++) {
				linked_obj *header = (linked_obj *) (objs
						+ obj_size * i);
				*header = linked_obj();
				tmp_list.add(header);
			}
			pthread_spin_lock(&lock);
			alloc_bufs.push_back(objs);
			list.add_list(&tmp_list);
			pthread_spin_unlock(&lock);
		}
		else {
			pthread_spin_unlock(&lock);
			// If we can't allocate all objects, then free all objects that
			// have been allocated, and return 0.
			free(objs, num);
			return 0;
		}
	}
	return nobjs;
#endif
}

slab_allocator::~slab_allocator()
{
	pthread_spin_destroy(&lock);
}

void slab_allocator::free(char **objs, int nobjs) {
#ifdef MEMCHECK
	for (int i = 0; i < nobjs; i++)
		allocator.dealloc(objs[i]);
#else
	linked_obj_list tmp_list;
	for (int i = 0; i < nobjs; i++) {
		linked_obj *o = (linked_obj *) objs[i];
		*o = linked_obj();
		tmp_list.add(o);
	}
	pthread_spin_lock(&lock);
	list.add_list(&tmp_list);
	pthread_spin_unlock(&lock);
#endif
}

