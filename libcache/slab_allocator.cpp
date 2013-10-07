#include <numa.h>
#include <sys/mman.h>

#include "slab_allocator.h"

void slab_allocator::free(char *obj)
{
	if (local_buf_size == 0) {
		slab_allocator::free(&obj, 1);
	}
	else {
		fifo_queue<char *> *local_free_refs
			= (fifo_queue<char *> *) pthread_getspecific(local_free_key);
		if (local_free_refs == NULL) {
			local_free_refs = fifo_queue<char *>::create(node_id, local_buf_size);
			pthread_setspecific(local_free_key, local_free_refs);
		}
		if (local_free_refs->is_full()) {
			char *objs[local_buf_size];
			int num = local_free_refs->fetch(objs, local_buf_size);
			slab_allocator::free(objs, num);
		}
		local_free_refs->push_back(obj);
	}
}

char *slab_allocator::alloc()
{
	if (local_buf_size == 0) {
		char *obj;
		int num = alloc(&obj, 1);
		if (num == 0)
			return NULL;
		else
			return obj;
	}
	else {
		fifo_queue<char *> *local_buf_refs
			= (fifo_queue<char *> *) pthread_getspecific(local_buf_key);
		if (local_buf_refs == NULL) {
			assert(node_id >= 0);
			local_buf_refs = fifo_queue<char *>::create(node_id, local_buf_size);
			pthread_setspecific(local_buf_key, local_buf_refs);
		}

		if (local_buf_refs->is_empty()) {
			char *objs[local_buf_size];
			int num = alloc(objs, local_buf_size);
			if (num == 0)
				return NULL;
			int num_added = local_buf_refs->add(objs, num);
			assert(num_added == num);
		}
		return local_buf_refs->pop_front();
	}
}

int slab_allocator::alloc(char **objs, int nobjs) {
#ifdef MEMCHECK
	for (int i = 0; i < nobjs; i++) {
		objs[i] = (char *) allocator.alloc(obj_size);
		if (init)
			memset(objs[i], 0, obj_size);
	}
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
#ifdef USE_IOAT
			if (pinned) {
				int ret = mlock(objs, increase_size);
				if (ret < 0)
					perror("mlock");
				assert(ret == 0);
			}
#endif
			assert(((long) objs) % PAGE_SIZE == 0);
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
	for (unsigned i = 0; i < alloc_bufs.size(); i++) {
#ifdef USE_IOAT
		if (pinned) {
#ifdef DEBUG
			printf("unpin buf %p of %ld bytes\n", alloc_bufs[i], increase_size);
#endif
			munlock(alloc_bufs[i], increase_size);
		}
#endif
		numa_free(alloc_bufs[i], increase_size);
	}
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

