#ifndef __SLAB_ALLOCATOR_H__
#define __SLAB_ALLOCATOR_H__

#include <stdlib.h>
#include <assert.h>

class slab_allocator
{
	class linked_obj {
		linked_obj *next;
	public:
		linked_obj() {
			next = NULL;
		}

		void add(linked_obj *obj) {
			obj->next = this->next;
			this->next = obj;
		}

		linked_obj *get_next() {
			return next;
		}

		void set_next(linked_obj *obj) {
			next = obj;
		}
	};

	class linked_obj_list {
		linked_obj head;
		int size;
	public:
		linked_obj_list() {
			size = 0;
		}

		void add(linked_obj *obj) {
			head.add(obj);
			size++;
		}

		linked_obj *pop() {
			if (size == 0)
				return NULL;
			else {
				linked_obj *ret = head.get_next();
				head.set_next(ret->get_next());
				size--;
				return ret;
			}
		}

		int get_size() {
			return size;
		}
	};

	int obj_size;
	linked_obj_list list;
	// the maximal size of all objects
	long max_size;
	// the current size of memory used by the allocator.
	long curr_size;
	// the size to increase each time there aren't enough objects
	long increase_size;
public:
	slab_allocator(int obj_size, long increase_size, long max_size) {
		this->obj_size = obj_size;
		this->max_size = max_size;
		this->curr_size = 0;
		this->increase_size = increase_size;
		assert((unsigned) obj_size >= sizeof(linked_obj));
	}

	virtual ~slab_allocator();

	int alloc(char **objs, int num);

	void free(char **objs, int num);

	long get_max_size() {
		return max_size;
	}
};

const static long MAX_SIZE = 0x7fffffffffffffffL;

template<class T>
class obj_allocator: public slab_allocator {
public:
	obj_allocator(long increase_size,
			long max_size = MAX_SIZE): slab_allocator(sizeof(T),
				increase_size, max_size) {
		assert(increase_size <= max_size);
	}

	int alloc_objs(T **objs, int num) {
		int ret = slab_allocator::alloc((char **) objs, num);
		for (int i = 0; i < ret; i++) {
			*objs[i] = T();
		}
		return ret;
	}

	T *alloc_obj() {
		T *obj = NULL;
		int ret = alloc_objs(&obj, 1);
		return ret ? obj : NULL;
	}

	void free(T **objs, int num) {
		slab_allocator::free((char **) objs, num);
	}

	void free(T *obj) {
		slab_allocator::free((char **) &obj, 1);
	}
};

#endif
