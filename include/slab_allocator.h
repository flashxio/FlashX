#ifndef __SLAB_ALLOCATOR_H__
#define __SLAB_ALLOCATOR_H__

#include <stdlib.h>
#include <assert.h>

#include "concurrency.h"
#include "aligned_allocator.h"
#include "parameters.h"
#include "container.h"

class slab_allocator
{
public:
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
		linked_obj *end;		// point to the last element or NULL.
		int size;
	public:
		linked_obj_list() {
			size = 0;
			end = NULL;
		}

		void add(linked_obj *obj) {
			if (get_size() == 0) {
				end = obj;
				assert(head.get_next() == NULL);
			}
			head.add(obj);
			size++;
		}

		void add_list(linked_obj_list *list) {
			if (list->get_size() == 0)
				return;
			if (get_size() == 0) {
				size = list->get_size();
				head = list->head;
				end = list->end;
				list->head.set_next(NULL);
				list->end = NULL;
				list->size = 0;
			}
			else {
				list->end->set_next(head.get_next());
				head = list->head;
				size += list->size;
				list->head.set_next(NULL);
				list->end = NULL;
				list->size = 0;
			}
		}

		/**
		 * It removes the specified number of objects from the list
		 * and returns them in a form of linked list..
		 */
		linked_obj *pop(int num) {
			linked_obj *obj = head.get_next();
			int num_pop = 0;
			linked_obj *prev_obj = NULL;
			for (int i = 0; i < num && obj != NULL; i++) {
				prev_obj = obj;
				obj = obj->get_next();
				num_pop++;
			}
			linked_obj *ret = head.get_next();
			if (prev_obj)
				prev_obj->set_next(NULL);

			head.set_next(obj);
			size -= num_pop;
			assert(size >= 0);
			if (size == 0)
				end = NULL;
			return ret;
		}

		int get_size() {
			return size;
		}
	};

private:
	const int obj_size;
	// the size to increase each time there aren't enough objects
	const long increase_size;
	// the maximal size of all objects
	const long max_size;
	const int node_id;
	const int local_buf_size;

	linked_obj_list list;
	// the current size of memory used by the allocator.
	long curr_size;
	bool init;
	bool pinned;

	std::vector<char *> alloc_bufs;

	pthread_spinlock_t lock;
	// The buffers pre-allocated to serve allocation requests
	// from the local threads.
	pthread_key_t local_buf_key;
	// The buffers freed in the local threads, which hasn't been
	// added the main buffer.
	pthread_key_t local_free_key;

	std::string name;
	static atomic_integer alloc_counter;

#ifdef MEMCHECK
	aligned_allocator allocator;
#endif
public:
	slab_allocator(const std::string &name, int _obj_size, long _increase_size,
			// We allow pages to be pinned when allocated.
			long _max_size, int _node_id, bool init = false, bool pinned = false,
			int _local_buf_size = LOCAL_BUF_SIZE): obj_size(
				_obj_size), increase_size(ROUNDUP_PAGE(_increase_size)),
			max_size(_max_size), node_id(_node_id), local_buf_size(_local_buf_size)
#ifdef MEMCHECK
		, allocator(obj_size)
#endif
	{
		// To make each name unique.
		this->name = name + "-" + itoa(alloc_counter.inc(1));
		this->init = init;
		this->pinned = pinned;
		curr_size = 0;
		assert((unsigned) obj_size >= sizeof(linked_obj));
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		pthread_key_create(&local_buf_key, NULL);
		pthread_key_create(&local_free_key, NULL);
	}

	virtual ~slab_allocator();

	int get_obj_size() const {
		return obj_size;
	}

	int alloc(char **objs, int num);

	void free(char **objs, int num);

	char *alloc();

	void free(char *obj);

	// Use it carefully. It's not thread-safe.
	bool contains(const char *addr) const {
		for (unsigned i = 0; i < alloc_bufs.size(); i++) {
			if (addr >= alloc_bufs[i] && addr < alloc_bufs[i]
					+ increase_size)
				return true;
		}
		return false;
	}

	long get_max_size() const {
		return max_size;
	}
};

const static long MAX_SIZE = 0x7fffffffffffffffL;

template<class T>
class obj_initiator
{
public:
	virtual void init(T *obj) = 0;
};

template<class T>
class default_obj_initiator: public obj_initiator<T>
{
public:
	void init(T *obj) {
		T tmp;
		*obj = tmp;
	}
};

template<class T>
class obj_allocator: public slab_allocator
{
	obj_initiator<T> *initiator;
public:
	obj_allocator(const std::string &name, int node_id, long increase_size,
			long max_size = MAX_SIZE,
			obj_initiator<T> *initiator = new default_obj_initiator<T>(
				// leave some space for linked_obj, so the values in an object
				// won't be modified.
				)): slab_allocator(name, sizeof(T) + sizeof(slab_allocator::linked_obj),
				increase_size, max_size, node_id, true) {
		assert(increase_size <= max_size);
		this->initiator = initiator;
	}

	virtual int alloc_objs(T **objs, int num) {
		char *addrs[num];
		int ret = slab_allocator::alloc(addrs, num);
		for (int i = 0; i < ret; i++) {
			objs[i] = (T *) (addrs + sizeof(slab_allocator::linked_obj));
			initiator->init(objs[i]);
		}
		return ret;
	}

	virtual T *alloc_obj() {
		char *addr = slab_allocator::alloc();
		T *obj = (T *) (addr + sizeof(slab_allocator::linked_obj));
		initiator->init(obj);
		return obj;
	}

	virtual void free(T **objs, int num) {
		char *addrs[num];
		for (int i = 0; i < num; i++) {
			addrs[i] = ((char *) objs[i]) - sizeof(slab_allocator::linked_obj);
		}
		slab_allocator::free(addrs, num);
	}

	virtual void free(T *obj) {
		char *addr = ((char *) obj) - sizeof(slab_allocator::linked_obj);
		slab_allocator::free(addr);
	}
};

#endif
