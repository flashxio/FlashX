#ifndef __PAGE_QUEUE_H__
#define __PAGE_QUEUE_H__

#include "slab_allocator.h"

class linked_page_queue;

class linked_obj
{
	linked_obj *prev, *next;
	void *payload;
public:
	linked_obj() {
		prev = next = this;
		payload = NULL;
	}

	linked_obj(void *payload) {
		prev = next = this;
		this->payload = payload;
	}

	void set_payload(void *payload) {
		this->payload = payload;
	}

	void *get_payload() const {
		return payload;
	}

	/* Add an obj behind the obj in the list. */
	void add_front(linked_obj *obj) {
		linked_obj *next = this->next;
		obj->next = next;
		obj->prev = this;
		this->next = obj;
		next->prev = obj;
	}

	/* Add an obj before the obj in the list. */
	void add_back(linked_obj *obj) {
		linked_obj *prev = this->prev;
		obj->next = this;
		obj->prev = prev;
		this->prev = obj;
		prev->next = obj;
	}

	void remove_from_list() {
		linked_obj *prev = this->prev;
		linked_obj *next = this->next;
		prev->next = next;
		next->prev = prev;
		this->next = this;
		this->prev = this;
	}

	bool is_empty() const {
		return this->next == this;
	}

	linked_obj *front() const {
		return next;
	}

	linked_obj *back() const {
		return prev;
	}

	friend class linked_page_queue;
};

/**
 * The queue is formed as a linked list, so it only supports sequential access.
 * The queue is designed to reduce writes to the shared memory when we need to
 * reorganize the page list. The idea is to split the page metadata from the 
 * data structure that forms the linked list. As long as the linked page list
 * isn't shared by multiple CPUs, any modification on the page list occurs
 * on the local memory.
 */
class linked_page_queue {
	linked_obj head;
	int _size;

	bool local_allocator;
	obj_allocator<linked_obj> *allocator;

protected:
	void remove(linked_obj *obj) {
		obj->remove_from_list();
		_size--;
		allocator->free(obj);
	}
public:
	class iterator {
		linked_obj *curr_loc;
		linked_page_queue *queue;
		int num_iter;	// number of pages that have been accessed.

		iterator(linked_obj *head, linked_page_queue *queue) {
			this->curr_loc = head;
			this->queue = queue;
			num_iter = 0;
		}
	public:
		iterator() {
			curr_loc = NULL;
			queue = NULL;
			num_iter = 0;
		}

		bool has_next() const {
			if (curr_loc == NULL || queue == NULL)
				return false;
			return num_iter < queue->size();
		}

		/* move to the next object and return the next object. */
		frame *next() {
			assert(curr_loc != NULL && queue != NULL);
			curr_loc = curr_loc->front();
			num_iter++;
			return (frame *) curr_loc->get_payload();
		}

		/* 
		 * These methods are extensions of the basic iterator.
		 * They can only be called after next() is called at least once.
		 */

		/*
		 * return the current object.
		 * return NULL if it's pointing to the head of the queue.
		 */
		frame *curr() {
			assert(curr_loc != NULL && queue != NULL);
			if (curr_loc == &queue->head)
				return NULL;
			return (frame *) curr_loc->get_payload();
		}

		void set(frame *payload) {
			assert(curr_loc != NULL && queue != NULL);
			if (curr_loc != &queue->head)
				curr_loc->set_payload(payload);
		}

		/* 
		 * remove the current frame in the queue
		 * and move to the next page automatically.
		 */
		void remove() {
			assert(curr_loc != NULL && queue != NULL);
			if (queue->size() <= 0 && curr_loc != &queue->head)
				return;
			linked_obj *tmp = curr_loc;
			curr_loc = curr_loc->back();
			num_iter--;
			queue->remove(tmp);
		}

		// for test
		linked_page_queue *owner() const {
			return queue;
		}

		friend class linked_page_queue;
	};

	virtual ~linked_page_queue() {
		if (local_allocator)
			delete allocator;
	}

	linked_page_queue(obj_allocator<linked_obj> *allocator) {
		_size = 0;
		this->allocator = allocator;
		local_allocator = false;
	}

	linked_page_queue() {
		_size = 0;
		allocator = new obj_allocator<linked_obj>(-1, PAGE_SIZE);
		local_allocator = true;
	}

	virtual iterator begin() {
		return iterator(&head, this);
	}

	virtual linked_obj *const push_back(frame *pg) {
		linked_obj *obj = allocator->alloc_obj();
		assert(obj);
		obj->set_payload(pg);
		head.add_back(obj);
		_size++;
		return obj;
	}

	virtual void pop_front() {
		if (size() == 0)
			return;
		assert(size() > 0);

		linked_obj *obj = head.front();
		obj->remove_from_list();
		_size--;
		allocator->free(obj);
	}

	bool empty() const {
		return head.is_empty();
	}

	frame *front() {
		if (size() <= 0)
			return NULL;
		return (frame *) head.front()->get_payload();
	}

	frame *back() {
		if (size() <= 0)
			return NULL;
		return (frame *) head.back()->get_payload();
	}

	int size() const {
		return _size;
	}

	void merge(linked_page_queue *list) {
		if (list->empty())
			return;

		linked_obj *list_begin = list->head.front();
		linked_obj *list_end = list->head.back();
		linked_obj *this_end = this->head.back();
		
		/* Remove all pages in the list. */
		list->head.next = &list->head;
		list->head.prev = &list->head;
		this->_size += list->_size;
		list->_size = 0;

		/* Add `list' to the end of this list. */
		this_end->next = list_begin;
		list_begin->prev = this_end;

		/* Link the end of `list' to the end of this list. */
		list_end->next = &this->head;
		this->head.prev = list_end;
	}

	// for test
	void print() const {
		linked_obj *f = head.front();
		printf("queue size: %d\n", _size);
		while (f != &head) {
			printf("%ld\t", ((frame *) f->get_payload())->get_offset());
			f = f->front();
		}
		printf("\n");
	}

	void remove(int idx) {
		int i = 0;
		for (iterator it = begin(); it.has_next(); i++) {
			it.next();
			if (idx == i)
				it.remove();
		}
	}
};

#endif
