#include "gclock.h"

void *page::data_start;

void testEnhancedGClock()
{
	const int BUF_SIZE = 5;
	enhanced_gclock_buffer buffer(BUF_SIZE);
	LF_gclock_buffer lf_buffer(BUF_SIZE);
	for (int i = 0; i < BUF_SIZE; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		frame *f1 = new frame(i * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		f->incrWC(r);
		f1->incrWC(r);
		frame *ret = buffer.add(f);
		lf_buffer.add(f1);
		assert(ret == NULL);
	}

//	buffer.print();

	for (int i = 0; i < BUF_SIZE * 4000000; i++) {
		frame *f = new frame((i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		frame *f1 = new frame((i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		f->incrWC(r);
		f1->incrWC(r);
		frame *ret = buffer.add(f);
		frame *tmp = lf_buffer.add(f1);
//		printf("the evicted frame from buffer is offset %ld, hits: %d\n",
//				ret->get_offset(), ret->getWC());
//		printf("the evicted frame from lf_buffer is offset %ld, hits: %d\n",
//				tmp->get_offset(), tmp->getWC());
		assert(ret);
		assert(ret->get_offset() == tmp->get_offset());
		delete ret;
		delete tmp;
//		buffer.print();
//		lf_buffer.print();
	}
}

void testLinkedPageList()
{
	const int LIST_SIZE = 8;
	linked_page_queue list;

	/* Add 8 frames to the list */
	for (int i = 0; i < LIST_SIZE; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		list.push_back(f);
	}
	list.print();

	frame *tmp = new frame (LIST_SIZE * PAGE_SIZE, (char *) NULL);
	printf("replace frame %ld with frame %ld\n", list.front()->get_offset(), tmp->get_offset());
	list.replace(list.front(), tmp);
	list.print();
	tmp = new frame ((LIST_SIZE + 1) * PAGE_SIZE, (char *) NULL);
	frame *replaced = list.front()->front()->front();
	printf("replace frame %ld with frame %ld\n", replaced->get_offset(), tmp->get_offset());
	list.replace(replaced, tmp);
	list.print();

	/* Randomly remove a frame from the list. */
	int r = random() % list.size();
	printf("remove %dth frame\n", r);
	list.remove(r);
	list.print();

	/* Create another list with 8 frames, and merge it to the first list. */
	linked_page_queue list1;
	for (int i = LIST_SIZE + 2; i < LIST_SIZE * 2; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		list1.push_back(f);
	}
	list1.print();
	list.merge(&list1);
	list.print();
	list1.print();
}

void testEnhanced1GClock()
{
	const int BUF_SIZE = 5;
	int ranges[] = {0};
	/*
	 * When there is only one range, enhanced_gclock_buffer1 should behave
	 * exactly the same as enhanced_gclock_buffer. 
	 */
	enhanced_gclock_buffer1 buffer(BUF_SIZE, ranges, 1);
	LF_gclock_buffer lf_buffer(BUF_SIZE);
	for (int i = 0; i < BUF_SIZE; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		frame *f1 = new frame(i * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		frame *ret = buffer.add(f);
		lf_buffer.add(f1);
		f->incrWC(r);
		f1->incrWC(r);
		assert(ret == NULL);
	}

	for (int i = 0; i < BUF_SIZE * 4000000; i++) {
		frame *f = new frame((i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		frame *f1 = new frame((i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		frame *ret = buffer.add(f);
		frame *tmp = lf_buffer.add(f1);
		f->incrWC(r);
		f1->incrWC(r);
		assert(ret);
		if (ret->get_offset() != tmp->get_offset()) {
			printf("ret: %ld, tmp: %ld\n", ret->get_offset(), tmp->get_offset());
		}
		assert(ret->get_offset() == tmp->get_offset());
		delete ret;
		delete tmp;
	}
}

int main()
{
//	testEnhancedGClock();
//	testLinkedPageList();
	testEnhanced1GClock();
}
