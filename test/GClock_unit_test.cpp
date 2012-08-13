#include <sys/time.h>

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

/**
 * This is to test the linked page list. 
 * It prints the page list and I have to manually check the correctness
 * of the status of the list.
 */
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
	linked_page_queue::iterator it = list.begin();
	it.next();
	it.set(tmp);
	list.print();
	printf("\n");

	tmp = new frame ((LIST_SIZE + 1) * PAGE_SIZE, (char *) NULL);
	it.next();
	frame *replaced = it.next();
	printf("replace frame %ld with frame %ld\n", replaced->get_offset(), tmp->get_offset());
	it.set(tmp);
	list.print();
	printf("\n");

	/* Randomly remove a frame from the list. */
	int r = random() % list.size();
	printf("remove %dth frame\n", r);
	list.remove(r);
	list.print();
	printf("\n");

	/* Create another list with 8 frames, and merge it to the first list. */
	linked_page_queue list1;
	for (int i = LIST_SIZE + 2; i < LIST_SIZE * 2; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		list1.push_back(f);
	}
	printf("a new created list\n");
	list1.print();
	list.merge(&list1);
	printf("the merged list\n");
	list.print();
	printf("the new created list after merging\n");
	list1.print();
	printf("\n");

	int i = 0;
	printf("print and remove the second page simultaneously\n");
	for (linked_page_queue::iterator it = list.begin(); it.has_next(); i++) {
		it.next();
		if (2 == i)
			it.remove();
		printf("%ld\t", it.curr()->get_offset());
	}
	printf("\nthe list has %d pages\n", list.size());
}

/**
 * This is to test the behavior of the enhanced gclock with one range.
 * Since the result of the enhanced gclock with one range is exactly 
 * the same as original gclock, the test is to compare the result from
 * the two versions of gclock.
 */
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

	for (int i = 0; i < 1000000; i++) {
		frame *f = new frame(((long) i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		frame *f1 = new frame(((long) i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		frame *ret = buffer.add(f);
		frame *tmp = lf_buffer.add(f1);
		f->incrWC(r);
		f1->incrWC(r);
		assert(ret);
		printf("ret: %ld, tmp: %ld\n", ret->get_offset(), tmp->get_offset());
		assert(ret->get_offset() == tmp->get_offset());
		delete ret;
		delete tmp;
	}
}

/**
 * This is to test the behavior of the enhanced gclock with multiple ranges.
 * Its behavior is different from the original version, so this test is to run
 * many times and print the buffer each time a page is added to the buffer.
 */
void test1Enhanced1GClock()
{
	const int BUF_SIZE = 5;
	int ranges[] = {0, 5};
	enhanced_gclock_buffer1 buffer(BUF_SIZE, ranges, 2);
	for (int i = 0; i < BUF_SIZE; i++) {
		frame *f = new frame(i * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		frame *ret = buffer.add(f);
		f->incrWC(r);
		assert(ret == NULL);
	}
	
	buffer.print();

	for (int i = 0; i < 1000000; i++) {
		frame *f = new frame((long) (i + BUF_SIZE) * PAGE_SIZE, (char *) NULL);
		int r = random() % 40;
		frame *ret = buffer.add(f);
		f->incrWC(r);
		buffer.print();
		assert(ret);
		delete ret;
	}
}

/**
 * This is to measure the performance of three versions of gclock.
 * I only measure the performance in terms of the number of writes
 * to the shared memory instead of running time.
 */
void testPerformance()
{
	const int BUF_SIZE = 1000;
	int ranges[] = {0, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
	enhanced_gclock_buffer1 buffer1(BUF_SIZE, ranges, 11);
	enhanced_gclock_buffer buffer(BUF_SIZE);
	LF_gclock_buffer lf_buffer(BUF_SIZE);
	struct timeval start, end;
	srandom(0);
	gettimeofday(&start, NULL);
	for (int i = 0; i < BUF_SIZE * 1000; i++) {
		frame *f1 = new frame((long) i * PAGE_SIZE, (char *) NULL);
		int r = random() % 2048;
		f1->incrWC(r);
		frame *ret = buffer1.add(f1);
		if (ret)
			delete ret;
	}
	gettimeofday(&end, NULL);
	printf("it takes %.3f seconds\n",
			end.tv_sec - start.tv_sec + ((double) end.tv_usec - start.tv_usec) / 1000000);
	buffer1.print(true);
	printf("\n");
	srandom(0);
	gettimeofday(&start, NULL);
	for (int i = 0; i < BUF_SIZE * 1000; i++) {
		frame *f = new frame((long) i * PAGE_SIZE, (char *) NULL);
		int r = random() % 2048;
		f->incrWC(r);
		frame *ret = buffer.add(f);
		if (ret)
			delete ret;
	}
	gettimeofday(&end, NULL);
	printf("it takes %.3f seconds\n",
			end.tv_sec - start.tv_sec + ((double) end.tv_usec - start.tv_usec) / 1000000);
	buffer.print(true);
	printf("\n");
	srandom(0);
	gettimeofday(&start, NULL);
	for (int i = 0; i < BUF_SIZE * 1000; i++) {
		frame *lf_f = new frame((long) i * PAGE_SIZE, (char *) NULL);
		int r = random() % 2048;
		lf_f->incrWC(r);
		frame *ret = lf_buffer.add(lf_f);
		if (ret)
			delete ret;
	}
	gettimeofday(&end, NULL);
	printf("it takes %.3f seconds\n",
			end.tv_sec - start.tv_sec + ((double) end.tv_usec - start.tv_usec) / 1000000);
	lf_buffer.print(true);
}

int main()
{
//	testEnhancedGClock();
	printf("test the linked page list. I have to manually check its correctness\n");
	testLinkedPageList();
	printf("press any key to go to the next test\n");
	getchar();
	printf("test the enhanced gclock with one range. This test is automic.\n");
	testEnhanced1GClock();
	printf("press any key to go to the next test\n");
	getchar();
	printf("test the enhanced gclock with multiple ranges.\n");
	printf("It prints the buffer each time a page is added\n");
	test1Enhanced1GClock();
	printf("press any key to go to the next test\n");
	getchar();
	printf("measure the performance of three versions of gclock in terms of number of writes to shared memory\n");
	testPerformance();
}
