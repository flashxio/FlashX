#include "gclock.h"

void testEnhancedGClock()
{
	const int BUF_SIZE = 5;
	enhanced_gclock_buffer buffer(BUF_SIZE);
	LF_gclock_buffer lf_buffer(BUF_SIZE);
	for (int i = 0; i < BUF_SIZE; i++) {
		frame *f = new frame(i * PAGE_SIZE, NULL);
		frame *f1 = new frame(i * PAGE_SIZE, NULL);
		int r = random() % 40;
		f->incrWC(r);
		f1->incrWC(r);
		frame *ret = buffer.add(f);
		lf_buffer.add(f1);
		assert(ret == NULL);
	}

//	buffer.print();

	for (int i = 0; i < BUF_SIZE * 4000000; i++) {
		frame *f = new frame((i + BUF_SIZE) * PAGE_SIZE, NULL);
		frame *f1 = new frame((i + BUF_SIZE) * PAGE_SIZE, NULL);
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

int main()
{
	testEnhancedGClock();
}
