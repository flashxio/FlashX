#ifndef __MMAP_PRIVATE_H__
#define __MMAP_PRIVATE_H__

#include "thread_private.h"

extern long npages;

class mmap_private: public thread_private
{
	void *addr;
	volatile static int first[NUM_THREADS];

public:
	mmap_private(const char *new_name, int idx,
			int entry_size): thread_private(idx, entry_size) {
		static void *addr = NULL;
		static const char *file_name = NULL;
		/* if we are mapping to a different file, do the real mapping. */
		if (file_name == NULL || strcmp(file_name, new_name)) {
			int fd = open(new_name, O_RDONLY);
			int ret;

			if (fd < 0) {
				perror("open");
				exit (1);
			}
			addr = mmap(NULL, ((ssize_t) npages) * PAGE_SIZE,
					PROT_READ, MAP_PRIVATE, fd, 0);
			if (addr == NULL) {
				perror("mmap");
				exit(1);
			}
			ret = madvise(addr, ((ssize_t) npages) * PAGE_SIZE, MADV_RANDOM);
			if (ret < 0) {
				perror("madvise");
				exit(1);
			}
			file_name = new_name;
		}
		this->addr = addr;
	}

	int thread_init() {
		return 0;
	}

	ssize_t access(char *buf, off_t offset, ssize_t size, int access_method) {
		char *page = (char *) addr + offset;
		/* I try to avoid gcc optimization eliminating the code below,
		 * and it works. */
		first[idx] = page[0];
		__asm__ __volatile__("" : : "g"(&first[idx]));
		return size;
	}
};
volatile int mmap_private::first[NUM_THREADS];

#endif
