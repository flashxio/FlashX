#include <stdio.h>

#include "mem_tracker.h"

extern size_t get_alloc_objs();
extern size_t get_alloc_bytes();

struct test
{
	char data[10];
};

int main()
{
	test *t = new test;
	printf("alloc %p\n", t);
	for (int i = 0; i < sizeof(t->data); i++)
		printf("data %d:%p: %d\n", i, &t->data[i], t->data[i]);
	printf("alloc %ld objs and %ld bytes\n",
			get_alloc_objs(), get_alloc_bytes());
	delete t;
	printf("alloc %ld objs and %ld bytes\n",
			get_alloc_objs(), get_alloc_bytes());

	t = new test[3];
	printf("alloc %p\n", t);
	for (int j = 0; j < 3; j++)
		for (int i = 0; i < sizeof(t->data); i++)
			printf("data %d:%d:%p: %d\n", j, i, &t[j].data[i], t[j].data[i]);
	printf("alloc %ld objs and %ld bytes\n",
			get_alloc_objs(), get_alloc_bytes());
	delete [] t;
	printf("alloc %ld objs and %ld bytes\n",
			get_alloc_objs(), get_alloc_bytes());
}
