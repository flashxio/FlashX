#include <stdio.h>

#include "BeladyAlgo.h"

int main()
{
	int offs[] = {3, 2, 4, 3, 4, 2, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4, 5, 6, 7, 2, 1};
	int length = sizeof(offs) / sizeof(int);
	printf("There are %d accesses\n", length);
	belady_algo algo(3);
	indexed_offset_scanner scanner(offs, length);
	int nhits = algo.access(scanner);
	printf("There are %d hits among %d accesses\n", nhits, length);
}
