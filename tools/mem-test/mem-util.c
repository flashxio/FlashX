#define ENTRY_SIZE 32

long sum1(char *array, int start_i, int end_i)
{
	int i;
	long res = 0;

	for (i = start_i; i < end_i; i++) {
		res += *((long *) (array + i * ENTRY_SIZE));
	}
	return res;
}

long sum2(char *array, int start_i, int end_i)
{
	int i;
	long res = 0;

	for (i = start_i; i < end_i; i++) {
		res += *((long *) (array + i * ENTRY_SIZE));
	}
	return res;
}
