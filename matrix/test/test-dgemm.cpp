#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <malloc.h>

#ifdef USE_MKL
#define MKL_ILP64
#define FLEXCOMPLEX
#define MKL_INT size_t
#include <mkl.h>
#else
#include <cblas.h>
#endif

inline static float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

void test_dgemm(size_t m, size_t n, size_t k)
{
	struct timeval start, end;
	double *A, *B, *C;

	printf("test tall col-wise matrix in MKL BLAS: M(%ld x %ld) * M(%ld %ld)\n",
			m, k, k, n);
	A = (double *) memalign(64, m*k*sizeof( double ));
	B = (double *) memalign(64, k*n*sizeof( double ));
	if (A == NULL || B == NULL) {
		printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
		free(A);
		free(B);
		return;
	}
#pragma omp parallel for
	for (size_t i = 0; i < (m*k); i++) {
		A[i] = (double)(i+1);
	}

	for (size_t i = 0; i < (k*n); i++) {
		B[i] = (double)(-i-1);
	}

//#pragma omp parallel for
//	for (size_t i = 0; i < (m*n); i++) {
//		C[i] = 0.0;
//	}

	gettimeofday(&start, NULL);
	C = (double *) memalign(64, m*n*sizeof( double ));
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A,
			m, B, k, 0, C, m);
	gettimeofday(&end, NULL);
	printf("It takes %.3f seconds to multiply row matrix in parallel\n",
			time_diff(start, end));
}

int main()
{
	test_dgemm(1024 * 1024 * 100, 20, 20);
}
