#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#include <string>
#include <vector>

#include "crs_header.h"

#define MKL_ILP64
#define FLEXCOMPLEX
#define MKL_INT crs_idx_t
#include <mkl.h>

inline static float time_diff(struct timeval time1, struct timeval time2)
{
	return time2.tv_sec - time1.tv_sec
			+ ((float)(time2.tv_usec - time1.tv_usec))/1000000;
}

void test_spmv(const std::vector<crs_idx_t> &row_idxs,
		const std::vector<crs_idx_t> &col_vec, const std::vector<double> &val_vec,
		size_t num_rows, size_t num_cols)
{
	printf("test SpMV\n");
	std::vector<double> in_vec(num_cols);
	printf("init input vector of %ld entries\n", in_vec.size());
#pragma omp parallel for
	for (size_t i = 0; i < in_vec.size(); i++)
		in_vec[i] = 1;
	std::vector<double> out_vec(num_rows);
	printf("init output vector of %ld entries\n", out_vec.size());

	crs_idx_t nrowA = num_rows;
	assert(col_vec.size() == val_vec.size());
	assert(row_idxs.back() == val_vec.size());
	struct timeval start, end;
	gettimeofday(&start, NULL);
	assert(nrowA == row_idxs.size() - 1);
	mkl_cspblas_dcsrgemv("N", &nrowA, val_vec.data(), row_idxs.data(), col_vec.data(),
			in_vec.data(), out_vec.data());
	gettimeofday(&end, NULL);
	printf("SpMV takes %.3fs\n", time_diff(start, end));

	assert(out_vec.size() == row_idxs.size() - 1);
	for (size_t i = 0; i < out_vec.size(); i++) {
		assert(out_vec[i] == row_idxs[i + 1] - row_idxs[i]);
	}
}

void test_spmm(const std::vector<crs_idx_t> &row_idxs,
		const std::vector<crs_idx_t> &col_vec, const std::vector<double> &val_vec,
		size_t num_rows, size_t num_cols, int num_vecs)
{
	printf("test SpMM on %d vectors\n", num_vecs);
	std::vector<double> in_vec(num_cols * num_vecs);
	printf("init input vector of %ld entries\n", in_vec.size());
	// The input matrix is organized in row major.
#pragma omp parallel for
	for (size_t i = 0; i < in_vec.size(); i++) {
		size_t col_idx = i % num_vecs;
		in_vec[i] = col_idx + 1;
	}
	std::vector<double> out_vec(num_rows * num_vecs);
	printf("init output vector of %ld entries\n", out_vec.size());

	assert(col_vec.size() == val_vec.size());
	assert(row_idxs.back() == val_vec.size());
	struct timeval start, end;
	gettimeofday(&start, NULL);
	crs_idx_t nrowA = num_rows;
	crs_idx_t ncolC = num_vecs;
	crs_idx_t ncolA = num_cols;
	double alpha = 1;
	double beta = 0;
	mkl_dcsrmm("N", &nrowA, &ncolC, &ncolA, &alpha, "G  C", val_vec.data(),
			col_vec.data(), row_idxs.data(), row_idxs.data() + 1, in_vec.data(),
			&ncolC, &beta, out_vec.data(), &ncolC);
	gettimeofday(&end, NULL);
	printf("SpMM takes %.3fs\n", time_diff(start, end));

	assert(out_vec.size() == (row_idxs.size() - 1) * num_vecs);
	// The output matrix is organized in row major.
#pragma omp parallel for
	for (size_t i = 0; i < out_vec.size(); i++) {
		size_t col_idx = i % num_vecs;
		size_t row_idx = i / num_vecs;
		assert(out_vec[i] == (row_idxs[row_idx + 1] - row_idxs[row_idx]) * (col_idx + 1));
	}
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test-mkl_multiply crs_file\n");
		exit(1);
	}

	std::string crs_file = argv[1];
	FILE *f = fopen(crs_file.c_str(), "r");
	if (f == NULL) {
		fprintf(stderr, "can't open %s: %s\n", crs_file.c_str(), strerror(errno));
		exit(1);
	}

	crs_header header;
	size_t ret = fread(&header, sizeof(header), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read header: %s\n", strerror(errno));
		exit(1);
	}

	size_t num_rows = header.get_num_rows();
	size_t num_cols = header.get_num_cols();
	size_t nnz = header.get_num_non_zeros();
	printf("There are %ld rows, %ld cols and %ld nnz\n", num_rows,
			num_cols, nnz);

	std::vector<crs_idx_t> row_idxs(num_rows + 1);
	std::vector<crs_idx_t> col_vec(nnz);
	printf("read CRS\n");
	ret = fread(row_idxs.data(), row_idxs.size() * sizeof(row_idxs[0]), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read row idxs: %s\n", strerror(errno));
		exit(1);
	}

	ret = fread(col_vec.data(), col_vec.size() * sizeof(col_vec[0]), 1, f);
	if (ret == 0) {
		fprintf(stderr, "can't read col idxs: %s\n", strerror(errno));
		exit(1);
	}

	for (size_t i = 0; i < col_vec.size(); i++)
		assert(col_vec[i] < num_cols);

	std::vector<double> val_vec(nnz);
	printf("Create value vector of %ld entries for the adj matrix\n", val_vec.size());
#pragma omp parallel for
	for (size_t i = 0; i < val_vec.size(); i++)
		val_vec[i] = 1;

	for (size_t i = 0; i < 5; i++)
		test_spmv(row_idxs, col_vec, val_vec, num_rows, num_cols);
	for (int num_vecs = 2; num_vecs <= 128; num_vecs *= 2) {
		for (size_t i = 0; i < 5; i++)
			test_spmm(row_idxs, col_vec, val_vec, num_rows, num_cols, num_vecs);
	}
}
