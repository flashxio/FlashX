#include <cblas.h>
#include <set>

#include "bulk_operate.h"
#include "mem_matrix_store.h"
#include "mem_vec_store.h"

using namespace fm;
using namespace fm::detail;

void test_reset1(std::shared_ptr<mem_matrix_store> store)
{
	store->reset_data();
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == 0);
}

void test_reset(size_t long_dim)
{
	printf("test reset on mem matrix, long dim: %ld\n", long_dim);
	// Test on mem matrix.
	test_reset1(mem_col_matrix_store::create(long_dim, 10,
				get_scalar_type<int>()));
	test_reset1(mem_row_matrix_store::create(long_dim, 10,
				get_scalar_type<int>()));

	// Test on mem sub matrix.
	mem_col_matrix_store::ptr col_store = mem_col_matrix_store::create(
				long_dim, 10, get_scalar_type<int>());
	mem_row_matrix_store::ptr row_store = mem_row_matrix_store::create(
				long_dim, 10, get_scalar_type<int>());
	std::vector<off_t> cols;
	for (size_t i = 0; i < cols.size(); i++) {
		if (i % 2 == 0)
			cols.push_back(i);
	}
	std::vector<off_t> rows;
	for (size_t i = 0; i < rows.size(); i++) {
		if (i % 2 == 0)
			rows.push_back(i);
	}
	test_reset1(mem_sub_col_matrix_store::create(*col_store, cols));
	test_reset1(mem_sub_row_matrix_store::create(*row_store, rows));
}

class set_col_operate: public type_set_operate<int>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

class set_row_operate: public type_set_operate<int>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

void verify_set(std::shared_ptr<mem_matrix_store> store)
{
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);
}

void test_set1(std::shared_ptr<mem_matrix_store> store)
{
	if (store->store_layout() == matrix_layout_t::L_COL)
		store->set_data(set_col_operate(store->get_num_cols()));
	else
		store->set_data(set_row_operate(store->get_num_cols()));
	verify_set(store);
}

void test_set(size_t long_dim)
{
	printf("test set on mem matrix, long dim: %ld\n", long_dim);
	// Test set on mem matrix.
	test_set1(mem_col_matrix_store::create(long_dim, 10, get_scalar_type<int>()));
	test_set1(mem_row_matrix_store::create(long_dim, 10, get_scalar_type<int>()));

	// Test set on sub matrix.
	mem_col_matrix_store::ptr col_store = mem_col_matrix_store::create(
				long_dim, 10, get_scalar_type<int>());
	mem_row_matrix_store::ptr row_store = mem_row_matrix_store::create(
				long_dim, 10, get_scalar_type<int>());
	std::vector<off_t> cols;
	for (size_t i = 0; i < cols.size(); i++) {
		if (i % 2 == 0)
			cols.push_back(i);
	}
	std::vector<off_t> rows;
	for (size_t i = 0; i < rows.size(); i++) {
		if (i % 2 == 0)
			rows.push_back(i);
	}
	test_set1(mem_sub_col_matrix_store::create(*col_store, cols));
	test_set1(mem_sub_row_matrix_store::create(*row_store, rows));
}

void test_sub_col_matrix()
{
	printf("test submatrix of a column-wise matrix\n");
	mem_col_matrix_store::ptr store = mem_col_matrix_store::create(
				1000, 10, get_scalar_type<int>());
	store->set_data(set_col_operate(store->get_num_cols()));
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_col_matrix_store::ptr sub_store = mem_sub_col_matrix_store::create(
			*store, idxs);
	assert(sub_store != NULL);
	assert(sub_store->get_num_rows() == store->get_num_rows());
	assert(sub_store->get_num_cols() == idxs.size());
	assert(sub_store->store_layout() == store->store_layout());
	assert(sub_store->get_entry_size() == store->get_entry_size());
	assert(sub_store->get_type() == store->get_type());

	std::vector<off_t> idxs2(2);
	idxs2[0] = 0;
	idxs2[1] = 1;
	mem_col_matrix_store::const_ptr subsub_store = mem_col_matrix_store::cast(
			sub_store->get_cols(idxs2));
	assert(subsub_store != NULL);
	assert(subsub_store->get_num_cols() == idxs2.size());
}

void test_sub_row_matrix()
{
	printf("test submatrix of a row-wise matrix\n");
	mem_row_matrix_store::ptr store = mem_row_matrix_store::create(
				1000, 10, get_scalar_type<int>());
	store->set_data(set_row_operate(store->get_num_cols()));
	std::set<off_t> idxs_set;
	for (size_t i = 0; i < store->get_num_rows() / 2; i++)
		idxs_set.insert(random() % store->get_num_rows());
	std::vector<off_t> idxs(idxs_set.begin(), idxs_set.end());
	mem_row_matrix_store::ptr sub_store = mem_sub_row_matrix_store::create(
			*store, idxs);
	assert(sub_store != NULL);
	assert(sub_store->get_num_rows() == idxs.size());
	assert(sub_store->get_num_cols() == store->get_num_cols());
	assert(sub_store->store_layout() == store->store_layout());
	assert(sub_store->get_entry_size() == store->get_entry_size());
	assert(sub_store->get_type() == store->get_type());

	idxs_set.clear();
	for (size_t i = 0; i < sub_store->get_num_rows() / 2; i++)
		idxs_set.insert(random() % sub_store->get_num_rows());
	std::vector<off_t> idxs2(idxs_set.begin(), idxs_set.end());
	mem_row_matrix_store::const_ptr subsub_store = mem_row_matrix_store::cast(
			sub_store->get_rows(idxs2));
	assert(subsub_store != NULL);
	assert(subsub_store->get_num_rows() == idxs2.size());
}

void test_io()
{
	printf("test read/write matrix to a file\n");
	std::string out_file;
	mem_matrix_store::const_ptr read_mat;
	mem_matrix_store::ptr orig_mat;

	orig_mat = mem_row_matrix_store::create(
				999999, 10, get_scalar_type<int>());
	orig_mat->set_data(set_row_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file, -1);
	assert(read_mat);
	assert(read_mat->store_layout() == matrix_layout_t::L_ROW);
	for (size_t i = 0; i < orig_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < orig_mat->get_num_cols(); j++)
			assert(orig_mat->get<int>(i, j) == read_mat->get<int>(i, j));
	}
	unlink(out_file.c_str());

	orig_mat = mem_col_matrix_store::create(
				999999, 10, get_scalar_type<int>());
	orig_mat->set_data(set_col_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file, -1);
	assert(read_mat);
	assert(read_mat->store_layout() == matrix_layout_t::L_COL);
	for (size_t i = 0; i < orig_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < orig_mat->get_num_cols(); j++)
			assert(orig_mat->get<int>(i, j) == read_mat->get<int>(i, j));
	}
	unlink(out_file.c_str());

	orig_mat = mem_row_matrix_store::create(
				10, 999999, get_scalar_type<int>());
	orig_mat->set_data(set_row_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file, -1);
	assert(read_mat);
	assert(read_mat->store_layout() == matrix_layout_t::L_ROW);
	for (size_t i = 0; i < orig_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < orig_mat->get_num_cols(); j++)
			assert(orig_mat->get<int>(i, j) == read_mat->get<int>(i, j));
	}
	unlink(out_file.c_str());

	orig_mat = mem_col_matrix_store::create(
				10, 999999, get_scalar_type<int>());
	orig_mat->set_data(set_col_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file, -1);
	assert(read_mat);
	assert(read_mat->store_layout() == matrix_layout_t::L_COL);
	for (size_t i = 0; i < orig_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < orig_mat->get_num_cols(); j++)
			assert(orig_mat->get<int>(i, j) == read_mat->get<int>(i, j));
	}
	unlink(out_file.c_str());
}

void test_transpose()
{
	printf("test transpose\n");
	mem_col_matrix_store::ptr store = mem_col_matrix_store::create(
				1000, 10, get_scalar_type<int>());
	store->set_data(set_col_operate(store->get_num_cols()));
	std::vector<off_t> idxs(3);
	idxs[0] = 1;
	idxs[1] = 5;
	idxs[2] = 3;
	mem_col_matrix_store::ptr sub_store = mem_sub_col_matrix_store::create(
			*store, idxs);

	mem_matrix_store::const_ptr t_store = mem_matrix_store::cast(
			store->transpose());
	assert(t_store->get_num_rows() == store->get_num_cols());
	assert(t_store->get_num_cols() == store->get_num_rows());
	for (size_t i = 0; i < t_store->get_num_rows(); i++)
		for (size_t j = 0; j < t_store->get_num_cols(); j++)
			assert(t_store->get<int>(i, j) == j * t_store->get_num_rows() + i);

	mem_matrix_store::const_ptr t_sub_store = mem_matrix_store::cast(
			sub_store->transpose());
	assert(t_sub_store->get_num_rows() == sub_store->get_num_cols());
	assert(t_sub_store->get_num_cols() == sub_store->get_num_rows());
	for (size_t i = 0; i < t_sub_store->get_num_rows(); i++)
		for (size_t j = 0; j < t_sub_store->get_num_cols(); j++)
			assert(t_sub_store->get<int>(i, j)
					== j * t_store->get_num_rows() + idxs[i]);
}

bool is_symmetric(mem_matrix_store::const_ptr store)
{
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = i + 1; j < store->get_num_cols(); j++)
			if (memcmp(store->get(i, j), store->get(j, i),
						store->get_entry_size()) != 0)
				return false;
	return true;
}

void test_symmetrize()
{
	printf("test symmetrize\n");
	bool ret;
	mem_matrix_store::ptr mat = mem_matrix_store::create(1000, 1000,
			matrix_layout_t::L_ROW, get_scalar_type<size_t>(), -1);
	scalar_variable_impl<size_t> start(1);
	scalar_variable_impl<size_t> stride(1);
	auto set = get_scalar_type<size_t>().get_set_seq(start, stride,
			mat->get_num_rows(), mat->get_num_cols(), true,
			mat->store_layout());
	mat->set_data(*set);
	ret = mat->symmetrize(true);
	assert(ret);
	assert(is_symmetric(mat));

	ret = mat->symmetrize(false);
	assert(ret);
	assert(is_symmetric(mat));

	mat = mem_matrix_store::create(1000, 1000, matrix_layout_t::L_COL,
			get_scalar_type<size_t>(), -1);
	mat->set_data(*set);
	ret = mat->symmetrize(true);
	assert(ret);
	assert(is_symmetric(mat));

	ret = mat->symmetrize(false);
	assert(ret);
	assert(is_symmetric(mat));
}

void wide_dsyrk_row(const std::pair<size_t, size_t> &Asize,
		const double *Amat, double *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.second;
	size_t ldc = Asize.first;
	cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
}

void wide_dsyrk_col(const std::pair<size_t, size_t> &Asize,
		const double *Amat, double *res_mat)
{
	size_t n = Asize.first;
	size_t k = Asize.second;
	size_t lda = Asize.first;
	size_t ldc = Asize.first;
	cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, n, k, 1, Amat, lda,
			0, res_mat, ldc);
}

void wide_dgemm_row(const std::pair<size_t, size_t> &Asize,
		const double *Amat, const std::pair<size_t, size_t> &Bsize,
		const double *Bmat, double *res_mat, size_t out_num_cols)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.second,
			Bmat, Bsize.second, 0, res_mat, out_num_cols);
}

void wide_dgemm_col(const std::pair<size_t, size_t> &Asize,
		const double *Amat, const std::pair<size_t, size_t> &Bsize,
		const double *Bmat, double *res_mat, size_t out_num_rows)
{
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
			Asize.first, Bsize.second, Asize.second, 1, Amat, Asize.first,
			Bmat, Bsize.first, 0, res_mat, out_num_rows);
}

bool same_data(mem_matrix_store::const_ptr mat1, mem_matrix_store::const_ptr mat2)
{
	for (size_t i = 0; i < mat1->get_num_rows(); i++)
		for (size_t j = 0; j < mat1->get_num_cols(); j++)
			if (memcmp(mat1->get(i, j), mat2->get(i, j), mat1->get_entry_size()))
				return false;
	return true;
}

void test_multiply(matrix_layout_t layout)
{
	printf("test multiply\n");
	mem_matrix_store::ptr mat1 = mem_matrix_store::create(10, 1000, layout,
			get_scalar_type<double>(), -1);
	scalar_variable_impl<double> start(1);
	scalar_variable_impl<double> stride(1);
	auto set = get_scalar_type<double>().get_set_seq(start, stride,
			mat1->get_num_rows(), mat1->get_num_cols(), true,
			mat1->store_layout());
	mat1->set_data(*set);
	mem_matrix_store::ptr mat2 = mem_matrix_store::create(1000, 10, layout,
			get_scalar_type<double>(), -1);
	set = get_scalar_type<double>().get_set_seq(start, stride,
			mat2->get_num_rows(), mat2->get_num_cols(), false,
			mat2->store_layout());
	mat2->set_data(*set);
	assert(same_data(mat1, std::dynamic_pointer_cast<const mem_matrix_store>(
					mat2->transpose())));
	mem_matrix_store::ptr res1 = mem_matrix_store::create(10, 10, layout,
			get_scalar_type<double>(), -1);
	res1->reset_data();
	mem_matrix_store::ptr res2 = mem_matrix_store::create(10, 10, layout,
			get_scalar_type<double>(), -1);
	res2->reset_data();

	std::pair<size_t, size_t> Asize(mat1->get_num_rows(), mat1->get_num_cols());
	std::pair<size_t, size_t> Bsize(mat2->get_num_rows(), mat2->get_num_cols());
	size_t out_size
		= res1->get_num_rows() * res1->get_num_cols() * res1->get_entry_size();
	if (layout == matrix_layout_t::L_COL) {
		wide_dgemm_col(Asize, (double *) mat1->get_raw_arr(), Bsize,
				(double *) mat2->get_raw_arr(), (double *) res1->get_raw_arr(),
				res1->get_num_rows());
		wide_dsyrk_col(Asize, (double *) mat1->get_raw_arr(),
				(double *) res2->get_raw_arr());
		res2->symmetrize(true);
		assert(memcmp(res1->get_raw_arr(), res2->get_raw_arr(), out_size) == 0);
	}
	else {
		wide_dgemm_row(Asize, (double *) mat1->get_raw_arr(), Bsize,
				(double *) mat2->get_raw_arr(), (double *) res1->get_raw_arr(),
				res1->get_num_cols());
		wide_dsyrk_row(Asize, (double *) mat1->get_raw_arr(),
				(double *) res2->get_raw_arr());
		res2->symmetrize(true);
		assert(memcmp(res1->get_raw_arr(), res2->get_raw_arr(), out_size) == 0);
	}
}

void test_resize(mem_matrix_store::ptr mat1)
{
	mat1->set_data(set_col_operate(mat1->get_num_cols()));
	mem_matrix_store::ptr mat2 = mem_matrix_store::create(mat1->get_num_rows(),
			mat1->get_num_cols(), mat1->store_layout(), mat1->get_type(), -1);
	mat2->set_data(set_col_operate(mat1->get_num_cols()));
	mat2->resize(random() % mat2->get_num_rows(), mat2->get_num_cols());
	for (size_t i = 0; i < mat2->get_num_rows(); i++)
		for (size_t j = 0; j < mat2->get_num_cols(); j++)
			assert(mat1->get<int>(i, j) == mat2->get<int>(i, j));
}

void test_resize()
{
	printf("test resize\n");
	mem_matrix_store::ptr mat1 = mem_matrix_store::create(100000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<int>(), -1);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(100000, 10, matrix_layout_t::L_COL,
			get_scalar_type<int>(), -1);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(10, 100000, matrix_layout_t::L_ROW,
			get_scalar_type<int>(), -1);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(10, 100000, matrix_layout_t::L_COL,
			get_scalar_type<int>(), -1);
	test_resize(mat1);
}

int main()
{
	test_resize();
	test_multiply(matrix_layout_t::L_COL);
	test_multiply(matrix_layout_t::L_ROW);
	test_symmetrize();
	test_reset(1000);
	test_reset(1000000);
	test_set(1000);
	test_set(1000000);
	test_sub_col_matrix();
	test_sub_row_matrix();
	test_io();
	test_transpose();
}
