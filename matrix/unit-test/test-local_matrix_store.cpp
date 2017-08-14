#include <cblas.h>

#include "bulk_operate.h"
#include "local_matrix_store.h"
#include "dense_matrix.h"
#include "col_vec.h"
#include "factor.h"
#include "mapply_matrix_store.h"

using namespace fm;
using namespace fm::detail;

class ltest_col_matrix_store: public lvirtual_col_matrix_store
{
	local_col_matrix_store::ptr buf;
public:
	ltest_col_matrix_store(local_col_matrix_store::ptr buf): lvirtual_col_matrix_store(
			buf->get_global_start_row(), buf->get_global_start_col(),
			buf->get_num_rows(), buf->get_num_cols(), buf->get_type(),
			buf->get_node_id()) {
		this->buf = buf;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		buf->resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		buf->reset_size();
	}

	using lvirtual_col_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		return buf->get_raw_arr();
	}

	using lvirtual_col_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		return matrix_store::const_ptr();
	}

	using lvirtual_col_matrix_store::get_col;
	virtual const char *get_col(size_t col) const {
		return buf->get_col(col);
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		return buf->get_portion(local_start_row, local_start_col,
				num_rows, num_cols);
	}
};

class ltest_row_matrix_store: public lvirtual_row_matrix_store
{
	local_row_matrix_store::ptr buf;
public:
	ltest_row_matrix_store(local_row_matrix_store::ptr buf): lvirtual_row_matrix_store(
			buf->get_global_start_row(), buf->get_global_start_col(),
			buf->get_num_rows(), buf->get_num_cols(), buf->get_type(),
			buf->get_node_id()) {
		this->buf = buf;
	}

	virtual bool resize(off_t local_start_row, off_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) {
		buf->resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
		return local_matrix_store::resize(local_start_row, local_start_col,
				local_num_rows, local_num_cols);
	}
	virtual void reset_size() {
		buf->reset_size();
	}

	using lvirtual_row_matrix_store::get_raw_arr;
	virtual const char *get_raw_arr() const {
		return buf->get_raw_arr();
	}

	using lvirtual_row_matrix_store::transpose;
	virtual matrix_store::const_ptr transpose() const {
		return matrix_store::const_ptr();
	}

	using lvirtual_row_matrix_store::get_row;
	virtual const char *get_row(size_t row) const {
		return buf->get_row(row);
	}

	virtual local_matrix_store::const_ptr get_portion(
			size_t local_start_row, size_t local_start_col, size_t num_rows,
			size_t num_cols) const {
		return buf->get_portion(local_start_row, local_start_col,
				num_rows, num_cols);
	}
};

void test_reset1(std::shared_ptr<local_matrix_store> store)
{
	assert(!store->read_only());
	store->reset_data();
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == 0);
}

void test_reset(size_t long_dim)
{
	printf("test reset on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_reset1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_reset1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_reset1(std::shared_ptr<local_matrix_store>(
				new local_ref_contig_col_matrix_store(col_store->get_raw_arr(), 0, 0,
					long_dim, 10, get_scalar_type<int>(), -1)));
	test_reset1(std::shared_ptr<local_matrix_store>(
				new local_ref_contig_row_matrix_store(row_store->get_raw_arr(), 0, 0,
					long_dim, 10, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	col_store = std::shared_ptr<local_col_matrix_store>(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	row_store = std::shared_ptr<local_row_matrix_store>(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_reset1(std::shared_ptr<local_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_reset1(std::shared_ptr<local_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0, rows.size(),
					row_store->get_num_cols(), get_scalar_type<int>(), -1)));
}

template<class T>
class set_col_operate: public type_set_operate<T>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(T *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

template<class T>
class set_row_operate: public type_set_operate<T>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(T *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

class set1_col_operate: public type_set_operate<int>
{
public:
	set1_col_operate() { }

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1;
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

class set1_row_operate: public type_set_operate<int>
{
public:
	set1_row_operate() { }

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1;
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

void verify_set(std::shared_ptr<local_matrix_store> store)
{
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	if (store->is_whole()) {
		size_t orig_num_rows = store->get_num_rows();
		size_t orig_num_cols = store->get_num_cols();
		size_t start_row = 2;
		size_t start_col = 0;

		store->resize(start_row, start_col, store->get_num_rows() - start_row,
				store->get_num_cols() - start_col);
		for (size_t i = 0; i < store->get_num_rows(); i++)
			for (size_t j = 0; j < store->get_num_cols(); j++)
				assert(store->get<int>(i, j)
						== (i + start_row) * orig_num_cols + (j + start_col));

		store->reset_size();
		start_row = 0;
		start_col = 2;
		store->resize(start_row, start_col, store->get_num_rows() - start_row,
				store->get_num_cols() - start_col);
		for (size_t i = 0; i < store->get_num_rows(); i++)
			for (size_t j = 0; j < store->get_num_cols(); j++)
				assert(store->get<int>(i, j)
						== (i + start_row) * orig_num_cols + (j + start_col));

		store->reset_size();
		start_row = 2;
		start_col = 2;
		store->resize(start_row, start_col, store->get_num_rows() - start_row,
				store->get_num_cols() - start_col);
		for (size_t i = 0; i < store->get_num_rows(); i++)
			for (size_t j = 0; j < store->get_num_cols(); j++)
				assert(store->get<int>(i, j)
						== (i + start_row) * orig_num_cols + (j + start_col));
	}
}

void test_set1(std::shared_ptr<local_matrix_store> store)
{
	assert(store->is_whole());
	assert(!store->read_only());
	if (store->store_layout() == matrix_layout_t::L_COL)
		store->set_data(set_col_operate<int>(store->get_num_cols()));
	else
		store->set_data(set_row_operate<int>(store->get_num_cols()));
	verify_set(store);
}

void test_set(size_t long_dim)
{
	printf("test set on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_set1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_set1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_set1(std::shared_ptr<local_matrix_store>(
				new local_ref_contig_col_matrix_store(col_store->get_raw_arr(), 0, 0,
					long_dim, 10, get_scalar_type<int>(), -1)));
	test_set1(std::shared_ptr<local_matrix_store>(
				new local_ref_contig_row_matrix_store(row_store->get_raw_arr(), 0, 0,
					long_dim, 10, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	col_store = std::shared_ptr<local_col_matrix_store>(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	row_store = std::shared_ptr<local_row_matrix_store>(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_set1(std::shared_ptr<local_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_set1(std::shared_ptr<local_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0, rows.size(),
					row_store->get_num_cols(), get_scalar_type<int>(), -1)));

	std::shared_ptr<local_matrix_store> tmp;
	// Test on local const reference matrix to a matrix stored contiguously.
	tmp = std::shared_ptr<local_matrix_store>(
			new local_cref_contig_col_matrix_store(col_store->get_raw_arr(), 0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	verify_set(tmp);
	tmp = std::shared_ptr<local_matrix_store>(
			new local_cref_contig_row_matrix_store(row_store->get_raw_arr(), 0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	verify_set(tmp);

	// Test on local const reference matrix to a matrix stored non-contiguously.
	std::vector<const char *> const_cols(col_store->get_num_cols());
	for (size_t i = 0; i < const_cols.size(); i++)
		const_cols[i] = col_store->get_col(i);
	std::vector<const char *> const_rows(row_store->get_num_rows());
	for (size_t i = 0; i < const_rows.size(); i++)
		const_rows[i] = row_store->get_row(i);
	tmp = std::shared_ptr<local_matrix_store>(
			new local_cref_col_matrix_store(const_cols, 0, 0,
				row_store->get_num_rows(), const_cols.size(),
				get_scalar_type<int>(), -1));
	verify_set(tmp);
	tmp = std::shared_ptr<local_matrix_store>(
			new local_cref_row_matrix_store(const_rows, 0, 0,
				const_rows.size(), row_store->get_num_cols(),
				get_scalar_type<int>(), -1));
	verify_set(tmp);
}

void test_agg1(std::shared_ptr<local_matrix_store> store)
{
	if (store->store_layout() == matrix_layout_t::L_COL)
		store->set_data(set_col_operate<int>(store->get_num_cols()));
	else
		store->set_data(set_row_operate<int>(store->get_num_cols()));

	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			store->get_type().get_basic_ops().get_add());
	for (size_t i = 0; i < 2; i++) {
		agg_operate::const_ptr op;
		if (i == 0)
			op = agg_operate::create(add, add);
		else
			op = agg_operate::create(add, bulk_operate::const_ptr());

		int sum = 0;
		local_ref_contig_col_matrix_store res((char *) &sum, 0, 0, 1, 1,
				get_scalar_type<int>(), -1);
		aggregate(*store, *op, matrix_margin::BOTH,
				store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
				res);
		int num_eles = store->get_num_rows() * store->get_num_cols();
		assert(sum == (num_eles - 1) * num_eles / 2);

		local_buf_col_matrix_store res1(0, 0, store->get_num_rows(), 1,
				op->get_output_type(), -1);
		aggregate(*store, *op, matrix_margin::MAR_ROW,
				store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
				res1);
		for (size_t i = 0; i < res1.get_num_rows(); i++) {
			int base = i * store->get_num_cols();
			assert(*(int *) res1.get(i, 0) == base * store->get_num_cols()
					+ (store->get_num_cols() - 1) * store->get_num_cols() / 2);
		}

		local_buf_col_matrix_store res2(0, 0, store->get_num_cols(), 1,
				op->get_output_type(), -1);
		aggregate(*store, *op, matrix_margin::MAR_COL,
				store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
				res2);
		for (size_t i = 0; i < res2.get_num_rows(); i++) {
			size_t n = store->get_num_rows();
			size_t m = store->get_num_cols();
			assert(*(int *) res2.get(i, 0) == i * n + (n - 1) * m * n / 2);
		}
	}
}

void test_aggregate(size_t long_dim)
{
	printf("test aggregate on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_agg1(std::shared_ptr<local_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_agg1(std::shared_ptr<local_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0,
					rows.size(), row_store->get_num_cols(),
					get_scalar_type<int>(), -1)));
}

void test_mapply21(std::shared_ptr<local_matrix_store> store)
{
	std::shared_ptr<local_matrix_store> res;
	if (store->store_layout() == matrix_layout_t::L_COL) {
		store->set_data(set_col_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}
	else {
		store->set_data(set_row_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}

	const bulk_operate &op = store->get_type().get_basic_ops().get_add();
	mapply2(*store, *store, op,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) * 2 == res->get<int>(i, j));

	if (store->store_layout() == matrix_layout_t::L_COL)
		store = local_matrix_store::ptr(new ltest_col_matrix_store(
					local_col_matrix_store::cast(store)));
	else
		store = local_matrix_store::ptr(new ltest_row_matrix_store(
					local_row_matrix_store::cast(store)));

	mapply2(*store, *store, op,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) * 2 == res->get<int>(i, j));
}

void test_mapply2(size_t long_dim)
{
	printf("test mapply2 on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_mapply21(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_mapply21(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_mapply21(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));
	test_mapply21(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_mapply21(std::shared_ptr<local_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_mapply21(std::shared_ptr<local_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0,
					rows.size(), row_store->get_num_cols(),
					get_scalar_type<int>(), -1)));
}

void test_sapply1(std::shared_ptr<local_matrix_store> store)
{
	std::shared_ptr<local_matrix_store> res;
	if (store->store_layout() == matrix_layout_t::L_COL) {
		store->set_data(set_col_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}
	else {
		store->set_data(set_row_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}

	const bulk_uoperate &op = *store->get_type().get_basic_uops().get_op(
			basic_uops::op_idx::NEG);
	sapply(*store, op,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == -res->get<int>(i, j));

	if (store->store_layout() == matrix_layout_t::L_COL)
		store = local_matrix_store::ptr(new ltest_col_matrix_store(
					local_col_matrix_store::cast(store)));
	else
		store = local_matrix_store::ptr(new ltest_row_matrix_store(
					local_row_matrix_store::cast(store)));

	sapply(*store, op,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == -res->get<int>(i, j));
}

void test_sapply(size_t long_dim)
{
	printf("test sapply on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_sapply1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_sapply1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_sapply1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));
	test_sapply1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_sapply1(std::shared_ptr<local_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_sapply1(std::shared_ptr<local_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0,
					rows.size(), row_store->get_num_cols(),
					get_scalar_type<int>(), -1)));
}

std::shared_ptr<local_matrix_store> naive_multiply(const local_matrix_store &m1,
		const local_matrix_store &m2)
{
	std::shared_ptr<local_matrix_store> res(new local_buf_col_matrix_store(
				0, 0, m1.get_num_rows(), m2.get_num_cols(),
				get_scalar_type<int>(), -1));
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m2.get_num_cols(); j++) {
			int sum = 0;
			for (size_t k = 0; k < m1.get_num_cols(); k++) {
				sum += m1.get<int>(i, k) * m2.get<int>(k, j);
			}
			res->set<int>(i, j, sum);
		}
	}
	return res;
}

void test_inner_prod1(std::shared_ptr<local_matrix_store> m1,
		std::shared_ptr<local_matrix_store> m2, matrix_layout_t res_layout)
{
	std::shared_ptr<local_matrix_store> res;
	if (m1->store_layout() == matrix_layout_t::L_COL)
		m1->set_data(set1_col_operate());
	else
		m1->set_data(set1_row_operate());

	if (res_layout == matrix_layout_t::L_COL)
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, m1->get_num_rows(), m2->get_num_cols(),
				get_scalar_type<int>(), -1));
	else
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, m1->get_num_rows(), m2->get_num_cols(),
				get_scalar_type<int>(), -1));

	if (m2->store_layout() == matrix_layout_t::L_COL)
		m2->set_data(set1_col_operate());
	else
		m2->set_data(set1_row_operate());

	// Multiply
	const scalar_type &type = m1->get_type();
	res->reset_data();
	if (m1->is_wide())
		inner_prod_wide(*m1, *m2, type.get_basic_ops().get_multiply(),
				type.get_basic_ops().get_add(), *res);
	else
		inner_prod_tall(*m1, *m2, type.get_basic_ops().get_multiply(),
				type.get_basic_ops().get_add(), *res);
	std::shared_ptr<local_matrix_store> res1 = naive_multiply(*m1, *m2);
	for (size_t i = 0; i < res->get_num_rows(); i++)
		for (size_t j = 0; j < res->get_num_cols(); j++)
			assert(res->get<int>(i, j) == res1->get<int>(i, j));
}

void test_inner_prod(size_t long_dim)
{
	printf("test inner product on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	std::shared_ptr<local_matrix_store> m1;
	std::shared_ptr<local_matrix_store> m2;

	printf("inner_prod tall col\n");
	m1 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	printf("inner_prod tall col\n");
	m1 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	printf("inner_prod tall row\n");
	m1 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	printf("inner_prod wide row\n");
	m1 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());
}

void verify_transpose(const local_matrix_store &m1, const local_matrix_store &m2)
{
	assert(m1.get_num_rows() == m2.get_num_cols());
	assert(m1.get_num_cols() == m2.get_num_rows());
	assert(m1.get_global_start_row() == m2.get_global_start_col());
	assert(m1.get_global_start_col() == m2.get_global_start_row());
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<int>(i, j) == m2.get<int>(j, i));
}

template<class MATRIX_TYPE>
void test_transpose1(std::shared_ptr<MATRIX_TYPE> m1)
{
	verify_transpose(*m1, *local_matrix_store::cast(m1->transpose()));
	m1->resize(2, 0, m1->get_num_rows() - 2, m1->get_num_cols());
	verify_transpose(*m1, *local_matrix_store::cast(m1->transpose()));
	m1->resize(0, 2, m1->get_num_rows(), m1->get_num_cols() - 2);
	verify_transpose(*m1, *local_matrix_store::cast(m1->transpose()));
	m1->resize(2, 2, m1->get_num_rows() - 2, m1->get_num_cols() - 2);
	verify_transpose(*m1, *local_matrix_store::cast(m1->transpose()));
	m1->reset_size();
}

void test_transpose(size_t long_dim)
{
	printf("test transpose on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	std::shared_ptr<local_col_matrix_store> m1;
	std::shared_ptr<local_row_matrix_store> m2;
	std::shared_ptr<local_col_matrix_store> m3;
	std::shared_ptr<local_row_matrix_store> m4;

	m1 = std::shared_ptr<local_col_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_transpose1(m1);

	m2 = std::shared_ptr<local_row_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	test_transpose1(m2);

	m3 = std::shared_ptr<local_col_matrix_store>(
			new local_ref_contig_col_matrix_store(m1->get_raw_arr(), 0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	test_transpose1(m3);

	m4 = std::shared_ptr<local_row_matrix_store>(
			new local_ref_contig_row_matrix_store(m1->get_raw_arr(), 0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	test_transpose1(m4);

	{
		std::vector<char *> cols(m1->get_num_cols());
		for (size_t i = 0; i < cols.size(); i++)
			cols[i] = m1->get_col(i);
		m3 = std::shared_ptr<local_col_matrix_store>(
				new local_ref_col_matrix_store(cols, 0, 0,
					m1->get_num_rows(), cols.size(), get_scalar_type<int>(), -1));
		test_transpose1(m3);
	}

	{
		std::vector<char *> rows(m2->get_num_rows());
		for (size_t i = 0; i < rows.size(); i++)
			rows[i] = m2->get_row(i);
		m4 = std::shared_ptr<local_row_matrix_store>(
				new local_ref_row_matrix_store(rows, 0, 0,
					rows.size(), m2->get_num_cols(), get_scalar_type<int>(), -1));
		test_transpose1(m4);
	}

	m3 = std::shared_ptr<local_col_matrix_store>(
			new local_cref_contig_col_matrix_store(m1->get_raw_arr(), 0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	test_transpose1(m3);

	m4 = std::shared_ptr<local_row_matrix_store>(
			new local_cref_contig_row_matrix_store(m1->get_raw_arr(), 0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	test_transpose1(m4);

	{
		std::vector<const char *> cols(m1->get_num_cols());
		for (size_t i = 0; i < cols.size(); i++)
			cols[i] = m1->get_col(i);
		m3 = std::shared_ptr<local_col_matrix_store>(
				new local_cref_col_matrix_store(cols, 0, 0,
					m1->get_num_rows(), cols.size(), get_scalar_type<int>(), -1));
		test_transpose1(m3);
	}

	{
		std::vector<const char *> rows(m2->get_num_rows());
		for (size_t i = 0; i < rows.size(); i++)
			rows[i] = m2->get_row(i);
		m4 = std::shared_ptr<local_row_matrix_store>(
				new local_cref_row_matrix_store(rows, 0, 0,
					rows.size(), m2->get_num_cols(), get_scalar_type<int>(), -1));
		test_transpose1(m4);
	}
}

void test_copy_from1(const local_matrix_store &m1, local_matrix_store &m2)
{
	m2.copy_from(m1);
	for (size_t i = 0; i < m1.get_num_rows(); i++) {
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<int>(i, j) == m2.get<int>(i, j));
	}
}

local_matrix_store::ptr get_cols(local_col_matrix_store::ptr store,
		size_t num_cols)
{
	assert(store->get_num_cols() >= num_cols);
	std::vector<char *> cols(num_cols);
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = store->get_col(i);
	return local_matrix_store::ptr(new local_ref_col_matrix_store(cols, 0, 0,
				store->get_num_rows(), cols.size(),
				store->get_type(), store->get_node_id()));
}

local_matrix_store::ptr get_rows(local_row_matrix_store::ptr store,
		size_t num_rows)
{
	assert(store->get_num_rows() >= num_rows);
	std::vector<char *> rows(num_rows);
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = store->get_row(i);
	return local_matrix_store::ptr(new local_ref_row_matrix_store(rows, 0, 0,
				rows.size(), store->get_num_cols(),
				store->get_type(), store->get_node_id()));
}

void test_copy_from(size_t long_dim, bool is_tall)
{
	printf("test copy from on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	local_matrix_store::ptr m1;
	local_matrix_store::ptr m2;
	local_matrix_store::ptr sub_m1;
	local_matrix_store::ptr sub_m2;

	// For the entire col-major matrix.
	if (is_tall) {
		m1 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
		m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	}
	else {
		m1 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0,  10, long_dim,get_scalar_type<int>(), -1));
		m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	}
	m1->set_data(set_col_operate<int>(m1->get_num_cols()));
	test_copy_from1(*m1, *m2);

	// For the sub col-major matrix.
	sub_m1 = get_cols(local_col_matrix_store::cast(m1), m1->get_num_cols() / 2);
	sub_m2 = get_cols(local_col_matrix_store::cast(m2), m2->get_num_cols() / 2);
	test_copy_from1(*sub_m1, *sub_m2);

	// For col-major 2 row-major.
	if (is_tall)
		m2 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	else
		m2 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);

	// For the entire row-major matrix.
	if (is_tall) {
		m1 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
		m2 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	}
	else {
		m1 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
		m2 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	}
	m1->set_data(set_row_operate<int>(m1->get_num_cols()));
	test_copy_from1(*m1, *m2);

	// For the sub row-major matrix.
	sub_m1 = get_rows(local_row_matrix_store::cast(m1), m1->get_num_rows() / 2);
	sub_m2 = get_rows(local_row_matrix_store::cast(m2), m2->get_num_rows() / 2);
	test_copy_from1(*sub_m1, *sub_m2);

	// For row-major 2 col-major
	if (is_tall)
		m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	else
		m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);
}

void test_copy_from(size_t long_dim)
{
	test_copy_from(long_dim, true);
	test_copy_from(long_dim, false);
}

void test_get_raw(size_t long_dim)
{
	// Test local buffer matrix store.
	local_matrix_store::ptr store(new local_buf_col_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	store->set_data(set_col_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(0, 1, store->get_num_rows(), store->get_num_cols() - 1);
	assert(store->get_raw_arr());

	store = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim, 1, get_scalar_type<int>(), -1));
	store->set_data(set_col_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);

	store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	store->set_data(set_row_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(1, 0, store->get_num_rows() - 1, store->get_num_cols());
	assert(store->get_raw_arr());

	store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				1, long_dim, get_scalar_type<int>(), -1));
	store->set_data(set_row_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);

	// Test local refernece matrix store.
	local_matrix_store::ptr buf(new local_buf_col_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	store = local_matrix_store::ptr(new local_ref_contig_col_matrix_store(
				buf->get_raw_arr(), 0, 0, buf->get_num_rows(), buf->get_num_cols(),
				buf->get_type(), buf->get_node_id()));
	store->set_data(set_col_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(0, 1, store->get_num_rows(), store->get_num_cols() - 1);
	assert(store->get_raw_arr());

	buf = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim, 1, get_scalar_type<int>(), -1));
	store = local_matrix_store::ptr(new local_ref_contig_col_matrix_store(
				buf->get_raw_arr(), 0, 0, buf->get_num_rows(), buf->get_num_cols(),
				buf->get_type(), buf->get_node_id()));
	store->set_data(set_col_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);

	buf = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	store = local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
				buf->get_raw_arr(), 0, 0, buf->get_num_rows(), buf->get_num_cols(),
				buf->get_type(), buf->get_node_id()));
	store->set_data(set_row_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(1, 0, store->get_num_rows() - 1, store->get_num_cols());
	assert(store->get_raw_arr());

	buf = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				1, long_dim, get_scalar_type<int>(), -1));
	store = local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
				buf->get_raw_arr(), 0, 0, buf->get_num_rows(), buf->get_num_cols(),
				buf->get_type(), buf->get_node_id()));
	store->set_data(set_row_operate<int>(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);
}

void test_conv_layout(size_t long_dim)
{
	// Tall matrices
	local_row_matrix_store::ptr row_store(new local_buf_row_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	row_store->set_data(set_row_operate<int>(row_store->get_num_cols()));
	local_col_matrix_store::ptr col_store(new local_buf_col_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	col_store->copy_from(*row_store);
	local_matrix_store::ptr store = col_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	row_store->reset_data();
	row_store->copy_from(*col_store);
	store = row_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	local_row_matrix_store::ptr row_buf = row_store;
	row_store = local_row_matrix_store::ptr(new local_ref_row_matrix_store(rows,
				0, 0, row_store->get_num_rows(), row_store->get_num_cols(),
				row_store->get_type(), -1));
	col_store->reset_data();
	col_store->copy_from(*row_store);
	store = col_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	row_store->reset_data();
	row_store->copy_from(*col_store);
	store = row_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	// Wide matrices
	row_store = local_row_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	row_store->set_data(set_row_operate<int>(row_store->get_num_cols()));
	col_store = local_col_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	col_store->copy_from(*row_store);
	store = col_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	row_store->reset_data();
	row_store->copy_from(*col_store);
	store = row_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	local_col_matrix_store::ptr col_buf = col_store;
	col_store = local_col_matrix_store::ptr(new local_ref_col_matrix_store(cols,
				0, 0, col_store->get_num_rows(), col_store->get_num_cols(),
				col_store->get_type(), -1));
	row_store->reset_data();
	row_store->copy_from(*col_store);
	store = row_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);

	col_store->reset_data();
	col_store->copy_from(*row_store);
	store = col_store;
	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == i * store->get_num_cols() + j);
}

void blas_col_multiply(const local_col_matrix_store &m1,
		const local_col_matrix_store &m2, local_col_matrix_store &res)
{
	assert(m1.get_type() == m2.get_type());
	local_col_matrix_store::ptr m1_buf;
	local_col_matrix_store::ptr m2_buf;
	local_col_matrix_store::ptr res_buf;
	const char *m1_raw = m1.get_raw_arr();
	if (m1_raw == NULL) {
		m1_buf = local_col_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				m1.get_num_rows(), m1.get_num_cols(), m1.get_type(), -1));
		m1_buf->copy_from(m1);
		m1_raw = m1_buf->get_raw_arr();
	}
	const char *m2_raw = m2.get_raw_arr();
	assert(m2_raw);
	char *res_raw = res.get_raw_arr();
	if (res_raw == NULL) {
		res_buf = local_col_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
					res.get_num_rows(), res.get_num_cols(), res.get_type(), -1));
		res_raw = res_buf->get_raw_arr();
	}
	if (m1.get_type() == get_scalar_type<double>()) {
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				m1.get_num_rows(), m2.get_num_cols(),
				m1.get_num_cols(), 1, (const double *) m1_raw,
				m1.get_num_rows(), (const double *) m2_raw,
				m2.get_num_rows(), 0, (double *) res_raw,
				res.get_num_rows());
	}
	else {
		assert(m1.get_type() == get_scalar_type<float>());
		cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				m1.get_num_rows(), m2.get_num_cols(),
				m1.get_num_cols(), 1, (const float *) m1_raw,
				m1.get_num_rows(), (const float *) m2_raw,
				m2.get_num_rows(), 0, (float *) res_raw,
				res.get_num_rows());
	}
	if (res_raw != res.get_raw_arr())
		res.copy_from(*res_buf);
}

template<class T>
void equal_mat(const local_matrix_store &m1, const local_matrix_store &m2)
{
	assert(m1.get_num_rows() == m2.get_num_rows());
	assert(m1.get_num_cols() == m2.get_num_cols());
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<T>(i, j) == m2.get<T>(i, j));
}

void blas_row_multiply(const local_row_matrix_store &m1,
		const local_row_matrix_store &m2, local_row_matrix_store &res)
{
	local_matrix_store::ptr col_m1 = m1.conv2(matrix_layout_t::L_COL);
	equal_mat<double>(*col_m1, m1);
	local_matrix_store::ptr col_m2 = m2.conv2(matrix_layout_t::L_COL);
	equal_mat<double>(*col_m2, m2);
	local_col_matrix_store::ptr col_res(new local_buf_col_matrix_store(0, 0,
				res.get_num_rows(), res.get_num_cols(), res.get_type(), -1));
	blas_col_multiply(dynamic_cast<const local_col_matrix_store &>(*col_m1),
			dynamic_cast<const local_col_matrix_store &>(*col_m2), *col_res);
	res.copy_from(*col_res);
}

void blas_multiply(const local_matrix_store &m1, const local_matrix_store &m2,
		local_matrix_store &res)
{
	if (m1.store_layout() == matrix_layout_t::L_ROW)
		blas_row_multiply(dynamic_cast<const local_row_matrix_store &>(m1),
				dynamic_cast<const local_row_matrix_store &>(m2),
				dynamic_cast<local_row_matrix_store &>(res));
	else
		blas_col_multiply(dynamic_cast<const local_col_matrix_store &>(m1),
				dynamic_cast<const local_col_matrix_store &>(m2),
				dynamic_cast<local_col_matrix_store &>(res));
}

void test_tall_multiply(size_t long_dim)
{
	local_matrix_store::ptr left_store, right_store, out1, out2;
	long_dim = 10;

	// test row matrix.
	printf("multiply tall row matrix (%ld)\n", long_dim);
	left_store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<double>(), -1));
	left_store->set_data(set_row_operate<double>(left_store->get_num_cols()));
	left_store->resize(0, 0, left_store->get_num_rows() - 1,
			left_store->get_num_cols());
	right_store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, 12, get_scalar_type<double>(), -1));
	right_store->set_data(set_row_operate<double>(right_store->get_num_cols()));
	out1 = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				long_dim, 12, get_scalar_type<double>(), -1));
	out1->resize(0, 0, out1->get_num_rows() - 1, out1->get_num_cols());
	out2 = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				long_dim - 1, 12, get_scalar_type<double>(), -1));

	std::pair<local_matrix_store::ptr, local_matrix_store::ptr> bufs;
	matrix_tall_multiply(*left_store, *right_store, *out1, bufs);
	blas_multiply(*left_store, *right_store, *out2);
	equal_mat<double>(*out1, *out2);

	printf("multiply tall col matrix (%ld)\n", long_dim);
	left_store = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<double>(), -1));
	left_store->resize(0, 0, left_store->get_num_rows() - 1,
			left_store->get_num_cols());
	left_store->set_data(set_col_operate<double>(left_store->get_num_cols()));
	right_store = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				10, 12, get_scalar_type<double>(), -1));
	right_store->set_data(set_col_operate<double>(right_store->get_num_cols()));
	out1 = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim, 12, get_scalar_type<double>(), -1));
	out1->resize(0, 0, out1->get_num_rows() - 1, out1->get_num_cols());
	out2 = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim - 1, 12, get_scalar_type<double>(), -1));

	matrix_tall_multiply(*left_store, *right_store, *out1, bufs);
	blas_multiply(*left_store, *right_store, *out2);
	equal_mat<double>(*out1, *out2);
}

void test_wide_multiply(size_t long_dim)
{
}

void test_multiply(size_t long_dim)
{
	test_tall_multiply(long_dim);
	test_wide_multiply(long_dim);
}

void verify_resize_tall(matrix_store::const_ptr store)
{
	local_matrix_store::const_ptr portion = store->get_portion(0);
	local_matrix_store *mutable_portion
		= const_cast<local_matrix_store *>(portion.get());
	bool ret = mutable_portion->resize(0, 0, 100, portion->get_num_cols());
	assert(ret);

	portion = store->get_portion(0);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	size_t portion_num_rows = portion->get_num_rows();
	size_t portion_num_cols = portion->get_num_cols();
	ret = mutable_portion->resize(0, 0, 100, portion->get_num_cols() - 1);
	assert(!ret);
	assert(portion->get_num_rows() == portion_num_rows);
	assert(portion->get_num_cols() == portion_num_cols);
	portion->materialize_self();

	portion = std::dynamic_pointer_cast<const local_matrix_store>(
			portion->transpose());
	assert(portion);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	ret = mutable_portion->resize(0, 0, portion->get_num_rows(), 100);
	assert(ret);

	portion = store->get_portion(0);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	portion_num_rows = portion->get_num_rows();
	portion_num_cols = portion->get_num_cols();
	ret = mutable_portion->resize(0, 0, portion->get_num_rows() - 1, 100);
	assert(!ret);
	assert(portion->get_num_rows() == portion_num_rows);
	assert(portion->get_num_cols() == portion_num_cols);
	portion->materialize_self();
}

void verify_resize_wide(matrix_store::const_ptr store)
{
	local_matrix_store::const_ptr portion = store->get_portion(0);
	local_matrix_store *mutable_portion
		= const_cast<local_matrix_store *>(portion.get());
	bool ret = mutable_portion->resize(0, 0, portion->get_num_rows(), 100);
	assert(ret);

	portion = store->get_portion(0);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	size_t portion_num_rows = portion->get_num_rows();
	size_t portion_num_cols = portion->get_num_cols();
	ret = mutable_portion->resize(0, 0, portion->get_num_rows() - 1, 100);
	assert(!ret);
	assert(portion->get_num_rows() == portion_num_rows);
	assert(portion->get_num_cols() == portion_num_cols);
	portion->materialize_self();

	portion = std::dynamic_pointer_cast<const local_matrix_store>(
			portion->transpose());
	assert(portion);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	ret = mutable_portion->resize(0, 0, 100, portion->get_num_cols());
	assert(ret);

	portion = store->get_portion(0);
	mutable_portion = const_cast<local_matrix_store *>(portion.get());
	portion_num_rows = portion->get_num_rows();
	portion_num_cols = portion->get_num_cols();
	ret = mutable_portion->resize(0, 0, 100, portion->get_num_cols() - 1);
	assert(!ret);
	assert(portion->get_num_rows() == portion_num_rows);
	assert(portion->get_num_cols() == portion_num_cols);
	portion->materialize_self();
}

void test_resize_multiply(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_COL);
	dense_matrix::ptr small = dense_matrix::create_randu<double>(0, 1,
			10, 9, matrix_layout_t::L_COL);

	printf("verify tall matrix multiply\n");
	dense_matrix::ptr res = mat->multiply(*small);
	verify_resize_tall(res->get_raw_store());
	res = res->transpose();
	verify_resize_wide(res->get_raw_store());

	printf("verify tall matrix inner product\n");
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			*mat->get_type().get_basic_ops().get_op(basic_ops::op_idx::ADD));
	bulk_operate::const_ptr mul = bulk_operate::conv2ptr(
			*mat->get_type().get_basic_ops().get_op(basic_ops::op_idx::MUL));
	res = mat->inner_prod(*small, mul, add);
	verify_resize_tall(res->get_raw_store());
	res = res->transpose();
	verify_resize_wide(res->get_raw_store());
}

void test_resize_mapply_rowcol(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_COL);
	col_vec::ptr vec = col_vec::create(dense_matrix::create_randu<double>(0, 1,
				mat->get_num_cols(), 1, matrix_layout_t::L_COL));
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			*mat->get_type().get_basic_ops().get_op(basic_ops::op_idx::ADD));

	printf("verify mapply rows\n");
	dense_matrix::ptr res = mat->mapply_rows(vec, add);
	verify_resize_tall(res->get_raw_store());

	printf("verify mapply cols\n");
	mat = mat->transpose();
	res = mat->mapply_cols(vec, add);
	verify_resize_wide(res->get_raw_store());
}

void test_resize_get_rowcols(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_ROW);
	mat = mat->add(*mat);
	std::vector<off_t> idxs(mat->get_num_cols() / 2);
	for (size_t i = 0; i < idxs.size(); i++)
		idxs[i] = random() % mat->get_num_cols();

	printf("verify get cols\n");
	dense_matrix::ptr res = mat->get_cols(idxs);
	verify_resize_tall(res->get_raw_store());

	printf("verify get rows\n");
	mat = mat->transpose();
	res = mat->get_rows(idxs);
	verify_resize_wide(res->get_raw_store());
}

void test_resize_repeat(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			10, 9, matrix_layout_t::L_COL);
	col_vec::ptr idxs = col_vec::create(dense_matrix::create_randu<size_t>(
				0, mat->get_num_rows() - 1, long_dim, 1,
				matrix_layout_t::L_COL));

	printf("verify repeat rows\n");
	dense_matrix::ptr res = mat->get_rows(idxs);
	verify_resize_tall(res->get_raw_store());

	printf("verify repeat cols\n");
	mat = mat->transpose();
	res = mat->get_cols(idxs);
	verify_resize_wide(res->get_raw_store());
}

class sum_apply_op: public arr_apply_operate
{
public:
	void run(const local_vec_store &in, local_vec_store &out) const {
		assert(in.get_type() == get_scalar_type<double>());
		assert(out.get_type() == get_scalar_type<double>());
		double res = 0;
		for (size_t i = 0; i < in.get_length(); i++)
			res += in.get<double>(i);
		out.set<double>(0, res);
	}

	const scalar_type &get_input_type() const {
		return get_scalar_type<double>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<double>();
	}
	virtual size_t get_num_out_eles(size_t num_input) const {
		return 1;
	}
};

void test_resize_apply(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			long_dim, 10, matrix_layout_t::L_COL);

	printf("verify apply rows\n");
	dense_matrix::ptr res = mat->apply(matrix_margin::MAR_ROW,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	verify_resize_tall(res->get_raw_store());

	printf("verify apply cols\n");
	mat = mat->transpose();
	res = mat->apply(matrix_margin::MAR_COL,
			arr_apply_operate::const_ptr(new sum_apply_op()));
	verify_resize_wide(res->get_raw_store());
}

void test_resize_groupby(size_t long_dim)
{
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			10, long_dim, matrix_layout_t::L_COL);
	factor_col_vector::ptr labels = factor_col_vector::create(factor(3),
			dense_matrix::create_randu<int>(0, 2, mat->get_num_rows(), 1,
				matrix_layout_t::L_COL));
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			*mat->get_type().get_basic_ops().get_op(basic_ops::op_idx::ADD));
	printf("verify group rows\n");
	dense_matrix::ptr res = mat->groupby_row(labels, add);
	verify_resize_wide(res->get_raw_store());
}

// This is to test if we can resize a local matrix as we expect.
void test_resize(size_t long_dim)
{
	test_resize_groupby(long_dim);
	test_resize_multiply(long_dim);
	test_resize_mapply_rowcol(long_dim);
	test_resize_get_rowcols(long_dim);
	test_resize_repeat(long_dim);
	test_resize_apply(long_dim);
}

class nonresize_copy_op: public detail::portion_mapply_op
{
	size_t resize_dim;
public:
	nonresize_copy_op(size_t num_rows, size_t num_cols, const scalar_type &type,
			size_t resize_dim): detail::portion_mapply_op(
				num_rows, num_cols, type) {
		this->resize_dim = resize_dim;
	}

	virtual bool is_resizable(size_t local_start_row, size_t local_start_col,
			size_t local_num_rows, size_t local_num_cols) const {
		if (resize_dim == 1)
			// non-resizable on cols.
			return local_num_rows == get_out_num_rows();
		else if (resize_dim == 2)
			// non-resizable on rows.
			return local_num_cols == get_out_num_cols();
		else
			return true;
	}

	virtual void run(
			const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		out.copy_from(*ins[0]);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		assert(0);
		return detail::portion_mapply_op::const_ptr();
	}
	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string();
	}
};

matrix_store::const_ptr get_data_mat(size_t nrow, size_t ncol,
		matrix_layout_t layout, size_t resize_dim)
{
	matrix_store::ptr in = mem_matrix_store::create(nrow, ncol, layout,
			get_scalar_type<double>(), -1);
	scalar_variable_impl<double> start(0);
	scalar_variable_impl<double> stride(1);
	auto set = get_scalar_type<double>().get_set_seq(start, stride, nrow, ncol,
			true, layout);
	in->set_data(*set);
	if (resize_dim == 1 || resize_dim == 2)
		return matrix_store::const_ptr(new mapply_matrix_store(
					std::vector<matrix_store::const_ptr>(1, in),
					portion_mapply_op::const_ptr(new nonresize_copy_op(nrow,
							ncol, in->get_type(), resize_dim)),
					layout));
	else
		return in;
}

void verify_equal(const matrix_store &m1, const matrix_store &m2)
{
	const mem_matrix_store &mem_m1 = dynamic_cast<const mem_matrix_store &>(m1);
	const mem_matrix_store &mem_m2 = dynamic_cast<const mem_matrix_store &>(m2);
	assert(memcmp(mem_m1.get_raw_arr(), mem_m2.get_raw_arr(),
				m1.get_num_rows() * m1.get_num_cols() * m1.get_entry_size()) == 0);
}

void test_resize_compute_multiply()
{
	matrix_store::const_ptr m1, m2;
	matrix_store::ptr correct_res, res;

	printf("verify resize compute on inner prod wide\n");
	// resizable local matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	const bulk_operate &add = m1->get_type().get_basic_ops().get_add();
	const bulk_operate &mul = m1->get_type().get_basic_ops().get_multiply();
	correct_res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_wide(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*correct_res->get_portion(0));

	// non-resizable on rows in the left matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_wide(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	// non-resizable on cols in the right matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 1);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_wide(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	// non-resizable on rows in the left matrix and cols in the right matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 1);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_wide(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	printf("verify resize compute on inner prod tall\n");
	// resizable local matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	correct_res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_tall(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*correct_res->get_portion(0));

	// non-resizable on rows in the left matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	inner_prod_tall(*m1->get_portion(0), *m2->get_portion(0), mul, add,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	std::pair<local_matrix_store::ptr, local_matrix_store::ptr> bufs;
	printf("verify resize compute on multiply tall\n");
	m1 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	correct_res = mem_matrix_store::create(500, 500, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	matrix_tall_multiply(*m1->get_portion(0), *m2->get_portion(0),
			*correct_res->get_portion(0), bufs);

	m1 = get_data_mat(500, 500, matrix_layout_t::L_COL, 1);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	matrix_tall_multiply(*m1->get_portion(0), *m2->get_portion(0),
			*res->get_portion(0), bufs);
	verify_equal(*correct_res, *res);
}

void test_resize_compute_agg()
{
	matrix_store::const_ptr m;
	matrix_store::ptr correct_res, res;

	printf("verify resize compute on agg both\n");
	// resizable local matrix.
	m = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			m->get_type().get_basic_ops().get_add());
	agg_operate::const_ptr agg_add = agg_operate::create(add, add);
	correct_res = mem_matrix_store::create(1, 1, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	aggregate(*m->get_portion(0), *agg_add, matrix_margin::BOTH,
			part_dim_t::PART_DIM1, *correct_res->get_portion(0));

	// non-resizable on cols.
	m = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	res = mem_matrix_store::create(1, 1, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	aggregate(*m->get_portion(0), *agg_add, matrix_margin::BOTH,
			part_dim_t::PART_DIM2, *res->get_portion(0));
	verify_equal(*correct_res, *res);

	// non-resizable on rows.
	m = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	res = mem_matrix_store::create(1, 1, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	aggregate(*m->get_portion(0), *agg_add, matrix_margin::BOTH,
			part_dim_t::PART_DIM1, *res->get_portion(0));
	verify_equal(*correct_res, *res);
}

void test_resize_compute_mapply2()
{
	matrix_store::const_ptr m1, m2;
	matrix_store::ptr correct_res, res;

	printf("verify resize compute on mapply2\n");
	// resizable local matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	const bulk_operate &add = m1->get_type().get_basic_ops().get_add();
	correct_res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	mapply2(*m1->get_portion(0), *m2->get_portion(0), add, part_dim_t::PART_NONE,
			*correct_res->get_portion(0));

	// non-resizable on cols.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	mapply2(*m1->get_portion(0), *m2->get_portion(0), add, part_dim_t::PART_DIM2,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	mapply2(*m1->get_portion(0), *m2->get_portion(0), add, part_dim_t::PART_DIM2,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	// non-resizable on rows.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	mapply2(*m1->get_portion(0), *m2->get_portion(0), add, part_dim_t::PART_DIM1,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	m2 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	mapply2(*m1->get_portion(0), *m2->get_portion(0), add, part_dim_t::PART_DIM1,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);
}

void test_resize_compute_sapply()
{
	matrix_store::const_ptr m1;
	matrix_store::ptr correct_res, res;

	printf("verify resize compute on sapply\n");
	// resizable local matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	const bulk_uoperate &neg = *m1->get_type().get_basic_uops().get_op(
			basic_uops::op_idx::NEG);
	correct_res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	sapply(*m1->get_portion(0), neg, part_dim_t::PART_NONE,
			*correct_res->get_portion(0));

	// non-resizable on cols.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	sapply(*m1->get_portion(0), neg, part_dim_t::PART_DIM2,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);

	// non-resizable on rows.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	res = mem_matrix_store::create(500, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	sapply(*m1->get_portion(0), neg, part_dim_t::PART_DIM1,
			*res->get_portion(0));
	verify_equal(*correct_res, *res);
}

void test_resize_compute_groupby()
{
	matrix_store::const_ptr m1;
	matrix_store::ptr correct_res, res;


	printf("verify resize compute on groupby\n");
	// groupby on rows, resizable local matrix.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 0);
	bulk_operate::const_ptr add = bulk_operate::conv2ptr(
			m1->get_type().get_basic_ops().get_add());
	agg_operate::const_ptr agg_add = agg_operate::create(add, add);
	dense_matrix::ptr labels = factor_col_vector::create(factor(3),
			dense_matrix::create_randu<int>(0, 2, m1->get_num_rows(), 1,
				matrix_layout_t::L_COL));
	std::vector<bool> agg_flags(m1->get_num_rows());
	correct_res = mem_matrix_store::create(3, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_ROW, part_dim_t::PART_NONE,
			*correct_res->get_portion(0), agg_flags);

	// non-resizable on cols.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 2);
	res = mem_matrix_store::create(3, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	agg_flags.clear();
	agg_flags.resize(m1->get_num_rows());
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_ROW, part_dim_t::PART_DIM2,
			*res->get_portion(0), agg_flags);
	verify_equal(*correct_res, *res);

	// non-resizable on rows.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_ROW, 1);
	res = mem_matrix_store::create(3, 500, matrix_layout_t::L_ROW,
			get_scalar_type<double>(), -1);
	agg_flags.clear();
	agg_flags.resize(m1->get_num_rows());
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_ROW, part_dim_t::PART_DIM1,
			*res->get_portion(0), agg_flags);
	verify_equal(*correct_res, *res);

	// groupby on cols
	m1 = get_data_mat(500, 500, matrix_layout_t::L_COL, 0);
	labels = labels->transpose();
	correct_res = mem_matrix_store::create(500, 3, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	agg_flags.clear();
	agg_flags.resize(m1->get_num_cols());
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_COL, part_dim_t::PART_DIM2,
			*correct_res->get_portion(0), agg_flags);

	// non-resizable on cols.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_COL, 2);
	res = mem_matrix_store::create(500, 3, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	agg_flags.clear();
	agg_flags.resize(m1->get_num_cols());
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_COL, part_dim_t::PART_DIM2,
			*res->get_portion(0), agg_flags);
	verify_equal(*correct_res, *res);

	// non-resizable on rows.
	m1 = get_data_mat(500, 500, matrix_layout_t::L_COL, 1);
	res = mem_matrix_store::create(500, 3, matrix_layout_t::L_COL,
			get_scalar_type<double>(), -1);
	agg_flags.clear();
	agg_flags.resize(m1->get_num_cols());
	groupby(*labels->get_raw_store()->get_portion(0), *m1->get_portion(0),
			*agg_add, matrix_margin::MAR_COL, part_dim_t::PART_DIM1,
			*res->get_portion(0), agg_flags);
	verify_equal(*correct_res, *res);
}

// This is to test if we can still get the correct results when resize fails.
void test_resize_compute()
{
	test_resize_compute_multiply();
	test_resize_compute_agg();
	test_resize_compute_mapply2();
	test_resize_compute_sapply();
	test_resize_compute_groupby();
}

void test_cum1(std::shared_ptr<local_matrix_store> store)
{
	printf("test cum on %ldx%ld matrix, layout: %s\n", store->get_num_rows(),
			store->get_num_cols(),
			store->store_layout() == matrix_layout_t::L_COL ? "col" : "row");
	std::shared_ptr<local_matrix_store> res;
	if (store->store_layout() == matrix_layout_t::L_COL) {
		store->set_data(set_col_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				store->get_type(), -1));
	}
	else {
		store->set_data(set_row_operate<int>(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				store->get_type(), -1));
	}

	const agg_operate &op = *store->get_type().get_agg_ops().get_op(
			agg_ops::op_idx::SUM);
	cum(*store, NULL, op, matrix_margin::MAR_ROW,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++) {
		for (size_t j = 1; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == res->get<int>(i, j) - res->get<int>(i, j - 1));
		assert(store->get<int>(i, 0) == res->get<int>(i, 0));
	}

	cum(*store, NULL, op, matrix_margin::MAR_COL,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t j = 0; j < store->get_num_cols(); j++) {
		for (size_t i = 1; i < store->get_num_rows(); i++)
			assert(store->get<int>(i, j) == res->get<int>(i, j) - res->get<int>(i - 1, j));
		assert(store->get<int>(0, j) == res->get<int>(0, j));
	}

	// Test on virtual local matrix
	if (store->store_layout() == matrix_layout_t::L_COL)
		store = local_matrix_store::ptr(new ltest_col_matrix_store(
					local_col_matrix_store::cast(store)));
	else
		store = local_matrix_store::ptr(new ltest_row_matrix_store(
					local_row_matrix_store::cast(store)));
	cum(*store, NULL, op, matrix_margin::MAR_ROW,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t i = 0; i < store->get_num_rows(); i++) {
		for (size_t j = 1; j < store->get_num_cols(); j++)
			assert(store->get<int>(i, j) == res->get<int>(i, j) - res->get<int>(i, j - 1));
		assert(store->get<int>(i, 0) == res->get<int>(i, 0));
	}

	cum(*store, NULL, op, matrix_margin::MAR_COL,
			store->is_wide() ? detail::part_dim_t::PART_DIM2 : detail::part_dim_t::PART_DIM1,
			*res);
	for (size_t j = 0; j < store->get_num_cols(); j++) {
		for (size_t i = 1; i < store->get_num_rows(); i++)
			assert(store->get<int>(i, j) == res->get<int>(i, j) - res->get<int>(i - 1, j));
		assert(store->get<int>(0, j) == res->get<int>(0, j));
	}
}

void test_cum(size_t long_dim)
{
	printf("test cumsum on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_cum1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_cum1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_cum1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));
	test_cum1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1)));

	// Test on local reference matrix to a matrix stored non-contiguously.
	std::shared_ptr<local_col_matrix_store> col_store(new local_buf_col_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::shared_ptr<local_row_matrix_store> row_store(new local_buf_row_matrix_store(
				0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	std::vector<char *> cols(col_store->get_num_cols());
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = col_store->get_col(i);
	std::vector<char *> rows(row_store->get_num_rows());
	for (size_t i = 0; i < rows.size(); i++)
		rows[i] = row_store->get_row(i);
	test_cum1(std::shared_ptr<local_matrix_store>(new local_ref_col_matrix_store(
					cols, 0, 0, col_store->get_num_rows(), cols.size(),
					get_scalar_type<int>(), -1)));
	test_cum1(std::shared_ptr<local_matrix_store>(new local_ref_row_matrix_store(
					rows, 0, 0, rows.size(), row_store->get_num_cols(),
					get_scalar_type<int>(), -1)));
}

int main()
{
	test_cum(10000);
	test_resize(10000);
	test_resize_compute();
	test_multiply(1000);
	test_multiply(10000);
	test_conv_layout(1000);
	test_get_raw(1000);
	test_reset(1000);
	test_reset(10000);
	test_set(1000);
	test_set(10000);
	test_aggregate(1000);
	test_aggregate(10000);
	test_mapply2(1000);
	test_mapply2(10000);
	test_sapply(1000);
	test_sapply(10000);
	test_inner_prod(1000);
	test_inner_prod(10000);
	test_transpose(10000);
	test_copy_from(1000);
}
