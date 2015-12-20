#include "bulk_operate.h"
#include "local_matrix_store.h"
#include "dense_matrix.h"

using namespace fm;
using namespace fm::detail;

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
};

class set1_col_operate: public type_set_operate<int>
{
public:
	set1_col_operate() { }

	void set(int *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = 1;
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
		store->set_data(set_col_operate(store->get_num_cols()));
	else
		store->set_data(set_row_operate(store->get_num_cols()));
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
		store->set_data(set_col_operate(store->get_num_cols()));
	else
		store->set_data(set_row_operate(store->get_num_cols()));

	const bulk_operate &op = store->get_type().get_basic_ops().get_add();
	int sum = 0;
	local_ref_vec_store res((char *) &sum, 0, 1, get_scalar_type<int>(), -1);
	aggregate(*store, op, matrix_margin::BOTH, res);
	int num_eles = store->get_num_rows() * store->get_num_cols();
	assert(sum == (num_eles - 1) * num_eles / 2);
}

void test_aggregate(size_t long_dim)
{
	printf("test aggregate on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));
	test_agg1(std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1)));

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
		store->set_data(set_col_operate(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}
	else {
		store->set_data(set_row_operate(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}

	const bulk_operate &op = store->get_type().get_basic_ops().get_add();
	mapply2(*store, *store, op, *res);
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
		store->set_data(set_col_operate(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}
	else {
		store->set_data(set_row_operate(store->get_num_cols()));
		res = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
				0, 0, store->get_num_rows(), store->get_num_cols(),
				get_scalar_type<int>(), -1));
	}

	const bulk_uoperate &op = *store->get_type().get_basic_uops().get_op(
			basic_uops::op_idx::NEG);
	sapply(*store, op, *res);
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
	inner_prod(*m1, *m2, type.get_basic_ops().get_multiply(),
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

	m1 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	m1 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	m1 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	m1 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, m1->store_layout());

	m1 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, 10, long_dim, get_scalar_type<int>(), -1));
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_inner_prod1(m1, m2, matrix_layout_t::L_ROW);
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

void test_copy_from(size_t long_dim)
{
	printf("test copy from on local matrix, long dim: %ld\n", long_dim);
	// Test on local buffer matrix.
	local_matrix_store::ptr m1;
	local_matrix_store::ptr m2;
	local_matrix_store::ptr sub_m1;
	local_matrix_store::ptr sub_m2;

	// For the entire col-major matrix.
	m1 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m1->set_data(set_col_operate(m1->get_num_cols()));
	m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);

	// For the sub col-major matrix.
	sub_m1 = get_cols(local_col_matrix_store::cast(m1), m1->get_num_cols() / 2);
	sub_m2 = get_cols(local_col_matrix_store::cast(m2), m2->get_num_cols() / 2);
	test_copy_from1(*sub_m1, *sub_m2);

	// For col-major 2 row-major.
	m2 = std::shared_ptr<local_matrix_store>(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);

	// For the entire row-major matrix.
	m1 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	m1->set_data(set_row_operate(m1->get_num_cols()));
	m2 = local_matrix_store::ptr(new local_buf_row_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);

	// For the sub row-major matrix.
	sub_m1 = get_rows(local_row_matrix_store::cast(m1), m1->get_num_rows() / 2);
	sub_m2 = get_rows(local_row_matrix_store::cast(m2), m2->get_num_rows() / 2);
	test_copy_from1(*sub_m1, *sub_m2);

	// For row-major 2 col-major
	m2 = local_matrix_store::ptr(new local_buf_col_matrix_store(
					0, 0, long_dim, 10, get_scalar_type<int>(), -1));
	test_copy_from1(*m1, *m2);
}

void test_get_raw(size_t long_dim)
{
	// Test local buffer matrix store.
	local_matrix_store::ptr store(new local_buf_col_matrix_store(0, 0,
				long_dim, 10, get_scalar_type<int>(), -1));
	store->set_data(set_col_operate(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(0, 1, store->get_num_rows(), store->get_num_cols() - 1);
	assert(store->get_raw_arr());

	store = local_matrix_store::ptr(new local_buf_col_matrix_store(0, 0,
				long_dim, 1, get_scalar_type<int>(), -1));
	store->set_data(set_col_operate(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);

	store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	store->set_data(set_row_operate(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr() == NULL);
	store->reset_size();
	store->resize(1, 0, store->get_num_rows() - 1, store->get_num_cols());
	assert(store->get_raw_arr());

	store = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				1, long_dim, get_scalar_type<int>(), -1));
	store->set_data(set_row_operate(store->get_num_cols()));
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
	store->set_data(set_col_operate(store->get_num_cols()));
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
	store->set_data(set_col_operate(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(100, 0, store->get_num_rows() - 100, store->get_num_cols());
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);

	buf = local_matrix_store::ptr(new local_buf_row_matrix_store(0, 0,
				10, long_dim, get_scalar_type<int>(), -1));
	store = local_matrix_store::ptr(new local_ref_contig_row_matrix_store(
				buf->get_raw_arr(), 0, 0, buf->get_num_rows(), buf->get_num_cols(),
				buf->get_type(), buf->get_node_id()));
	store->set_data(set_row_operate(store->get_num_cols()));
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
	store->set_data(set_row_operate(store->get_num_cols()));
	assert(store->get_raw_arr());
	store->resize(0, 100, store->get_num_rows(), store->get_num_cols() - 100);
	assert(store->get_raw_arr());
	assert(*(int *) store->get_raw_arr() == 100);
}

int main()
{
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
