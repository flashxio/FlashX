#include <set>

#include "bulk_operate.h"
#include "mem_matrix_store.h"

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
	for (size_t i = 0; i < idxs.size(); i++)
		assert(memcmp(sub_store->get_col(i), store->get_col(idxs[i]),
				sub_store->get_entry_size() * sub_store->get_num_rows()) == 0);

	std::vector<off_t> idxs2(2);
	idxs2[0] = 0;
	idxs2[1] = 1;
	mem_col_matrix_store::const_ptr subsub_store = mem_col_matrix_store::cast(
			sub_store->get_cols(idxs2));
	assert(subsub_store != NULL);
	assert(subsub_store->get_num_cols() == idxs2.size());
	for (size_t i = 0; i < idxs2.size(); i++)
		assert(memcmp(subsub_store->get_col(i), store->get_col(idxs[idxs2[i]]),
					store->get_entry_size() * store->get_num_rows()) == 0);
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
	for (size_t i = 0; i < idxs.size(); i++)
		assert(memcmp(sub_store->get_row(i), store->get_row(idxs[i]),
				sub_store->get_entry_size() * sub_store->get_num_cols()) == 0);

	idxs_set.clear();
	for (size_t i = 0; i < sub_store->get_num_rows() / 2; i++)
		idxs_set.insert(random() % sub_store->get_num_rows());
	std::vector<off_t> idxs2(idxs_set.begin(), idxs_set.end());
	mem_row_matrix_store::const_ptr subsub_store = mem_row_matrix_store::cast(
			sub_store->get_rows(idxs2));
	assert(subsub_store != NULL);
	assert(subsub_store->get_num_rows() == idxs2.size());
	for (size_t i = 0; i < idxs2.size(); i++)
		assert(memcmp(subsub_store->get_row(i), store->get_row(idxs[idxs2[i]]),
					store->get_entry_size() * store->get_num_cols()) == 0);
}

void test_io()
{
	printf("test read/write matrix to a file\n");
	std::string out_file;
	mem_matrix_store::ptr read_mat;
	mem_matrix_store::ptr orig_mat;

	orig_mat = mem_row_matrix_store::create(
				1000, 10, get_scalar_type<int>());
	orig_mat->set_data(set_row_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file);
	assert(read_mat);
	assert(read_mat->store_layout() == matrix_layout_t::L_ROW);
	for (size_t i = 0; i < orig_mat->get_num_rows(); i++) {
		for (size_t j = 0; j < orig_mat->get_num_cols(); j++)
			assert(orig_mat->get<int>(i, j) == read_mat->get<int>(i, j));
	}
	unlink(out_file.c_str());

	orig_mat = mem_col_matrix_store::create(
				1000, 10, get_scalar_type<int>());
	orig_mat->set_data(set_col_operate(orig_mat->get_num_cols()));
	out_file = tmpnam(NULL);
	orig_mat->write2file(out_file);
	read_mat = mem_matrix_store::load(out_file);
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

int main()
{
	test_reset(1000);
	test_reset(1000000);
	test_set(1000);
	test_set(1000000);
	test_sub_col_matrix();
	test_sub_row_matrix();
	test_io();
	test_transpose();
}
