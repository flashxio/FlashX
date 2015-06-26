#include <malloc.h>

#include "EM_dense_matrix.h"
#include "sparse_matrix.h"
#include "dense_matrix.h"

using namespace fm;

class set_col_operate: public type_set_operate<long>
{
	size_t num_cols;
public:
	set_col_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(long *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = (row_idx + i) * num_cols + col_idx;
		}
	}
};

class set_row_operate: public type_set_operate<long>
{
	size_t num_cols;
public:
	set_row_operate(size_t num_cols) {
		this->num_cols = num_cols;
	}

	void set(long *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++) {
			arr[i] = row_idx * num_cols + col_idx + i;
		}
	}
};

void verify_portion(const detail::local_matrix_store &portion,
		size_t tot_num_cols)
{
	if (portion.store_layout() == matrix_layout_t::L_ROW) {
		const detail::local_row_matrix_store &rows
			= dynamic_cast<const detail::local_row_matrix_store &>(portion);
		for (size_t i = 0; i < rows.get_num_rows(); i++) {
			const long *row = (const long *) rows.get_row(i);
			for (size_t j = 0; j < rows.get_num_cols(); j++) {
				if (row[j] != (rows.get_global_start_row()
							+ i) * tot_num_cols + rows.get_global_start_col() + j)
					printf("%ld, %ld: %ld\n", rows.get_global_start_row() + i,
							rows.get_global_start_col() + j, row[j]);
				assert(row[j] == (rows.get_global_start_row()
							+ i) * tot_num_cols + rows.get_global_start_col() + j);
			}
		}
	}
	else {
		const detail::local_col_matrix_store &cols
			= dynamic_cast<const detail::local_col_matrix_store &>(portion);
		for (size_t i = 0; i < cols.get_num_cols(); i++) {
			const long *col = (const long *) cols.get_col(i);
			for (size_t j = 0; j < cols.get_num_rows(); j++)
				assert(col[j] == (cols.get_global_start_row()
							+ j) * tot_num_cols + cols.get_global_start_col() + i);
		}
	}
}

void verify_matrix(detail::EM_matrix_store::ptr store)
{
	std::pair<size_t, size_t> portion_size = store->get_portion_size();
	if (store->is_wide()) {
		for (size_t off = 0; off < store->get_num_cols();
				off += portion_size.second) {
			size_t num_cols = std::min(portion_size.second,
					store->get_num_cols() - off);
			detail::local_matrix_store::ptr portion = store->get_portion(
					0, off, store->get_num_rows(), num_cols);
			verify_portion(*portion, store->get_num_cols());
		}
	}
	else {
		for (size_t off = 0; off < store->get_num_rows();
				off += portion_size.first) {
			size_t num_rows = std::min(portion_size.first,
					store->get_num_rows() - off);
			detail::local_matrix_store::ptr portion = store->get_portion(
					off, 0, num_rows, store->get_num_cols());
			verify_portion(*portion, store->get_num_cols());
		}
	}
}

void test_set_data()
{
	detail::EM_matrix_store::ptr mat;

	printf("test setdata a row tall matrix\n");
	mat = detail::EM_matrix_store::create(9999999, 10, matrix_layout_t::L_ROW,
			get_scalar_type<long>());
	mat->set_data(set_row_operate(mat->get_num_cols()));
	verify_matrix(mat);

	printf("test setdata a column tall matrix\n");
	mat = detail::EM_matrix_store::create(9999999, 10, matrix_layout_t::L_COL,
			get_scalar_type<long>());
	mat->set_data(set_col_operate(mat->get_num_cols()));
	verify_matrix(mat);

	printf("test setdata a row wide matrix\n");
	mat = detail::EM_matrix_store::create(10, 9999999, matrix_layout_t::L_ROW,
			get_scalar_type<long>());
	mat->set_data(set_row_operate(mat->get_num_cols()));
	verify_matrix(mat);

	printf("test setdata a column wide matrix\n");
	mat = detail::EM_matrix_store::create(10, 9999999, matrix_layout_t::L_COL,
			get_scalar_type<long>());
	mat->set_data(set_col_operate(mat->get_num_cols()));
	verify_matrix(mat);
}

void test_get_portion()
{
	// Test getting unaligned portions.
}

void test_inner_prod()
{
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	test_set_data();
	test_get_portion();
	test_inner_prod();

	destroy_flash_matrix();
}
