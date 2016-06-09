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

void verify_matrix(detail::EM_matrix_store::const_ptr store)
{
	std::pair<size_t, size_t> portion_size = store->get_portion_size();
	if (store->is_wide()) {
		for (size_t off = 0; off < store->get_num_cols();
				off += portion_size.second) {
			size_t num_cols = std::min(portion_size.second,
					store->get_num_cols() - off);
			detail::local_matrix_store::const_ptr portion = store->get_portion(
					0, off, store->get_num_rows(), num_cols);
			verify_portion(*portion, store->get_num_cols());
		}
	}
	else {
		for (size_t off = 0; off < store->get_num_rows();
				off += portion_size.first) {
			size_t num_rows = std::min(portion_size.first,
					store->get_num_rows() - off);
			detail::local_matrix_store::const_ptr portion = store->get_portion(
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

void test_get_row_col()
{
	printf("accessing a row/column\n");
	detail::EM_matrix_store::ptr mat;
	off_t idx = 1;

	mat = detail::EM_matrix_store::create(9999999, 10, matrix_layout_t::L_COL,
			get_scalar_type<long>());
	mat->set_data(set_col_operate(mat->get_num_cols()));
	detail::smp_vec_store::const_ptr col = detail::smp_vec_store::cast(
			mat->get_col_vec(idx));
	assert(col->get_length() == mat->get_num_rows());
	assert(col->get_type() == mat->get_type());
	for (size_t i = 0; i < col->get_length(); i++)
		assert(col->get<long>(i) == i * mat->get_num_cols() + idx);

	mat = detail::EM_matrix_store::create(10, 9999999, matrix_layout_t::L_ROW,
			get_scalar_type<long>());
	mat->set_data(set_row_operate(mat->get_num_cols()));
	detail::smp_vec_store::const_ptr row = detail::smp_vec_store::cast(
			mat->get_row_vec(idx));
	assert(row->get_length() == mat->get_num_cols());
	for (size_t i = 0; i < row->get_length(); i++)
		assert(row->get<long>(i) == mat->get_num_cols() * idx + i);

	mat = detail::EM_matrix_store::create(9999, 10, matrix_layout_t::L_ROW,
			get_scalar_type<long>());
	mat->set_data(set_row_operate(mat->get_num_cols()));
	row = detail::smp_vec_store::cast(mat->get_row_vec(idx));
	assert(row->get_length() == mat->get_num_cols());
	for (size_t i = 0; i < row->get_length(); i++)
		assert(row->get<long>(i) == mat->get_num_cols() * idx + i);
}

void _test_stream(size_t num_rows, size_t num_cols, matrix_layout_t layout)
{
	printf("stream to EM matrix (%ld,%ld) layout: %d\n", num_rows, num_cols,
			layout);
	detail::EM_matrix_store::ptr mat = detail::EM_matrix_store::create(
			num_rows, num_cols, layout, get_scalar_type<long>());
	detail::mem_matrix_store::ptr mem_mat = detail::mem_matrix_store::create(
			num_rows, num_cols, layout, mat->get_type(), -1);
	for (size_t part_size = detail::mem_matrix_store::CHUNK_SIZE / 2;
			part_size <= detail::mem_matrix_store::CHUNK_SIZE / 2; part_size *= 2) {
		detail::EM_matrix_stream::ptr stream = detail::EM_matrix_stream::create(mat);
		size_t num_parts;
		if (mat->is_wide())
			num_parts = div_ceil(mat->get_num_cols(), part_size);
		else
			num_parts = div_ceil(mat->get_num_rows(), part_size);
		std::vector<size_t> part_ids(num_parts);
		for (size_t i = 0; i < part_ids.size(); i++)
			part_ids[i] = i;
		while (!part_ids.empty()) {
			size_t idx = random() % std::min(40UL, part_ids.size());
			size_t part_id = part_ids[idx];
			part_ids.erase(part_ids.begin() + idx);
			size_t actual_part_size;
			if (mat->is_wide())
				actual_part_size = std::min(part_size,
						mat->get_num_cols() - part_id * part_size);
			else
				actual_part_size = std::min(part_size,
						mat->get_num_rows() - part_id * part_size);
			detail::local_matrix_store::ptr buf;
			if (mat->is_wide() && mat->store_layout() == matrix_layout_t::L_COL) {
				buf = detail::local_matrix_store::ptr(
						new detail::local_buf_col_matrix_store(
							0, part_id * part_size, mat->get_num_rows(),
							actual_part_size, mat->get_type(), -1));
				buf->set_data(set_col_operate(mat->get_num_cols()));
			}
			else if (mat->is_wide()) {
				buf = detail::local_matrix_store::ptr(
						new detail::local_buf_row_matrix_store(
							0, part_id * part_size, mat->get_num_rows(),
							actual_part_size, mat->get_type(), -1));
				buf->set_data(set_row_operate(mat->get_num_cols()));
			}
			else if (mat->store_layout() == matrix_layout_t::L_COL) {
				buf = detail::local_matrix_store::ptr(
						new detail::local_buf_col_matrix_store(
							part_id * part_size, 0, actual_part_size,
							mat->get_num_cols(), mat->get_type(), -1));
				buf->set_data(set_col_operate(mat->get_num_cols()));
			}
			else {
				buf = detail::local_matrix_store::ptr(
						new detail::local_buf_row_matrix_store(
							part_id * part_size, 0, actual_part_size,
							mat->get_num_cols(), mat->get_type(), -1));
				buf->set_data(set_row_operate(mat->get_num_cols()));
			}
			mem_mat->write_portion_async(buf, buf->get_global_start_row(),
					buf->get_global_start_col());
			stream->write_async(buf, buf->get_global_start_row(),
					buf->get_global_start_col());
		}
		stream->flush();
		mat->wait4complete();

		dense_matrix::ptr tmp1 = dense_matrix::create(mat);
		dense_matrix::ptr tmp2 = dense_matrix::create(mem_mat);
		scalar_variable::ptr sum1 = tmp1->sum();
		scalar_variable::ptr sum2 = tmp2->sum();
		dense_matrix::ptr diff = tmp1->minus(*tmp2)->abs();
		scalar_variable::ptr sum3 = diff->sum();
		printf("sum1: %ld, sum2: %ld\n", scalar_variable::get_val<long>(*sum1),
				scalar_variable::get_val<long>(*sum2));
		printf("sum of diff: %ld\n", scalar_variable::get_val<long>(*sum3));
		assert(scalar_variable::get_val<long>(*sum3) == 0);
	}
}

void test_stream()
{
	_test_stream(9999999, 10, matrix_layout_t::L_ROW);
	_test_stream(10, 9999999, matrix_layout_t::L_COL);
	_test_stream(10, 9999999, matrix_layout_t::L_ROW);
	_test_stream(9999999, 10, matrix_layout_t::L_COL);
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

	test_stream();
	test_get_row_col();
	test_set_data();
	test_get_portion();

	destroy_flash_matrix();
}
