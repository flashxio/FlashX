#include <stdlib.h>

#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"
#include "local_matrix_store.h"
#include "sparse_matrix.h"
#include "data_io.h"

using namespace fm;
using namespace fm::detail;

int num_nodes = 1;
int num_threads = 8;

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
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

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
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

class set_rand: public type_set_operate<long>
{
public:
	void set(long *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random();
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

void test_init1(NUMA_matrix_store &store)
{
	store.reset_data();
#pragma omp parallel for
	for (size_t i = 0; i < store.get_num_rows(); i++) {
		for (size_t j = 0; j < store.get_num_cols(); j++)
			assert(store.get<long>(i, j) == 0);
	}

	if (store.store_layout() == matrix_layout_t::L_ROW)
		store.set_data(set_row_operate(store.get_num_cols()));
	else
		store.set_data(set_col_operate(store.get_num_cols()));
#pragma omp parallel for
	for (size_t i = 0; i < store.get_num_rows(); i++) {
		for (size_t j = 0; j < store.get_num_cols(); j++)
			assert(store.get<long>(i, j) == i * store.get_num_cols() + j);
	}
}

void test_init()
{
	NUMA_matrix_store::ptr mat;

	printf("test init on row tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_init1(*mat);

	printf("test init on col tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_init1(*mat);

	printf("test init on row wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_init1(*mat);

	printf("test init on col wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_init1(*mat);
}

void test_portion1(NUMA_matrix_store &store)
{
	std::pair<size_t, size_t> portion_size = store.get_portion_size();
	if (store.store_layout() == matrix_layout_t::L_ROW)
		store.set_data(set_row_operate(store.get_num_cols()));
	else
		store.set_data(set_col_operate(store.get_num_cols()));
#pragma omp parallel for
	for (size_t k = 0; k < store.get_num_portions(); k++) {
		local_matrix_store::ptr portion = store.get_portion(k);
		if (store.is_wide()) {
			assert(portion->get_global_start_row() == 0);
			assert(portion->get_global_start_col() == k * portion_size.second);
		}
		else {
			assert(portion->get_global_start_row() == k * portion_size.first);
			assert(portion->get_global_start_col() == 0);
		}
		for (size_t i = 0; i < portion->get_num_rows(); i++)
			for (size_t j = 0; j < portion->get_num_cols(); j++) {
				long expected = (portion->get_global_start_row()
						+ i) * store.get_num_cols()
					+ portion->get_global_start_col() + j;
				assert(portion->get<long>(i, j) == expected);
			}
	}
}

void test_portion()
{
	NUMA_matrix_store::ptr mat;

	printf("test get portion on row tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_portion1(*mat);

	printf("test get portion on col tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_portion1(*mat);

	printf("test get portion on row wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_portion1(*mat);

	printf("test get portion on col wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_portion1(*mat);
}

void test_transpose1(const NUMA_matrix_store &m1, const NUMA_matrix_store &m2)
{
	assert(m1.get_num_rows() == m2.get_num_cols());
	assert(m1.get_num_cols() == m2.get_num_rows());
#pragma omp parallel for
	for (size_t i = 0; i < m1.get_num_rows(); i++)
		for (size_t j = 0; j < m1.get_num_cols(); j++)
			assert(m1.get<int>(i, j) == m2.get<int>(j, i));
}

void test_transpose()
{
	// Test on local buffer matrix.
	NUMA_matrix_store::ptr mat;

	printf("test transpose on row tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_transpose1(*mat, *NUMA_matrix_store::cast(mat->transpose()));

	printf("test transpose on col tall\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_transpose1(*mat, *NUMA_matrix_store::cast(mat->transpose()));

	printf("test transpose on row wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_transpose1(*mat, *NUMA_matrix_store::cast(mat->transpose()));

	printf("test transpose on col wide\n");
	mat = NUMA_matrix_store::create(10, 10000000, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_transpose1(*mat, *NUMA_matrix_store::cast(mat->transpose()));
}

void test_write2file1(mem_matrix_store::ptr mat)
{
	char *tmp_file_name = tempnam(".", "tmp.mat");
	if (mat->store_layout() == matrix_layout_t::L_ROW)
		mat->set_data(set_row_operate(mat->get_num_cols()));
	else
		mat->set_data(set_col_operate(mat->get_num_cols()));
	bool ret = mat->write2file(tmp_file_name);
	assert(ret);

	mem_matrix_store::const_ptr read_mat = mem_matrix_store::load(tmp_file_name,
			mat->get_num_nodes());
	assert(read_mat);
	assert(read_mat->get_num_rows() == mat->get_num_rows());
	assert(read_mat->get_num_cols() == mat->get_num_cols());
	assert(read_mat->get_type() == mat->get_type());
	assert(read_mat->store_layout() == mat->store_layout());
#pragma omp parallel for
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++)
			assert(mat->get<long>(i, j) == read_mat->get<long>(i, j));
	}

	unlink(tmp_file_name);
}

void test_write2text_file(mem_matrix_store::ptr mat)
{
	char *tmp_file_name = tempnam(".", "tmp_mat.txt");
	if (mat->store_layout() == matrix_layout_t::L_ROW)
		mat->set_data(set_row_operate(mat->get_num_cols()));
	else
		mat->set_data(set_col_operate(mat->get_num_cols()));
	bool ret = mat->write2file(tmp_file_name, true);
	assert(ret);

	printf("read %s\n", tmp_file_name);
	dense_matrix::ptr read_mat = read_matrix(
			std::vector<std::string>(1, tmp_file_name),
			true, true, "L", " ", mat->get_num_cols());
	assert(read_mat);
	read_mat = read_mat->conv2(matrix_layout_t::L_ROW);
	read_mat->materialize_self();
	printf("read mat: %ld, %ld\n", read_mat->get_num_rows(),
			read_mat->get_num_cols());
	assert(read_mat->get_num_rows() == mat->get_num_rows());
	assert(read_mat->get_num_cols() == mat->get_num_cols());
	assert(read_mat->get_type() == mat->get_type());
	assert(read_mat->store_layout() == mat->store_layout());
	mem_matrix_store::const_ptr read_store
		= std::dynamic_pointer_cast<const mem_matrix_store>(
				read_mat->get_raw_store());
#pragma omp parallel for
	for (size_t i = 0; i < mat->get_num_rows(); i++) {
		for (size_t j = 0; j < mat->get_num_cols(); j++) {
			if (mat->get<long>(i, j) != read_store->get<long>(i, j))
				printf("%ld,%ld: %ld, %ld\n", i, j, mat->get<long>(i, j),
						read_store->get<long>(i, j));
			assert(mat->get<long>(i, j) == read_store->get<long>(i, j));
		}
	}

	unlink(tmp_file_name);
}

void test_write2file()
{
	mem_matrix_store::ptr mat;

	printf("write a tall row matrix to text\n");
	mat = mem_matrix_store::create(1000000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_write2text_file(mat);

	printf("write a tall row matrix\n");
	mat = mem_matrix_store::create(10000000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_write2file1(mat);
	printf("write a tall column matrix\n");
	mat = mem_matrix_store::create(10000000, 10,
			matrix_layout_t::L_COL, get_scalar_type<long>(), num_nodes);
	test_write2file1(mat);

	printf("write a wide row matrix\n");
	mat = mem_matrix_store::create(10, 10000000,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_write2file1(mat);
	printf("write a wide column matrix\n");
	mat = mem_matrix_store::create(10, 10000000,
			matrix_layout_t::L_COL, get_scalar_type<long>(), num_nodes);
	test_write2file1(mat);
}

void test_resize(mem_matrix_store::ptr mat1)
{
	mat1->set_data(set_col_operate(mat1->get_num_cols()));
	mem_matrix_store::ptr mat2 = mem_matrix_store::create(mat1->get_num_rows(),
			mat1->get_num_cols(), mat1->store_layout(), mat1->get_type(),
			mat1->get_num_nodes());
	mat2->set_data(set_col_operate(mat1->get_num_cols()));
	if (!mat1->is_wide()) {
		size_t new_num_rows = random() % mat2->get_num_rows();
		bool ret = mat2->resize(new_num_rows, mat2->get_num_cols());
		assert(ret);
		assert(mat2->get_num_rows() == new_num_rows);
	}
	else {
		size_t new_num_cols = random() % mat2->get_num_cols();
		bool ret = mat2->resize(mat2->get_num_rows(), new_num_cols);
		assert(ret);
		assert(mat2->get_num_cols() == new_num_cols);
	}
	for (size_t i = 0; i < mat2->get_num_rows(); i++)
		for (size_t j = 0; j < mat2->get_num_cols(); j++)
			assert(mat1->get<long>(i, j) == mat2->get<long>(i, j));
}

void test_resize()
{
	printf("test resize\n");
	mem_matrix_store::ptr mat1 = mem_matrix_store::create(100000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(100000, 10, matrix_layout_t::L_COL,
			get_scalar_type<long>(), num_nodes);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(10, 100000, matrix_layout_t::L_ROW,
			get_scalar_type<long>(), num_nodes);
	test_resize(mat1);

	mat1 = mem_matrix_store::create(10, 100000, matrix_layout_t::L_COL,
			get_scalar_type<long>(), num_nodes);
	test_resize(mat1);
}

void test_tall_write(mem_matrix_store::ptr store)
{
	size_t portion_size = store->get_portion_size().first * 2;
	for (size_t start_row = 0; start_row < store->get_num_rows(); ) {
		size_t num_rows = (random() % portion_size) + 100;
		num_rows = std::min(num_rows, store->get_num_rows() - start_row);
		local_matrix_store::ptr lstore;
		if (store->store_layout() == matrix_layout_t::L_ROW) {
			lstore = local_matrix_store::ptr(new local_buf_row_matrix_store(
						start_row, 0, num_rows, store->get_num_cols(),
						store->get_type(), -1, false));
			lstore->set_data(set_row_operate(store->get_num_cols()));
		}
		else {
			lstore = local_matrix_store::ptr(new local_buf_col_matrix_store(
						start_row, 0, num_rows, store->get_num_cols(),
						store->get_type(), -1, false));
			lstore->set_data(set_col_operate(store->get_num_cols()));
		}
		store->write_portion_async(lstore, lstore->get_global_start_row(),
				lstore->get_global_start_col());
		start_row += num_rows;
	}

	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++) {
			assert(store->get<long>(i, j) == i * store->get_num_cols() + j);
		}
}

void test_wide_write(mem_matrix_store::ptr store)
{
	size_t portion_size = store->get_portion_size().second * 2;
	for (size_t start_col = 0; start_col < store->get_num_cols(); ) {
		size_t num_cols = (random() % portion_size) + 100;
		num_cols = std::min(num_cols, store->get_num_cols() - start_col);
		local_matrix_store::ptr lstore;
		if (store->store_layout() == matrix_layout_t::L_ROW) {
			lstore = local_matrix_store::ptr(new local_buf_row_matrix_store(
						0, start_col, store->get_num_rows(), num_cols,
						store->get_type(), -1, false));
			lstore->set_data(set_row_operate(store->get_num_cols()));
		}
		else {
			lstore = local_matrix_store::ptr(new local_buf_col_matrix_store(
						0, start_col, store->get_num_rows(), num_cols,
						store->get_type(), -1, false));
			lstore->set_data(set_col_operate(store->get_num_cols()));
		}
		store->write_portion_async(lstore, lstore->get_global_start_row(),
				lstore->get_global_start_col());
		start_col += num_cols;
	}

	for (size_t i = 0; i < store->get_num_rows(); i++)
		for (size_t j = 0; j < store->get_num_cols(); j++) {
			if (store->get<long>(i, j) != i * store->get_num_cols() + j)
				printf("%ld,%ld: %ld,%ld\n", i, j, store->get<long>(i, j),
						i * store->get_num_cols() + j);
			assert(store->get<long>(i, j) == i * store->get_num_cols() + j);
		}
}

void test_write_portion()
{
	printf("test write portion\n");
	mem_matrix_store::ptr mat1 = mem_matrix_store::create(100000, 10,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_tall_write(mat1);

	mat1 = mem_matrix_store::create(10, 100000,
			matrix_layout_t::L_ROW, get_scalar_type<long>(), num_nodes);
	test_wide_write(mat1);
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
	num_nodes = safs::params.get_num_nodes();

	test_write_portion();
	test_write2file();
	test_resize();
	test_portion();
	test_init();
	test_transpose();

	destroy_flash_matrix();
}
