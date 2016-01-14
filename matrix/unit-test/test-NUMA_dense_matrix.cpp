#include <stdlib.h>

#include "NUMA_dense_matrix.h"
#include "mem_worker_thread.h"
#include "local_matrix_store.h"
#include "sparse_matrix.h"

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
};

class set_rand: public type_set_operate<long>
{
public:
	void set(long *arr, size_t num_eles, off_t row_idx, off_t col_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random();
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

void test_write2file1(NUMA_matrix_store::ptr mat)
{
	char *tmp_file_name = tempnam(".", "tmp.mat");
	if (mat->store_layout() == matrix_layout_t::L_ROW)
		mat->set_data(set_row_operate(mat->get_num_cols()));
	else
		mat->set_data(set_col_operate(mat->get_num_cols()));
	bool ret = mat->write2file(tmp_file_name);
	assert(ret);

	mem_matrix_store::ptr read_mat = mem_matrix_store::load(tmp_file_name);
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

void test_write2file()
{
	NUMA_matrix_store::ptr mat;
	printf("write a tall row matrix\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_ROW, get_scalar_type<long>());
	test_write2file1(mat);
	printf("write a tall column matrix\n");
	mat = NUMA_matrix_store::create(10000000, 10, num_nodes,
			matrix_layout_t::L_COL, get_scalar_type<long>());
	test_write2file1(mat);
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

	test_write2file();
	test_portion();
	test_init();
	test_transpose();

	destroy_flash_matrix();
}
