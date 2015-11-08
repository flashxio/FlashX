#include "NUMA_vector.h"
#include "NUMA_dense_matrix.h"
#include "matrix_config.h"
#include "mem_worker_thread.h"
#include "matrix_store.h"

const size_t num_eles = 1024 * 1024 * 10;
size_t num_nodes = 1;
size_t nthreads = 8;

using namespace fm;
using namespace detail;

class set_seq_vec: public type_set_vec_operate<long>
{
public:
	virtual void set(long *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = start_idx + i;
	}
};

class set_rand_vec: public type_set_vec_operate<long>
{
public:
	virtual void set(long *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random();
	}
};

void test_init()
{
	printf("test init NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	vec->set_data(set_seq_vec());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert((size_t) vec->get<long>(i) == i);
}

void test_mapping()
{
	printf("test mapping in NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::vector<size_t> lens = vec->get_mapper().cal_local_lengths(
			vec->get_length());
	size_t tot_len = 0;
	for (size_t i = 0; i < lens.size(); i++)
		tot_len += lens[i];
	assert(tot_len == vec->get_length());

	for (size_t i = 0; i < vec->get_length(); i++) {
		std::pair<int, size_t> phy_loc = vec->get_mapper().map2physical(i);
		size_t loc = vec->get_mapper().map2logical(phy_loc.first, phy_loc.second);
		assert(loc == i);
	}
}

void test_copy()
{
	printf("test copying NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::unique_ptr<long[]> raw_arr(new long[vec->get_length()]);
	for (size_t i = 0; i < vec->get_length(); i++)
		raw_arr[i] = random();
	vec->copy_from((char *) raw_arr.get(), vec->get_length() * sizeof(long));
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(raw_arr[i] == vec->get<long>(i));

	NUMA_vec_store::ptr vec1 = NUMA_vec_store::create(vec->get_length(), num_nodes,
			get_scalar_type<long>());
	bool ret = vec1->copy_from(*vec);
	assert(ret);
	assert(vec1->get_length() == vec->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(vec1->get<long>(i) == vec->get<long>(i));
}

void test_deep_copy()
{
	printf("test deep copying NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	vec->set_data(set_seq_vec());

	NUMA_vec_store::ptr vec1 = NUMA_vec_store::cast(vec->deep_copy());
	NUMA_vec_store::ptr vec2 = NUMA_vec_store::cast(vec->deep_copy());
	vec->set_data(set_rand_vec());
	assert(vec1->get_length() == vec2->get_length());
	for (size_t i = 0; i < vec1->get_length(); i++)
		assert(vec1->get<long>(i) == vec2->get<long>(i));
}

void test_sort()
{
	printf("test sorting NUMA vector store\n");
	NUMA_vec_store::ptr vec = NUMA_vec_store::create((random() % num_eles) + num_eles,
			num_nodes, get_scalar_type<long>());
	std::unique_ptr<long[]> raw_arr(new long[vec->get_length()]);
	for (size_t i = 0; i < vec->get_length(); i++)
		raw_arr[i] = random();
	vec->copy_from((char *) raw_arr.get(), vec->get_length() * sizeof(long));
	vec->sort();
	std::sort(raw_arr.get(), raw_arr.get() + vec->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(raw_arr[i] == vec->get<long>(i));
}

void test_conv2mat()
{
	printf("test converting a NUMA vector to a matrix\n");
	size_t nrow = (random() % num_eles) + num_eles;
	size_t ncol = 3;
	NUMA_vec_store::ptr vec = NUMA_vec_store::create(nrow * ncol, num_nodes,
			get_scalar_type<long>());
	matrix_store::const_ptr mat = vec->conv2mat(nrow * ncol, 1, true);
	assert(mat->get_num_rows() == nrow * ncol);
	assert(mat->get_num_cols() == 1);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);

	mat = vec->conv2mat(nrow * ncol, 1, false);
	assert(mat->get_num_rows() == nrow * ncol);
	assert(mat->get_num_cols() == 1);
	assert(mat->store_layout() == matrix_layout_t::L_COL);

	mat = vec->conv2mat(1, nrow * ncol, true);
	assert(mat->get_num_rows() == 1);
	assert(mat->get_num_cols() == nrow * ncol);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);

	mat = vec->conv2mat(1, nrow * ncol, false);
	assert(mat->get_num_rows() == 1);
	assert(mat->get_num_cols() == nrow * ncol);
	assert(mat->store_layout() == matrix_layout_t::L_COL);

	mat = vec->conv2mat(nrow, ncol, true);
	assert(mat->get_num_rows() == nrow);
	assert(mat->get_num_cols() == ncol);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);
	{
		NUMA_row_tall_matrix_store::const_ptr numa_mat
			= std::dynamic_pointer_cast<const NUMA_row_tall_matrix_store>(mat);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				assert(*(long *) numa_mat->get(i, j) == vec->get<long>(
							i * ncol + j));
	}

	mat = vec->conv2mat(nrow, ncol, false);
	assert(mat->get_num_rows() == nrow);
	assert(mat->get_num_cols() == ncol);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
	{
		NUMA_col_tall_matrix_store::const_ptr numa_mat
			= std::dynamic_pointer_cast<const NUMA_col_tall_matrix_store>(mat);
		for (size_t i = 0; i < ncol; i++)
			for (size_t j = 0; j < nrow; j++)
				assert(*(long *) numa_mat->get(j, i) == vec->get<long>(
							i * nrow + j));
	}

	mat = vec->conv2mat(ncol, nrow, true);
	assert(mat->get_num_rows() == ncol);
	assert(mat->get_num_cols() == nrow);
	assert(mat->store_layout() == matrix_layout_t::L_ROW);
	{
		NUMA_row_wide_matrix_store::const_ptr numa_mat
			= std::dynamic_pointer_cast<const NUMA_row_wide_matrix_store>(mat);
		for (size_t i = 0; i < ncol; i++)
			for (size_t j = 0; j < nrow; j++)
				assert(*(long *) numa_mat->get(i, j) == vec->get<long>(
							i * nrow + j));
	}

	mat = vec->conv2mat(ncol, nrow, false);
	assert(mat->get_num_rows() == ncol);
	assert(mat->get_num_cols() == nrow);
	assert(mat->store_layout() == matrix_layout_t::L_COL);
	{
		NUMA_col_wide_matrix_store::const_ptr numa_mat
			= std::dynamic_pointer_cast<const NUMA_col_wide_matrix_store>(mat);
		for (size_t i = 0; i < nrow; i++)
			for (size_t j = 0; j < ncol; j++)
				assert(*(long *) numa_mat->get(j, i) == vec->get<long>(
							i * ncol + j));
	}
}

int main(int argc, char *argv[])
{
	if (argc >= 3) {
		num_nodes = atoi(argv[1]);
		nthreads = atoi(argv[2]);
	}

	matrix_conf.set_num_nodes(num_nodes);
	matrix_conf.set_num_DM_threads(nthreads);
	detail::mem_thread_pool::init_global_mem_threads(num_nodes,
			nthreads / num_nodes);
	test_init();
	test_mapping();
	test_copy();
	test_deep_copy();
	test_sort();
	test_conv2mat();
}
