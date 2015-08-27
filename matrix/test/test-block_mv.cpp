#ifdef PROFILER
#include <google/profiler.h>
#endif

#include "eigensolver/block_dense_matrix.h"

using namespace fm;

size_t long_dim = 60 * 1024 * 1024;
size_t repeats = 1;

void test_gemm(eigen::block_multi_vector::ptr mv)
{
	bool in_mem = mv->get_block(0)->is_in_mem();
	int num_nodes = mv->get_block(0)->get_data().get_num_nodes();
	dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
			mv->get_num_cols(), mv->get_block_size(),
			matrix_layout_t::L_COL, -1, true);

	detail::mem_col_matrix_store::const_ptr B
		= detail::mem_col_matrix_store::cast(mat->get_raw_store());
	scalar_variable_impl<double> alpha(2);
	scalar_variable_impl<double> beta(0);

	struct timeval start, end;

#if 0
	eigen::block_multi_vector::ptr res0 = eigen::block_multi_vector::create(
			long_dim, mv->get_block_size(), mv->get_block_size(),
			get_scalar_type<double>(), true);
	res0->set_block(0, dense_matrix::create_const<double>(0, long_dim,
				mv->get_block_size(), matrix_layout_t::L_COL, num_nodes, in_mem));
	res0->set_multiply_blocks(1);
	gettimeofday(&start, NULL);
	res0 = res0->gemm(*mv, B, alpha, beta);
	assert(res0->get_num_blocks() == 1);
	dense_matrix::ptr res_mat0 = res0->get_block(0);
	res_mat0->materialize_self();
	gettimeofday(&end, NULL);
	printf("agg materialization takes %.3f seconds\n", time_diff(start, end));
#endif

	for (size_t i = 0; i < repeats; i++) {
		eigen::block_multi_vector::ptr res1 = eigen::block_multi_vector::create(
				long_dim, mv->get_block_size(), mv->get_block_size(),
				get_scalar_type<double>(), true);
		res1->set_block(0, dense_matrix::create_const<double>(0, long_dim,
					mv->get_block_size(), matrix_layout_t::L_COL, num_nodes, in_mem));
		res1->set_multiply_blocks(4);
		gettimeofday(&start, NULL);
		res1 = res1->gemm(*mv, B, alpha, beta);
		assert(res1->get_num_blocks() == 1);
		dense_matrix::ptr res_mat1 = res1->get_block(0);
		res_mat1->materialize_self();
		gettimeofday(&end, NULL);
		printf("hierarchical materialization takes %.3f seconds\n",
				time_diff(start, end));
	}

	for (size_t i = 0; i < repeats; i++) {
		eigen::block_multi_vector::ptr res2 = eigen::block_multi_vector::create(
				long_dim, mv->get_block_size(), mv->get_block_size(),
				get_scalar_type<double>(), true);
		res2->set_block(0, dense_matrix::create_const<double>(0, long_dim,
					mv->get_block_size(), matrix_layout_t::L_COL, num_nodes, in_mem));
		res2->set_multiply_blocks(mv->get_num_blocks());
		gettimeofday(&start, NULL);
		res2 = res2->gemm(*mv, B, alpha, beta);
		assert(res2->get_num_blocks() == 1);
		dense_matrix::ptr res_mat2 = res2->get_block(0);
		res_mat2->materialize_self();
		gettimeofday(&end, NULL);
		printf("flat materialization takes %.3f seconds\n", time_diff(start, end));
	}

#if 0
	dense_matrix::ptr diff = res_mat1->minus(*res_mat2);
	scalar_variable::ptr max_diff = diff->abs()->max();
	scalar_variable::ptr max1 = res_mat1->max();
	scalar_variable::ptr max2 = res_mat2->max();
	printf("max diff: %g, max mat1: %g, max mat2: %g\n",
			*(double *) max_diff->get_raw(), *(double *) max1->get_raw(),
			*(double *) max2->get_raw());
#endif
}

void test_gemm(bool in_mem, size_t block_size, size_t min_num_blocks,
		size_t max_num_blocks, size_t num_cached_blocks)
{
	std::vector<dense_matrix::ptr> mats(max_num_blocks);
	for (size_t i = 0; i < mats.size(); i++) {
		if (i < num_cached_blocks)
			mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
					block_size, matrix_layout_t::L_COL,
					matrix_conf.get_num_nodes(), true);
		else
			mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
					block_size, matrix_layout_t::L_COL,
					matrix_conf.get_num_nodes(), in_mem);
	}

	for (size_t num_blocks = min_num_blocks; num_blocks <= max_num_blocks;
			num_blocks *= 2) {
		eigen::block_multi_vector::ptr mv = eigen::block_multi_vector::create(
				long_dim, num_blocks * block_size, block_size,
				get_scalar_type<double>(), in_mem);
		printf("gemm on block multi-vector (block size: %ld, #blocks: %ld)\n",
				mv->get_block_size(), mv->get_num_blocks());
		for (size_t i = 0; i < mv->get_num_blocks(); i++)
			mv->set_block(i, mats[i]);
		test_gemm(mv);
	}
}

void test_MvTransMv(eigen::block_multi_vector::ptr mv1,
		eigen::block_multi_vector::ptr mv2)
{
	struct timeval start, end;

#if 0
	mv1->set_multiply_blocks(1);
	gettimeofday(&start, NULL);
	fm::dense_matrix::ptr res1 = mv1->MvTransMv(*mv2);
	gettimeofday(&end, NULL);
	printf("MvTransMv (1 block) takes %.3f seconds\n", time_diff(start, end));
#endif

	mv1->set_multiply_blocks(4);
#ifdef PROFILER
	ProfilerStart("MvTransMv.4B.prof");
#endif
	for (size_t i = 0; i < repeats; i++) {
		gettimeofday(&start, NULL);
		fm::dense_matrix::ptr res2 = mv1->MvTransMv(*mv2);
		gettimeofday(&end, NULL);
		printf("MvTransMv (4 blocks) takes %.3f seconds\n", time_diff(start, end));
	}
#ifdef PROFILER
	ProfilerStop();
#endif

	mv1->set_multiply_blocks(mv1->get_num_blocks());
#ifdef PROFILER
	ProfilerStart("MvTransMv.all.prof");
#endif
	for (size_t i = 0; i < repeats; i++) {
		gettimeofday(&start, NULL);
		fm::dense_matrix::ptr res3 = mv1->MvTransMv(*mv2);
		gettimeofday(&end, NULL);
		printf("MvTransMv (all blocks) takes %.3f seconds\n", time_diff(start, end));
	}
#ifdef PROFILER
	ProfilerStop();
#endif

#if 0
	dense_matrix::ptr diff = res2->minus(*res3);
	scalar_variable::ptr max_diff = diff->abs()->max();
	scalar_variable::ptr max1 = res2->max();
	scalar_variable::ptr max2 = res3->max();
	printf("max diff: %g, max mat1: %g, max mat2: %g\n",
			*(double *) max_diff->get_raw(), *(double *) max1->get_raw(),
			*(double *) max2->get_raw());
#endif
}

void test_MvTransMv(bool in_mem, size_t block_size,
		size_t min_num_blocks, size_t max_num_blocks, size_t num_cached_blocks)
{
	std::vector<dense_matrix::ptr> mats(max_num_blocks);
	for (size_t i = 0; i < mats.size(); i++) {
		if (i < num_cached_blocks)
			mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
					block_size, matrix_layout_t::L_COL,
					matrix_conf.get_num_nodes(), true);
		else
			mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
					block_size, matrix_layout_t::L_COL,
					matrix_conf.get_num_nodes(), in_mem);
	}
	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, block_size, block_size, get_scalar_type<double>(), in_mem);
	mv2->set_block(0, dense_matrix::create_randu<double>(0, 1, long_dim,
				mv2->get_block_size(), matrix_layout_t::L_COL,
				matrix_conf.get_num_nodes(), in_mem));

	for (size_t num_blocks = min_num_blocks; num_blocks <= max_num_blocks;
			num_blocks *= 2) {
		eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
				long_dim, num_blocks * block_size, block_size,
				get_scalar_type<double>(), in_mem);
		for (size_t i = 0; i < mv1->get_num_blocks(); i++)
			mv1->set_block(i, mats[i]);
		printf("MvTransMv on block MV (block size: %ld, #blocks: %ld)\n",
				block_size, mv1->get_num_blocks());
		test_MvTransMv(mv1, mv2);
	}
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

	test_gemm(true, 4, 1, 128, 0);
	test_gemm(true, 64, 1, 8, 0);
	test_MvTransMv(true, 4, 1, 128, 0);
	test_MvTransMv(true, 64, 1, 8, 0);

	test_gemm(false, 4, 8, 128, 4);
	test_gemm(false, 64, 1, 8, 0);
	test_MvTransMv(false, 4, 8, 128, 4);
	test_MvTransMv(false, 64, 1, 8, 0);

	destroy_flash_matrix();
}
