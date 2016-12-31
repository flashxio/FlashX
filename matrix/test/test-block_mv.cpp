#ifdef PROFILER
#include <google/profiler.h>
#endif

#include "RAID_config.h"

#include "EM_dense_matrix.h"
#include "eigensolver/block_dense_matrix.h"
#include "matrix_stats.h"

using namespace fm;

size_t long_dim = 60 * 1024 * 1024;
size_t repeats = 1;
size_t multiply_blocks = 8;

safs::safs_file_group::ptr group;

class rand_rotate_file_group: public safs::safs_file_group
{
	// The base permute is the permutation that other permutations are based
	// on. Other permutations just rotate the base permutation from a random
	// location. Every #disks files share the same base permutation.
	std::vector<std::vector<int> > base_permutes;
	std::vector<std::vector<int> > rand_rotates;
	size_t num_files;
	size_t num_same;
public:
	rand_rotate_file_group(const safs::RAID_config &conf, size_t num_same);
	std::vector<int> add_file(safs::safs_file &file);
	std::string get_name() const {
		return std::string("rand_rotate-") + itoa(num_same);
	}
};

static std::vector<int> shuffle_disks(int num_disks)
{
	std::vector<int> permute(num_disks);
	for (size_t i = 0; i < permute.size(); i++)
		permute[i] = i;
	random_shuffle(permute.begin(), permute.end());
	return permute;
}

rand_rotate_file_group::rand_rotate_file_group(const safs::RAID_config &conf,
		size_t num_same)
{
	this->num_same = num_same;
	int num_disks = conf.get_num_disks();
	assert(num_disks % num_same == 0);
	num_files = 0;
	base_permutes.push_back(shuffle_disks(num_disks));
	std::vector<int> base_rotate = shuffle_disks(num_disks);
	base_rotate.resize(num_disks / num_same);
	std::vector<int> rotate(num_disks);
	for (size_t i = 0; i < base_rotate.size(); i++) {
		for (size_t j = 0; j < num_same; j++) {
			rotate[i * num_same + j] = base_rotate[i];
			printf("rotate %ld: %d\n", i * num_same + j, base_rotate[i]);
		}
	}
	rand_rotates.push_back(rotate);
}

std::vector<int> rand_rotate_file_group::add_file(safs::safs_file &file)
{
	size_t num_disks = base_permutes.front().size();
	size_t base_idx = num_files / num_disks;
	if (base_idx >= base_permutes.size()) {
		base_permutes.push_back(shuffle_disks(num_disks));
		std::vector<int> base_rotate = shuffle_disks(num_disks);
		base_rotate.resize(num_disks / num_same);
		std::vector<int> rotate(num_disks);
		for (size_t i = 0; i < base_rotate.size(); i++) {
			for (size_t j = 0; j < num_same; j++) {
				rotate[i * num_same + j] = base_rotate[i];
				printf("rotate %ld: %d\n", i * num_same + j, base_rotate[i]);
			}
		}
		rand_rotates.push_back(rotate);
	}
	assert(base_permutes.size() > base_idx);

	std::vector<int> base = base_permutes[base_idx];
	std::vector<int> ret(num_disks);
	size_t rotate = rand_rotates[base_idx][num_files % num_disks];
	for (size_t i = 0; i < ret.size(); i++)
		ret[i] = base[(rotate + i) % num_disks];
	num_files++;
	return ret;
}

class rand_permute_file_group: public safs::safs_file_group
{
	size_t num_files;
	size_t num_disks;
	size_t num_same;
	std::vector<int> permute;
public:
	rand_permute_file_group(const safs::RAID_config &conf, size_t num_same) {
		this->num_same = num_same;
		num_disks = conf.get_num_disks();
		assert(num_disks % num_same == 0);
		num_files = 0;
	}

	std::vector<int> add_file(safs::safs_file &file) {
		if (num_files % num_same == 0)
			permute = shuffle_disks(num_disks);
		num_files++;
		return permute;
	}

	std::string get_name() const {
		return std::string("rand_permute-") + itoa(num_same);
	}
};

std::vector<dense_matrix::ptr> get_EM_matrices(size_t num_rows, size_t num_cols,
		size_t num_mats)
{
	std::vector<dense_matrix::ptr> mats(num_mats);
	for (size_t i = 0; i < mats.size(); i++) {
		std::string mat_name = (boost::format("test-%1%rows-%2%cols-%3%.mat")
				% num_rows % num_cols % i).str();
		detail::EM_matrix_store::const_ptr store
			= detail::EM_matrix_store::create(mat_name);
		if (store == NULL) {
			mats[i] = dense_matrix::create_randu<double>(0, 1, num_rows, num_cols,
					matrix_layout_t::L_COL, -1, false, group);
			detail::EM_matrix_store::const_ptr store
				= detail::EM_matrix_store::cast(mats[i]->get_raw_store());
			store->set_persistent(mat_name);
		}
		else
			mats[i] = dense_matrix::create(store);
	}
	return mats;
}

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
			get_scalar_type<double>(), true, false);
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
				get_scalar_type<double>(), true, false);
		res1->set_block(0, dense_matrix::create_const<double>(0, long_dim,
					mv->get_block_size(), matrix_layout_t::L_COL, num_nodes, in_mem));
		res1->set_multiply_blocks(multiply_blocks);
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
				get_scalar_type<double>(), true, false);
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
	std::vector<dense_matrix::ptr> EM_mats;
	if (!in_mem)
		EM_mats = get_EM_matrices(long_dim, block_size,
				mats.size() - num_cached_blocks);
	for (size_t i = 0; i < mats.size(); i++) {
		if (!in_mem && i >= num_cached_blocks)
			mats[i] = EM_mats[i - num_cached_blocks];
	}

	for (size_t num_blocks = min_num_blocks; num_blocks <= max_num_blocks;
			num_blocks *= 2) {
		eigen::block_multi_vector::ptr mv = eigen::block_multi_vector::create(
				long_dim, num_blocks * block_size, block_size,
				get_scalar_type<double>(), in_mem, false);
		printf("gemm on block multi-vector (block size: %ld, #blocks: %ld)\n",
				mv->get_block_size(), mv->get_num_blocks());
		for (size_t i = 0; i < mv->get_num_blocks(); i++) {
			if (mats[i] == NULL)
				mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
						block_size, matrix_layout_t::L_COL,
						matrix_conf.get_num_nodes(), true);
			mv->set_block(i, mats[i]);
		}
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

	mv1->set_multiply_blocks(multiply_blocks);
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
	std::vector<dense_matrix::ptr> EM_mats;
	if (!in_mem)
		EM_mats = get_EM_matrices(long_dim, block_size,
				mats.size() - num_cached_blocks + 1);
	for (size_t i = 0; i < mats.size(); i++) {
		if (!in_mem && i >= num_cached_blocks)
			mats[i] = EM_mats[i - num_cached_blocks];
	}
	eigen::block_multi_vector::ptr mv2 = eigen::block_multi_vector::create(
			long_dim, block_size, block_size, get_scalar_type<double>(),
			in_mem, false);
	if (in_mem)
		mv2->set_block(0, dense_matrix::create_randu<double>(0, 1, long_dim,
					block_size, matrix_layout_t::L_COL,
					matrix_conf.get_num_nodes(), true));
	else
		mv2->set_block(0, EM_mats.back());

	for (size_t num_blocks = min_num_blocks; num_blocks <= max_num_blocks;
			num_blocks *= 2) {
		eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
				long_dim, num_blocks * block_size, block_size,
				get_scalar_type<double>(), in_mem, false);
		printf("MvTransMv on block MV (block size: %ld, #blocks: %ld)\n",
				block_size, mv1->get_num_blocks());
		for (size_t i = 0; i < mv1->get_num_blocks(); i++) {
			if (mats[i] == NULL)
				mats[i] = dense_matrix::create_randu<double>(0, 1, long_dim,
						block_size, matrix_layout_t::L_COL,
						matrix_conf.get_num_nodes(), true);
			mv1->set_block(i, mats[i]);
		}
		test_MvTransMv(mv1, mv2);
	}
}

void test_gemm_simple(bool in_mem, size_t dim1, size_t dim2)
{
	dense_matrix::ptr mat1;
	if (in_mem)
		mat1 = dense_matrix::create_randu<double>(0, 1, long_dim, dim1,
				matrix_layout_t::L_COL, matrix_conf.get_num_nodes(), in_mem);
	else
		mat1 = get_EM_matrices(long_dim, dim1, 1).front();
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(0, 0, dim1,
			dim2, matrix_layout_t::L_COL, matrix_conf.get_num_nodes(), true);

	printf("simple gemm starts\n");
	struct timeval start, end;
	// Multiply on the entire matrix.
	gettimeofday(&start, NULL);

	detail::matrix_stats_t orig_stats = detail::matrix_stats;
	dense_matrix::ptr res = mat1->multiply(*mat2, matrix_layout_t::L_NONE);
	res->materialize_self();
	detail::matrix_stats.print_diff(orig_stats);

	gettimeofday(&end, NULL);
	printf("simple gemm takes %.3f seconds\n", time_diff(start, end));
}

void test_MvTransMv_simple(bool in_mem, size_t dim1, size_t dim2)
{
	dense_matrix::ptr mat1 = dense_matrix::create_randu<double>(0, 0, long_dim,
			dim1, matrix_layout_t::L_COL, matrix_conf.get_num_nodes(), in_mem);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(0, 0, long_dim,
			dim2, matrix_layout_t::L_COL, matrix_conf.get_num_nodes(), in_mem);

	struct timeval start, end;
	gettimeofday(&start, NULL);

	detail::matrix_stats_t orig_stats = detail::matrix_stats;
	mat1->transpose()->multiply(*mat2, matrix_layout_t::L_NONE);
	detail::matrix_stats.print_diff(orig_stats);

	gettimeofday(&end, NULL);
	printf("simple MvTransMv takes %.3f seconds\n", time_diff(start, end));
}

void test_gemm_simul(size_t block_size, size_t num_blocks)
{
	std::vector<dense_matrix::ptr> mats1 = get_EM_matrices(long_dim,
			block_size, num_blocks);
	dense_matrix::ptr mat2 = dense_matrix::create_randu<double>(0, 1,
			block_size * num_blocks, block_size, matrix_layout_t::L_COL,
			-1, true);
	detail::mem_col_matrix_store::const_ptr B
		= detail::mem_col_matrix_store::cast(mat2->get_raw_store());

	eigen::block_multi_vector::ptr mv = eigen::block_multi_vector::create(
			long_dim, block_size * num_blocks, block_size,
			get_scalar_type<double>(), false, false);
	for (size_t i = 0; i < mats1.size(); i++)
		mv->set_block(i, mats1[i]);

	eigen::block_multi_vector::ptr mv1 = eigen::block_multi_vector::create(
			long_dim, block_size, block_size,
			get_scalar_type<double>(), false, false);
	{
		dense_matrix::ptr mat = dense_matrix::create_randu<double>(0, 1,
				long_dim, block_size, matrix_layout_t::L_COL,
				matrix_conf.get_num_nodes(), true);

		detail::smp_vec_store::ptr vstore = detail::smp_vec_store::create(
				block_size, get_scalar_type<double>());
		for (size_t i = 0; i < vstore->get_length(); i++)
			vstore->set<double>(i, i);
		vector::ptr vec = vector::create(vstore);

		dense_matrix::ptr small_mat = dense_matrix::create_randu<double>(
				0, 1, block_size, block_size, matrix_layout_t::L_COL);

		mat = mat->scale_cols(col_vec::create(vec));
		mat = mat->multiply(*small_mat, matrix_layout_t::L_NONE);
		mat = mat->scale_cols(col_vec::create(vec));
		mv1->set_block(0, mat);
	}

	// Multiply on the merged matrix.
	struct timeval start, end;
	gettimeofday(&start, NULL);
	scalar_variable_impl<double> alpha(2);
	scalar_variable_impl<double> beta(1);
	detail::matrix_stats_t orig_stats2 = detail::matrix_stats;
	eigen::block_multi_vector::ptr res = mv1->gemm(*mv, B, alpha, beta);
	detail::matrix_stats.print_diff(orig_stats2);
	gettimeofday(&end, NULL);
	printf("simulated gemm takes %.3f seconds\n", time_diff(start, end));
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

	multiply_blocks = 8;
	group = safs::safs_file_group::create(safs::get_sys_RAID_conf(),
			safs::safs_file_group::group_t::NAIVE);
	printf("multiply block size: %ld\n", multiply_blocks);
	if (group)
		printf("file group type: %s\n", group->get_name().c_str());
	else
		printf("file group type: rand permute\n");

	test_gemm(true, 4, 1, 128, 0);
	test_gemm(true, 64, 1, 8, 0);
	test_MvTransMv(true, 4, 1, 128, 0);
	test_MvTransMv(true, 64, 1, 8, 0);

	test_gemm(false, 4, 8, 128, 4);
	test_gemm(false, 64, 1, 8, 0);
	test_MvTransMv(false, 4, 8, 128, 4);
	test_MvTransMv(false, 64, 1, 8, 0);

#if 0
	test_gemm_simple(false, 128 * 4, 4);
	test_MvTransMv_simple(false, 128 * 4, 4);
#endif
	destroy_flash_matrix();
}
