#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "project_matrix_store.h"

using namespace fm;
using namespace fm::detail;

void test_whole()
{
	sparse_project_matrix_store::ptr proj_store
		= sparse_project_matrix_store::create_sparse_rand(1000000, 100,
				matrix_layout_t::L_ROW, get_scalar_type<double>(), 0.001);
	size_t nnz = 0;
	double sum = 0;
	for (size_t i = 0; i < proj_store->get_num_portions(); i++) {
		matrix_store::const_ptr mat = proj_store;
		lsparse_row_matrix_store::const_ptr portion
			= std::dynamic_pointer_cast<const lsparse_row_matrix_store>(
					mat->get_portion(i));
		if (portion == NULL)
			continue;

		assert(portion->store_layout() == matrix_layout_t::L_ROW);
		std::vector<off_t> idxs;
		size_t local_nnz = 0;
		size_t orig_nrow = portion->get_num_rows();
		size_t num_portions = portion->get_num_rows() / 128;
		for (size_t pidx = 0; pidx < num_portions; pidx++) {
			size_t local_nrow = std::min(128UL, orig_nrow - pidx * 128);
			lsparse_row_matrix_store *mutable_portion
				= const_cast<lsparse_row_matrix_store *>(portion.get());
			mutable_portion->resize(pidx * 128, 0, local_nrow,
					portion->get_num_cols());
			size_t orig_local_nnz = local_nnz;
			for (size_t j = 0; j < portion->get_num_rows(); j++) {
				const char *row = portion->get_row_nnz(j, idxs);
				if (row == NULL) {
					assert(idxs.size() == 0);
					continue;
				}

				const double *drow = (const double *) row;
				local_nnz += idxs.size();
				for (size_t k = 0; k < idxs.size(); k++) {
					assert(idxs[k] >= 0 && idxs[k] < portion->get_num_rows());
					sum += drow[k];
				}
			}
		}
		assert(local_nnz == portion->get_nnz());
		nnz += portion->get_nnz();
	}
	assert(nnz == proj_store->get_nnz());
	printf("sum: %f\n", sum);

	dense_matrix::ptr mean = dense_matrix::create_randu<double>(
			0, 1, proj_store->get_num_rows(), 1, matrix_layout_t::L_COL);
	std::vector<dense_matrix::ptr> tmps(2);
	tmps[0] = mean;
	tmps[1] = dense_matrix::create(proj_store);
	dense_matrix::ptr proj = dense_matrix::cbind(tmps);

	dense_matrix::ptr data = dense_matrix::create_randu<double>(0, 1, 100,
			proj->get_num_rows(), matrix_layout_t::L_COL);
	dense_matrix::ptr res = data->multiply(*proj);
	res->materialize_self();

	dense_matrix::ptr dense_proj = dense_matrix::create(proj_store);
	dense_proj->materialize_self();
	tmps[0] = mean;
	tmps[1] = dense_proj;
	dense_proj = dense_matrix::cbind(tmps);

	scalar_variable::ptr proj_sum = dense_proj->sum();
	printf("dense proj sum: %f\n", scalar_variable::get_val<double>(*proj_sum));
	dense_matrix::ptr res2 = data->multiply(*dense_proj);
	dense_matrix::ptr diff = res->minus(*res2);
	scalar_variable::ptr diff_sum = diff->abs()->sum();
	assert(diff_sum->get_type() == get_scalar_type<double>());
	printf("diff: %f\n", scalar_variable::get_val<double>(*diff_sum));
}

void test_ltranspose()
{
	matrix_store::const_ptr proj_store
		= sparse_project_matrix_store::create_sparse_rand(1000000, 100,
				matrix_layout_t::L_ROW, get_scalar_type<double>(), 0.001);
	auto part = proj_store->get_portion(1);
	auto t = std::dynamic_pointer_cast<const local_matrix_store>(part->transpose());

	size_t num_bytes = t->get_num_rows() * t->get_num_cols() * t->get_entry_size();
	const char *arr = part->get_raw_arr();
	const char *tarr = t->get_raw_arr();
	assert(arr);
	assert(tarr);
	assert(memcmp(arr, tarr, num_bytes) == 0);
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

	test_whole();
	test_ltranspose();

	destroy_flash_matrix();
}
