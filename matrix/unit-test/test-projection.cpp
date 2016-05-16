#include "sparse_matrix.h"
#include "dense_matrix.h"
#include "project_matrix_store.h"

using namespace fm;
using namespace fm::detail;

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	sparse_project_matrix_store::ptr proj
		= sparse_project_matrix_store::create_sparse_rand(1000000, 100,
				matrix_layout_t::L_ROW, get_scalar_type<double>(), 0.001);
	size_t nnz = 0;
	double sum = 0;
	for (size_t i = 0; i < proj->get_num_portions(); i++) {
		matrix_store::const_ptr mat = proj;
		local_sparse_matrix_store::const_ptr portion
			= std::dynamic_pointer_cast<const local_sparse_matrix_store>(
					mat->get_portion(i));
		if (portion == NULL)
			continue;

		assert(portion->store_layout() == matrix_layout_t::L_ROW);
		std::vector<off_t> idxs;
		size_t local_nnz = 0;
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
		assert(local_nnz == portion->get_nnz());
		nnz += portion->get_nnz();
	}
	assert(nnz == proj->get_nnz());
	printf("sum: %f\n", sum);

	dense_matrix::ptr data = dense_matrix::create_randu<double>(0, 1, 100,
			proj->get_num_rows(), matrix_layout_t::L_COL);
	dense_matrix::ptr proj_mat = dense_matrix::create(proj);
	dense_matrix::ptr res = data->multiply(*proj_mat);
	res->materialize_self();

	dense_matrix::ptr dense_proj = dense_matrix::create(proj->conv_dense());
	scalar_variable::ptr proj_sum = dense_proj->sum();
	printf("dense proj sum: %f\n", scalar_variable::get_val<double>(*proj_sum));
	dense_matrix::ptr res2 = data->multiply(*dense_proj);
	dense_matrix::ptr diff = res->minus(*res2);
	scalar_variable::ptr diff_sum = diff->abs()->sum();
	assert(diff_sum->get_type() == get_scalar_type<double>());
	printf("diff: %f\n", scalar_variable::get_val<double>(*diff_sum));

	destroy_flash_matrix();
}
