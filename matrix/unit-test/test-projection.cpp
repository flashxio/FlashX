#include "dense_matrix.h"
#include "project_matrix_store.h"

using namespace fm;
using namespace fm::detail;

int main()
{
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
}
