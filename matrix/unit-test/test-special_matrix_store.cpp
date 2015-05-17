#include "mapply_matrix_store.h"
#include "local_matrix_store.h"

using namespace fm;

class mapply_add_op: public detail::portion_mapply_op
{
	const bulk_operate &op;
public:
	mapply_add_op(size_t out_num_rows, size_t out_num_cols,
			const scalar_type &type): detail::portion_mapply_op(
				out_num_rows, out_num_cols,
				type), op(type.get_basic_ops().get_add()) {
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		detail::mapply2(*ins[0], *ins[1], op, out);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr(new mapply_add_op(
					get_out_num_cols(), get_out_num_rows(), get_output_type()));
	}
};

class multiply_op: public detail::portion_mapply_op
{
	detail::local_matrix_store::const_ptr small_mat;
public:
	multiply_op(detail::local_matrix_store::const_ptr small_mat, size_t out_num_rows,
			const scalar_type &type): detail::portion_mapply_op(
				out_num_rows, small_mat->get_num_cols(), type) {
		this->small_mat = small_mat;
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		out.reset_data();
		const bulk_operate &left_op = get_output_type().get_basic_ops().get_multiply();
		const bulk_operate &right_op = get_output_type().get_basic_ops().get_add();
		detail::inner_prod(*ins[0], *small_mat, left_op, right_op, out);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const;
};

class t_multiply_op: public detail::portion_mapply_op
{
	multiply_op op;
public:
	t_multiply_op(const multiply_op &_op): detail::portion_mapply_op(
			_op.get_out_num_cols(), _op.get_out_num_rows(),
			_op.get_output_type()), op(_op) {
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		std::vector<detail::local_matrix_store::const_ptr> t_ins(1);
		t_ins[0] = detail::local_matrix_store::cast(ins[0]->transpose());
		detail::local_matrix_store::ptr t_out
			= detail::local_matrix_store::cast(out.transpose());
		op.run(t_ins, *t_out);
	}

	virtual detail::portion_mapply_op::const_ptr transpose() const {
		return detail::portion_mapply_op::const_ptr();
	}
};

detail::portion_mapply_op::const_ptr multiply_op::transpose() const
{
	return detail::portion_mapply_op::const_ptr(new t_multiply_op(*this));
}

void verify_result(mem_dense_matrix::ptr res1, mem_dense_matrix::ptr res2)
{
	assert(res1->get_num_rows() == res2->get_num_rows());
	assert(res1->get_num_cols() == res2->get_num_cols());
	res1->materialize_self();
	res2->materialize_self();
	for (size_t i = 0; i < res1->get_num_rows(); i++)
		for (size_t j = 0; j < res1->get_num_cols(); j++)
			assert(res1->get<double>(i, j) == res2->get<double>(i, j));
}

void test_mapply_matrix_store(size_t num_rows, size_t num_cols,
		matrix_layout_t layout)
{
	printf("mapply matrix store for %s %s-major matrix\n",
			num_rows > num_cols ? "tall" : "wide",
			layout == matrix_layout_t::L_ROW ? "row" : "col");
	mem_dense_matrix::ptr mat1 = mem_dense_matrix::create_rand<double>(
			-1, 1, num_rows, num_cols, layout, -1);
	mem_dense_matrix::ptr mat2 = mem_dense_matrix::create_rand<double>(
			-1, 1, num_rows, num_cols, layout, -1);

	// Create mapply matrix store.
	std::vector<detail::mem_matrix_store::const_ptr> in_stores(2);
	in_stores[0] = detail::mem_matrix_store::cast(mat1->get_raw_store());
	in_stores[1] = detail::mem_matrix_store::cast(mat2->get_raw_store());
	detail::portion_mapply_op::const_ptr op(new mapply_add_op(
				mat1->get_num_rows(), mat1->get_num_cols(), mat1->get_type()));
	detail::mapply_matrix_store::const_ptr mapply_store(new detail::mapply_matrix_store(
				in_stores, op, mat1->store_layout(),
				mat1->get_num_rows(), mat1->get_num_cols()));
	mem_dense_matrix::ptr res2 = mem_dense_matrix::create(
			detail::mem_matrix_store::cast(mapply_store->materialize()));
	for (size_t i = 0; i < res2->get_num_rows(); i++)
		for (size_t j = 0; j < res2->get_num_cols(); j++)
			assert(res2->get<double>(i, j) == mat1->get<double>(i,
						j) + mat2->get<double>(i, j));

	mem_dense_matrix::ptr t_res
		= mem_dense_matrix::create(detail::mem_matrix_store::cast(
				detail::virtual_matrix_store::cast(mapply_store->transpose())->materialize()));
	for (size_t i = 0; i < t_res->get_num_rows(); i++)
		for (size_t j = 0; j < t_res->get_num_cols(); j++)
			assert(t_res->get<double>(i, j) == mat1->get<double>(j,
						i) + mat2->get<double>(j, i));

	// Input of matrix multiplication.
	mem_dense_matrix::ptr small_mat = mem_dense_matrix::create_rand<double>(
			-1, 1, num_cols, num_cols, matrix_layout_t::L_COL, -1);
	assert(small_mat->get_data().get_num_portions() == 1);
	detail::local_matrix_store::const_ptr local_mat
		= dynamic_cast<const detail::mem_matrix_store &>(
				small_mat->get_data()).get_portion(0);

	// Matrix multiplication in a normal way.
	mem_dense_matrix::ptr large_mat = res2;
	mem_dense_matrix::ptr res1 = mem_dense_matrix::cast(
			large_mat->multiply(*small_mat, large_mat->store_layout()));

	large_mat = mem_dense_matrix::create(mapply_store);
	in_stores.resize(1);
	in_stores[0] = detail::mem_matrix_store::cast(large_mat->get_raw_store());
	op = detail::portion_mapply_op::const_ptr(new multiply_op(
				local_mat, large_mat->get_num_rows(), large_mat->get_type()));
	mapply_store = detail::mapply_matrix_store::const_ptr(new detail::mapply_matrix_store(
				in_stores, op, large_mat->store_layout(),
				large_mat->get_num_rows(), large_mat->get_num_cols()));
	res2 = mem_dense_matrix::create(
			detail::mem_matrix_store::cast(mapply_store->materialize()));
	verify_result(res1, res2);

	for (size_t k = 0; k < mapply_store->get_num_portions(); k++) {
		detail::local_matrix_store::const_ptr part = mapply_store->get_portion(k);
		size_t start_row = part->get_global_start_row();
		size_t start_col = part->get_global_start_col();
		for (size_t i = 0; i < part->get_num_rows(); i++)
			for (size_t j = 0; j < part->get_num_cols(); j++)
				assert(part->get<double>(i, j) == res1->get<double>(
							i + start_row, j + start_col));
	}

	t_res = mem_dense_matrix::create(detail::mem_matrix_store::cast(
				detail::virtual_matrix_store::cast(mapply_store->transpose())->materialize()));
	verify_result(mem_dense_matrix::cast(res1->transpose()), t_res);
}

int main()
{
	test_mapply_matrix_store(100000, 10, matrix_layout_t::L_COL);
	test_mapply_matrix_store(100000, 10, matrix_layout_t::L_ROW);
}
