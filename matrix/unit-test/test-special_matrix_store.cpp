#include "mapply_matrix_store.h"
#include "local_matrix_store.h"
#include "dense_matrix.h"

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

void verify_result(dense_matrix::ptr res1, dense_matrix::ptr res2)
{
	assert(res1->get_num_rows() == res2->get_num_rows());
	assert(res1->get_num_cols() == res2->get_num_cols());
	res1->materialize_self();
	res2->materialize_self();
	const detail::mem_matrix_store &mem_res1
		= dynamic_cast<const detail::mem_matrix_store &>(res1->get_data());
	const detail::mem_matrix_store &mem_res2
		= dynamic_cast<const detail::mem_matrix_store &>(res2->get_data());
	for (size_t i = 0; i < res1->get_num_rows(); i++)
		for (size_t j = 0; j < res1->get_num_cols(); j++) {
			assert(mem_res1.get<long>(i, j) == mem_res2.get<long>(i, j));
		}
}

void test_mapply_matrix_store(size_t num_rows, size_t num_cols,
		matrix_layout_t layout)
{
	printf("mapply matrix store for %s %s-major matrix\n",
			num_rows > num_cols ? "tall" : "wide",
			layout == matrix_layout_t::L_ROW ? "row" : "col");
	dense_matrix::ptr mat1 = dense_matrix::create_rand<long>(
			-1, 1, num_rows, num_cols, layout, -1);
	dense_matrix::ptr mat2 = dense_matrix::create_rand<long>(
			-1, 1, num_rows, num_cols, layout, -1);

	// Create mapply matrix store.
	std::vector<detail::matrix_store::const_ptr> in_stores(2);
	in_stores[0] = mat1->get_raw_store();
	in_stores[1] = mat2->get_raw_store();
	detail::portion_mapply_op::const_ptr op(new mapply_add_op(
				mat1->get_num_rows(), mat1->get_num_cols(), mat1->get_type()));
	detail::mapply_matrix_store::const_ptr mapply_store(new detail::mapply_matrix_store(
				in_stores, op, mat1->store_layout(),
				mat1->get_num_rows(), mat1->get_num_cols()));
	dense_matrix::ptr res2 = dense_matrix::create(mapply_store->materialize());
	{
		const detail::mem_matrix_store &mem_mat1
			= dynamic_cast<const detail::mem_matrix_store &>(mat1->get_data());
		const detail::mem_matrix_store &mem_mat2
			= dynamic_cast<const detail::mem_matrix_store &>(mat2->get_data());
		const detail::mem_matrix_store &mem_res
			= dynamic_cast<const detail::mem_matrix_store &>(res2->get_data());
		for (size_t i = 0; i < res2->get_num_rows(); i++)
			for (size_t j = 0; j < res2->get_num_cols(); j++)
				assert(mem_res.get<long>(i, j) == mem_mat1.get<long>(i,
							j) + mem_mat2.get<long>(i, j));
	}

	dense_matrix::ptr t_res = dense_matrix::create(
			detail::virtual_matrix_store::cast(mapply_store->transpose())->materialize());
	{
		const detail::mem_matrix_store &mem_mat1
			= dynamic_cast<const detail::mem_matrix_store &>(mat1->get_data());
		const detail::mem_matrix_store &mem_mat2
			= dynamic_cast<const detail::mem_matrix_store &>(mat2->get_data());
		const detail::mem_matrix_store &mem_res
			= dynamic_cast<const detail::mem_matrix_store &>(t_res->get_data());
		for (size_t i = 0; i < t_res->get_num_rows(); i++)
			for (size_t j = 0; j < t_res->get_num_cols(); j++)
				assert(mem_res.get<long>(i, j) == mem_mat1.get<long>(j,
							i) + mem_mat2.get<long>(j, i));
	}

	// Input of matrix multiplication.
	dense_matrix::ptr small_mat = dense_matrix::create_rand<long>(
			-1, 1, num_cols, num_cols, matrix_layout_t::L_COL, -1);
	assert(small_mat->get_data().get_num_portions() == 1);
	detail::local_matrix_store::const_ptr local_mat
		= dynamic_cast<const detail::mem_matrix_store &>(
				small_mat->get_data()).get_portion(0);

	// Matrix multiplication in a normal way.
	dense_matrix::ptr large_mat = res2;
	dense_matrix::ptr res1 = large_mat->multiply(*small_mat,
			large_mat->store_layout());

	large_mat = dense_matrix::create(mapply_store);
	in_stores.resize(1);
	in_stores[0] = large_mat->get_raw_store();
	op = detail::portion_mapply_op::const_ptr(new multiply_op(
				local_mat, large_mat->get_num_rows(), large_mat->get_type()));
	mapply_store = detail::mapply_matrix_store::const_ptr(new detail::mapply_matrix_store(
				in_stores, op, large_mat->store_layout(),
				large_mat->get_num_rows(), large_mat->get_num_cols()));
	res2 = dense_matrix::create(mapply_store->materialize());
	verify_result(res1, res2);

	const detail::mem_matrix_store &mem_res
		= dynamic_cast<const detail::mem_matrix_store &>(res1->get_data());
	for (size_t k = 0; k < mapply_store->get_num_portions(); k++) {
		detail::local_matrix_store::const_ptr part = mapply_store->get_portion(k);
		size_t start_row = part->get_global_start_row();
		size_t start_col = part->get_global_start_col();
		for (size_t i = 0; i < part->get_num_rows(); i++)
			for (size_t j = 0; j < part->get_num_cols(); j++)
				assert(part->get<long>(i, j) == mem_res.get<long>(
							i + start_row, j + start_col));
	}

	t_res = dense_matrix::create(detail::virtual_matrix_store::cast(
				mapply_store->transpose())->materialize());
	verify_result(res1->transpose(), t_res);
}

int main()
{
	test_mapply_matrix_store(100000, 10, matrix_layout_t::L_COL);
	test_mapply_matrix_store(100000, 10, matrix_layout_t::L_ROW);
}
