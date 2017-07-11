/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "factor.h"
#include "data_frame.h"
#include "mapply_matrix_store.h"

namespace fm
{

namespace
{

class factor_map
{
public:
	typedef std::shared_ptr<factor_map> ptr;

	virtual void map_data(const void *in, factor_value_t *vals,
			size_t num_eles) const = 0;
	virtual const scalar_type &get_input_type() const = 0;
	virtual size_t get_num_levels() const = 0;
};

template<class T>
class factor_map_impl: public factor_map
{
	std::unordered_map<T, factor_value_t> map;
public:
	static ptr create(detail::vec_store::const_ptr vec) {
		assert(vec->get_type() == get_scalar_type<T>());
		detail::mem_vec_store::const_ptr mem_vec
			= std::dynamic_pointer_cast<const detail::mem_vec_store>(vec);
		assert(mem_vec);
		factor_map_impl<T> *ret = new factor_map_impl<T>();
		size_t id = 0;
		for (size_t i = 0; i < vec->get_length(); i++) {
			auto ins_ret = ret->map.insert(std::pair<T, factor_value_t>(
						mem_vec->get<T>(i), id));
			// If we insert the value successfully, it means the value didn't
			// exist before. We assign an Id for the value and increase Id.
			if (ins_ret.second)
				id++;
		}
		assert(id == ret->map.size());
		return ptr(ret);
	}

	void map_data(const void *in, factor_value_t *vals, size_t num_eles) const {
		const T *t_in = reinterpret_cast<const T *>(in);
		for (size_t i = 0; i < num_eles; i++) {
			auto it = map.find(t_in[i]);
			assert(it != map.end());
			vals[i] = it->second;
		}
	}
	const scalar_type &get_input_type() const {
		return get_scalar_type<T>();
	}
	virtual size_t get_num_levels() const {
		return map.size();
	}
};

class remap_portion_op: public detail::portion_mapply_op
{
	factor_map::ptr map;
public:
	remap_portion_op(factor_map::ptr map, size_t num_rows,
			size_t num_cols): detail::portion_mapply_op(num_rows, num_cols,
				get_scalar_type<factor_value_t>()) {
		this->map = map;
	}

	virtual portion_mapply_op::const_ptr transpose() const {
		return portion_mapply_op::const_ptr(new remap_portion_op(map,
					get_out_num_cols(), get_out_num_rows()));
	}

	virtual void run(const std::vector<detail::local_matrix_store::const_ptr> &ins,
			detail::local_matrix_store &out) const {
		assert(ins.size() == 1);
		assert(ins[0]->get_type() == map->get_input_type());
		assert(out.get_type() == get_scalar_type<factor_value_t>());
		const char *addr = ins[0]->get_raw_arr();
		assert(addr);
		map->map_data(addr, reinterpret_cast<factor_value_t *>(out.get_raw_arr()),
				ins[0]->get_num_rows() * ins[0]->get_num_cols());
	}

	virtual std::string to_string(
			const std::vector<detail::matrix_store::const_ptr> &mats) const {
		return std::string("remap(") + mats[0]->get_name() + ")";
	}
};

}

factor_col_vector::ptr factor_col_vector::create(dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with more than one col to a vector";
		return ptr();
	}
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);

	// Let's cast the element type to some specific types that
	// we create a map on. Otherwise, the key type of the map will be
	// arbitrary.
	if (mat->get_type().is_floating_point())
		mat = mat->cast_ele_type(get_scalar_type<double>());
	else
		mat = mat->cast_ele_type(get_scalar_type<long>());

	// Create the map from the input matrix.
	col_vec::ptr cvec = col_vec::create(mat);
	auto count = cvec->get_type().get_agg_ops().get_count();
	// We collect all of the unique values in the vector,
	// but we don't need these values to be sorted in a particular order.
	data_frame::ptr df = cvec->groupby(count, true, false);
	factor_map::ptr map;
	detail::vec_store::ptr val_vec = df->get_vec(0);
	// We should sort the values so that smaller values will be assigned
	// smaller Id.
	val_vec->sort();
	if (df->get_vec(0)->get_type().is_floating_point())
		map = factor_map_impl<double>::create(val_vec);
	else
		map = factor_map_impl<long>::create(val_vec);

	// Construct the matrix with data mapped automatically.
	std::vector<detail::matrix_store::const_ptr> in_mats(1,
			mat->get_raw_store());
	detail::portion_mapply_op::const_ptr op(new remap_portion_op(map,
				mat->get_num_rows(), mat->get_num_cols()));
	detail::matrix_store::const_ptr store = detail::matrix_store::const_ptr(
			new detail::mapply_matrix_store(in_mats, op, matrix_layout_t::L_COL));
	auto ret = new factor_col_vector(factor(map->get_num_levels()), store);
	ret->uniq_vals = df->get_vec(0);
	ret->cnts = df->get_vec(1);
	return ptr(ret);
}

factor_col_vector::ptr factor_col_vector::create(const factor &f,
		dense_matrix::ptr mat)
{
	if (mat->get_num_cols() > 1) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't convert a matrix with more than one col to a vector";
		return ptr();
	}
	if (mat->get_type().is_floating_point()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The type of the input matrix can't be floating point";
		return ptr();
	}
	mat = mat->cast_ele_type(get_scalar_type<factor_value_t>());
	if (mat->get_data().store_layout() == matrix_layout_t::L_ROW)
		mat = mat->conv2(matrix_layout_t::L_COL);
	return ptr(new factor_col_vector(f, mat->get_raw_store()));
}

static detail::matrix_store::ptr create_factor_store(size_t len,
		int num_nodes, bool in_mem, const set_operate &op)
{
	detail::matrix_store::ptr ret = detail::matrix_store::create(len, 1,
			matrix_layout_t::L_COL, get_scalar_type<factor_value_t>(),
			num_nodes, in_mem);
	if (ret == NULL)
		return detail::matrix_store::ptr();

	ret->set_data(op);
	return ret;
}

factor_col_vector::factor_col_vector(const factor &_f, size_t len,
		int num_nodes, bool in_mem, const set_operate &op): col_vec(
			create_factor_store(len, num_nodes, in_mem, op)), f(_f) {
}

}
