/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
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

#include <unordered_map>

#include "matrix_ops.h"

using namespace fm;

namespace fmr
{

/*
 * R has only two data types in matrix multiplication: integer and numeric.
 * So we only need to predefine a small number of basic operations with
 * different types.
 */

static basic_ops_impl<int, int, int> R_basic_ops_II;
static basic_ops_impl<double, double, double> R_basic_ops_DD;

static basic_uops_impl<int, int> R_basic_uops_I;
static basic_uops_impl<double, double> R_basic_uops_D;

class generic_bulk_operate
{
	std::string name;
	// a bulk_operate for a different type.
	std::unordered_map<int, bulk_operate::const_ptr> ops;
public:
	generic_bulk_operate(const std::string &name) {
		this->name = name;
	}

	bulk_operate::const_ptr get_op(prim_type type) const {
		int key = type;
		auto it = ops.find(key);
		if (it == ops.end())
			return bulk_operate::const_ptr();
		else
			return it->second;
	}

	void add_op(bulk_operate::const_ptr op) {
		int key = op->get_left_type().get_type();
		ops.insert(std::pair<int, bulk_operate::const_ptr>(key, op));
	}

	const std::string &get_name() const {
		return name;
	}
};

class generic_bulk_uoperate
{
	std::string name;
	// a bulk_uoperate for a different type.
	std::vector<bulk_uoperate::const_ptr> ops;
public:
	generic_bulk_uoperate(const std::string &name) {
		this->name = name;
	}

	bulk_uoperate::const_ptr get_op(prim_type type) const {
		int off = type;
		if (type >= ops.size())
			return bulk_uoperate::const_ptr();
		else
			return ops[off];
	}

	void add_op(bulk_uoperate::const_ptr op) {
		int off = op->get_input_type().get_type();
		ops.resize(off + 1);
		ops[off] = op;
	}

	const std::string &get_name() const {
		return name;
	}
};

static std::vector<generic_bulk_operate> bulk_ops;
static std::vector<generic_bulk_uoperate> bulk_uops;

/*
 * Register a binary UDF.
 * A user has to provide UDFs for all different types.
 */
void register_udf(const std::vector<bulk_operate::const_ptr> &ops,
		const std::string &name)
{
	bulk_ops.emplace_back(name);
	for (size_t i = 0; i < ops.size(); i++)
		bulk_ops.back().add_op(ops[i]);
}

/*
 * Register a unary UDF.
 * A user has to provide UDFs for all different types.
 */
void register_udf(const std::vector<bulk_uoperate::const_ptr> &ops,
		const std::string &name)
{
	bulk_uops.emplace_back(name);
	for (size_t i = 0; i < ops.size(); i++)
		bulk_uops.back().add_op(ops[i]);
}

static bulk_operate::const_ptr _get_op(basic_ops::op_idx bo_idx, int noperands,
		prim_type type)
{
	if (noperands != 2) {
		fprintf(stderr, "This isn't a binary operator\n");
		return bulk_operate::const_ptr();
	}
	if (bo_idx < 0) {
		fprintf(stderr, "invalid operator index\n");
		return bulk_operate::const_ptr();
	}

	bulk_operate::const_ptr op;
	if (bo_idx < basic_ops::op_idx::NUM_OPS) {
		basic_ops *ops = NULL;
		if (type == prim_type::P_DOUBLE)
			ops = &R_basic_ops_DD;
		else if (type == prim_type::P_INTEGER)
			ops = &R_basic_ops_II;
		else {
			fprintf(stderr, "wrong type\n");
			return bulk_operate::const_ptr();
		}
		op = bulk_operate::conv2ptr(*ops->get_op(bo_idx));
		if (op == NULL) {
			fprintf(stderr, "invalid basic binary operator\n");
			return bulk_operate::const_ptr();
		}
	}
	else if ((size_t) (bo_idx - basic_ops::op_idx::NUM_OPS) < bulk_ops.size()) {
		size_t off = bo_idx - basic_ops::op_idx::NUM_OPS;
		op = bulk_ops[off].get_op(type);
		if (op == NULL) {
			fprintf(stderr,
					"can't find the specified operation with the right type\n");
			return bulk_operate::const_ptr();
		}
	}
	else {
		fprintf(stderr, "can't find the specified operation\n");
		return bulk_operate::const_ptr();
	}
	return op;
}

/*
 * Get a binary operator.
 */
bulk_operate::const_ptr get_op(SEXP pfun, prim_type type)
{
	Rcpp::S4 fun_obj(pfun);
	Rcpp::IntegerVector info = fun_obj.slot("info");
	return _get_op((basic_ops::op_idx) info[0], info[1], type);
}

/*
 * This construct an aggregation operator from binary operators.
 */
agg_operate::const_ptr get_agg_op(SEXP pfun, const scalar_type &mat_type)
{
	Rcpp::S4 sym_op(pfun);
	Rcpp::IntegerVector agg_info = sym_op.slot("agg");
	bulk_operate::const_ptr agg_op = _get_op((basic_ops::op_idx) agg_info[0],
			agg_info[1], mat_type.get_type());
	Rcpp::IntegerVector combine_info = sym_op.slot("combine");
	bulk_operate::const_ptr combine_op;
	if (combine_info[0] >= 0)
		combine_op = _get_op((basic_ops::op_idx) combine_info[0],
				combine_info[1], agg_op->get_output_type().get_type());
	return agg_operate::create(agg_op, combine_op);
}

/*
 * Get a unary operator.
 */
bulk_uoperate::const_ptr get_uop(SEXP pfun, prim_type type)
{
	Rcpp::S4 fun_obj(pfun);
	Rcpp::IntegerVector info = fun_obj.slot("info");
	basic_uops::op_idx bo_idx = (basic_uops::op_idx) info[0];
	int noperands = info[1];
	if (noperands != 1) {
		fprintf(stderr, "This isn't a unary operator\n");
		return NULL;
	}

	bulk_uoperate::const_ptr op;
	if (bo_idx < basic_uops::op_idx::NUM_OPS) {
		basic_uops *ops = NULL;
		if (type == prim_type::P_DOUBLE)
			ops = &R_basic_uops_D;
		else if (type == prim_type::P_INTEGER)
			ops = &R_basic_uops_I;
		else {
			fprintf(stderr, "wrong type\n");
			return NULL;
		}

		op = bulk_uoperate::conv2ptr(*ops->get_op(bo_idx));
		if (op == NULL) {
			fprintf(stderr, "invalid basic unary operator\n");
			return NULL;
		}
	}
	else if ((size_t) (bo_idx - basic_uops::op_idx::NUM_OPS) < bulk_uops.size()) {
		size_t off = bo_idx - basic_uops::op_idx::NUM_OPS;
		op = bulk_uops[off].get_op(type);
		if (op == NULL) {
			fprintf(stderr,
					"can't find the specified unary operation with right type\n");
			return bulk_uoperate::const_ptr();
		}
	}
	else {
		fprintf(stderr, "can't find the specified unary operation\n");
		return bulk_uoperate::const_ptr();
	}
	return op;
}

static int _get_op_id(const std::string &name)
{
	for (size_t i = 0; i < bulk_ops.size(); i++)
		if (bulk_ops[i].get_name() == name)
			return i + basic_ops::op_idx::NUM_OPS;
	return -1;
}

static int _get_uop_id(const std::string &name)
{
	for (size_t i = 0; i < bulk_uops.size(); i++)
		if (bulk_uops[i].get_name() == name)
			return i + basic_uops::op_idx::NUM_OPS;
	return -1;
}

op_id_t get_op_id(const std::string &name)
{
	// TODO I should use a hashtable.
	if (name == "add")
		return basic_ops::op_idx::ADD;
	else if (name == "+")
		return basic_ops::op_idx::ADD;
	else if (name == "sub")
		return basic_ops::op_idx::SUB;
	else if (name == "-")
		return basic_ops::op_idx::SUB;
	else if (name == "mul")
		return basic_ops::op_idx::MUL;
	else if (name == "*")
		return basic_ops::op_idx::MUL;
	else if (name == "div")
		return basic_ops::op_idx::DIV;
	else if (name == "/")
		return basic_ops::op_idx::DIV;
	else if (name == "min")
		return basic_ops::op_idx::MIN;
	else if (name == "max")
		return basic_ops::op_idx::MAX;
	else if (name == "pow")
		return basic_ops::op_idx::POW;
	else if (name == "eq")
		return basic_ops::op_idx::EQ;
	else if (name == "==")
		return basic_ops::op_idx::EQ;
	else if (name == "neq")
		return basic_ops::op_idx::NEQ;
	else if (name == "!=")
		return basic_ops::op_idx::NEQ;
	else if (name == "gt")
		return basic_ops::op_idx::GT;
	else if (name == ">")
		return basic_ops::op_idx::GT;
	else if (name == "ge")
		return basic_ops::op_idx::GE;
	else if (name == ">=")
		return basic_ops::op_idx::GE;
	else if (name == "lt")
		return basic_ops::op_idx::LT;
	else if (name == "<")
		return basic_ops::op_idx::LT;
	else if (name == "le")
		return basic_ops::op_idx::LE;
	else if (name == "<=")
		return basic_ops::op_idx::LE;
	else if (name == "|")
		return basic_ops::op_idx::OR;
	else if (name == "&")
		return basic_ops::op_idx::AND;
	else
		return _get_op_id(name);
}

op_id_t get_uop_id(const std::string &name)
{
	if (name == "neg")
		return basic_uops::op_idx::NEG;
	else if (name == "sqrt")
		return basic_uops::op_idx::SQRT;
	else if (name == "abs")
		return basic_uops::op_idx::ABS;
	else if (name == "not")
		return basic_uops::op_idx::NOT;
	else if (name == "ceil")
		return basic_uops::op_idx::CEIL;
	else if (name == "floor")
		return basic_uops::op_idx::FLOOR;
	else if (name == "round")
		return basic_uops::op_idx::ROUND;
	else if (name == "log")
		return basic_uops::op_idx::LOG;
	else if (name == "log2")
		return basic_uops::op_idx::LOG2;
	else if (name == "log10")
		return basic_uops::op_idx::LOG10;
	else
		return _get_uop_id(name);
}

template<class T>
class r_count_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}

	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		int *t_out = (int *) output;
		if (orig == NULL)
			t_out[0] = num_eles;
		else
			t_out[0] = (*(const int *) orig) + num_eles;
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "count";
	}
};

template<class T>
class r_which_max_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}

	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		int *t_out = (int *) output;
		const T *t_in = (const T *) in;
		if (num_eles == 0)
			return;
		T max = t_in[0];
		int idx = 0;
		for (size_t i = 1; i < num_eles; i++) {
			if (max < t_in[i]) {
				max = t_in[i];
				idx = i;
			}
		}
		t_out[0] = idx + 1;
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "which_max";
	}
};

template<class T>
class r_which_min_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		throw unsupported_exception();
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		throw unsupported_exception();
	}

	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		int *t_out = (int *) output;
		const T *t_in = (const T *) in;
		if (num_eles == 0)
			return;
		T min = t_in[0];
		int idx = 0;
		for (size_t i = 1; i < num_eles; i++) {
			if (min > t_in[i]) {
				min = t_in[i];
				idx = i;
			}
		}
		t_out[0] = idx + 1;
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual std::string get_name() const {
		return "which_min";
	}
};

template<class T>
class r_euclidean_operate: public bulk_operate
{
public:
	virtual void runAA(size_t num_eles, const void *left_arr,
			const void *right_arr, void *output_arr) const {
		const T *arr1 = (const T *) left_arr;
		const T *arr2 = (const T *) right_arr;
		T *out = (T *) output_arr;
		for (size_t i = 0; i < num_eles; i++)
			out[i] = (arr1[i] - arr2[i]) * (arr1[i] - arr2[i]);
	}
	virtual void runAE(size_t num_eles, const void *left_arr,
			const void *right, void *output_arr) const {
		const T *arr1 = (const T *) left_arr;
		const T v = *(const T *) right;
		T *out = (T *) output_arr;
		for (size_t i = 0; i < num_eles; i++)
			out[i] = (arr1[i] - v) * (arr1[i] - v);
	}
	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr, void *output_arr) const {
		const T v = *(const T *) left;
		const T *arr2 = (const T *) right_arr;
		T *out = (T *) output_arr;
		for (size_t i = 0; i < num_eles; i++)
			out[i] = (v - arr2[i]) * (v - arr2[i]);
	}

	virtual void runAgg(size_t num_eles, const void *in, const void *orig,
			void *output) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<T>();
	}
	virtual std::string get_name() const {
		return "euclidean";
	}
};

void init_udf_ext()
{
	std::vector<bulk_operate::const_ptr> ops;
	// Add count.
	ops.push_back(bulk_operate::const_ptr(new r_count_operate<bool>()));
	ops.push_back(bulk_operate::const_ptr(new r_count_operate<int>()));
	ops.push_back(bulk_operate::const_ptr(new r_count_operate<double>()));
	register_udf(ops, "count");

	ops.clear();
	ops.push_back(bulk_operate::const_ptr(new r_which_max_operate<int>()));
	ops.push_back(bulk_operate::const_ptr(new r_which_max_operate<double>()));
	register_udf(ops, "which.max");

	ops.clear();
	ops.push_back(bulk_operate::const_ptr(new r_which_min_operate<int>()));
	ops.push_back(bulk_operate::const_ptr(new r_which_min_operate<double>()));
	register_udf(ops, "which.min");

	ops.clear();
	ops.push_back(bulk_operate::const_ptr(new r_euclidean_operate<int>()));
	ops.push_back(bulk_operate::const_ptr(new r_euclidean_operate<double>()));
	register_udf(ops, "euclidean");

	std::vector<bulk_uoperate::const_ptr> uops;
	uops.push_back(bulk_uoperate::conv2ptr(
				get_scalar_type<bool>().get_type_cast(get_scalar_type<int>())));
	uops.push_back(bulk_uoperate::conv2ptr(
				get_scalar_type<double>().get_type_cast(get_scalar_type<int>())));
	register_udf(uops, "as.int");

	uops.clear();
	uops.push_back(bulk_uoperate::conv2ptr(
				get_scalar_type<bool>().get_type_cast(get_scalar_type<double>())));
	uops.push_back(bulk_uoperate::conv2ptr(
				get_scalar_type<int>().get_type_cast(get_scalar_type<double>())));
	register_udf(uops, "as.numeric");
}

}
