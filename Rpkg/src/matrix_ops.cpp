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
#include "mem_worker_thread.h"
#include "local_vec_store.h"
#include "bulk_operate_impl.h"
#include "dense_matrix.h"

using namespace fm;

R_type trans_FM2R(const scalar_type &type);

namespace fmr
{

/*
 * We need to redefine the basic operations for R.
 */

template<class T, bool is_logical>
bool R_is_na(T val)
{
	fprintf(stderr, "unknown type for NA\n");
	return false;
}

template<>
bool R_is_na<int, true>(int val)
{
	return val == NA_LOGICAL;
}

template<>
bool R_is_na<int, false>(int val)
{
	return val == NA_INTEGER;
}

template<>
bool R_is_na<double, true>(double val)
{
	// We can't check `val' against `NA_REAL'. Instead, we should use `ISNA'.
	//
	// The reason is explained below:
	// Basically, a double represents a rounded real number in the following
	// notation (see also the wikipedia article):
	//
	// \textrm{sign}\times 2^{e} \times 1.F .
	//
	// The sign is represented by 1 bit, the exponent e by 11 bits and
	// the mantissa F by 52 bits, so we have 64 bits in total. The special value
	// NaN (and also \pmInf) is coded using values of e that are not used to
	// represent numbers. NaN is represented by e=0x7ff (hexadecimal) and F\not=0.
	// The important thing is that it does not matter what the value of F is
	// when representing NaN. This leaves developers with lots of room in
	// the mantissa to give different meanings to NaN. In R the developers
	// chose F=1954 in the mantissa to represent NA. A C-level function called
	// R_IsNA detects the 1954 in NaN values.
	return ISNA(val);
}

template<>
bool R_is_na<double, false>(double val)
{
	return ISNA(val);
}

template<class T, bool is_logical>
T R_get_na()
{
	fprintf(stderr, "unknown type of NA\n");
	return 0;
}

template<>
int R_get_na<int, true>()
{
	return NA_LOGICAL;
}

template<>
int R_get_na<int, false>()
{
	return NA_INTEGER;
}

template<>
double R_get_na<double, true>()
{
	return NA_REAL;
}

template<>
double R_get_na<double, false>()
{
	return NA_REAL;
}

template<class T, bool is_logical>
R_type get_Rtype()
{
	fprintf(stderr, "unknown type\n");
	return R_type::R_NTYPES;
}

template<>
R_type get_Rtype<int, true>()
{
	return R_type::R_LOGICAL;
}

template<>
R_type get_Rtype<int, false>()
{
	return R_type::R_INT;
}

template<>
R_type get_Rtype<double, true>()
{
	return R_type::R_REAL;
}

template<>
R_type get_Rtype<double, false>()
{
	return R_type::R_REAL;
}

//////////////////////// binary operators /////////////////////////

template<class Type, bool is_logical>
struct min
{
	static std::string get_name() {
		return "min";
	}
	static Type get_agg_init() {
		return std::numeric_limits<Type>::max();
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return std::min(e1, e2);
	}
};

template<>
struct min<double, false>
{
	static std::string get_name() {
		return "min";
	}
	static double get_agg_init() {
		return std::numeric_limits<double>::max();
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::min(e1, e2);
	}
};

template<class Type, bool is_logical>
struct max
{
	static std::string get_name() {
		return "max";
	}
	static Type get_agg_init() {
		// min() for integers is used for NA in R.
		return std::numeric_limits<Type>::min() + 1;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return std::max(e1, e2);
	}
};

template<>
struct max<double, false>
{
	static std::string get_name() {
		return "max";
	}
	static double get_agg_init() {
		// We need to define the minimum float-point differently.
		return -std::numeric_limits<double>::max();
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::max(e1, e2);
	}
};

template<class Type, bool is_logical>
struct mod
{
	static std::string get_name() {
		return "%";
	}
	static Type get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return e1 % e2;
	}
};

template<>
struct mod<double, false>
{
	static std::string get_name() {
		return "%";
	}
	static double get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const double &e1, const double &e2) const {
		return std::fmod(e1, e2);
	}
};

template<class Type, bool is_logical>
struct add {
	static std::string get_name() {
		return "+";
	}
	static Type get_agg_init() {
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return e1 + e2;
	}
};

template<class Type, bool is_logical>
struct sub {
	static std::string get_name() {
		return "-";
	}
	static Type get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return e1 - e2;
	}
};

template<class Type, bool is_logical>
struct multiply {
	static std::string get_name() {
		return "*";
	}
	static Type get_agg_init() {
		return 1;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return e1 * e2;
	}
};

// Division is special. Its output should be float point.
// Therefore, we convert both input values to float point.
template<class Type, bool is_logical>
struct divide {
	static std::string get_name() {
		return "/";
	}
	static Type get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const Type &e1, const Type &e2) const {
		double d1 = e1;
		double d2 = e2;
		return d1 / d2;
	}
};

template<class Type, bool is_logical>
struct idiv {
	divide<Type, is_logical> div;
	static std::string get_name() {
		return "%/%";
	}
	static Type get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		return std::floor(div(e1, e2));
	}
};

template<class Type, bool is_logical>
struct pow {
	static std::string get_name() {
		return "pow";
	}
	static Type get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e1, const Type &e2) const {
		// If e1 is 1, whatever e2 is, the output is 1.
		// If e2 is 0, whatever e1 is, the output is 1.
		if (e1 == 1 || e2 == 0)
			return 1;
		else
			return std::pow(e1, e2);
	}
};

template<>
struct pow<double, false> {
	static std::string get_name() {
		return "pow";
	}
	static double get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const double &e1, const double &e2) const {
		if (e1 == 1 || e2 == 0)
			return 1;
		// If e1 is -Inf and e2 isn't an integer, C++ returns Inf,
		// but R wants NaN.
		else if (e1 == -std::numeric_limits<double>::infinity()
				&& floor(e2) != e2)
			return NAN;
		else
			return std::pow(e1, e2);
	}
};

template<class Type, bool is_logical>
struct eq {
	static std::string get_name() {
		return "==";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 == e2;
	}
};

template<class Type, bool is_logical>
struct neq {
	static std::string get_name() {
		return "!=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 != e2;
	}
};

template<class Type, bool is_logical>
struct gt {
	static std::string get_name() {
		return ">";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 > e2;
	}
};

template<class Type, bool is_logical>
struct ge {
	static std::string get_name() {
		return ">=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 >= e2;
	}
};

template<class Type, bool is_logical>
struct lt {
	static std::string get_name() {
		return "<";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 < e2;
	}
};

template<class Type, bool is_logical>
struct le {
	static std::string get_name() {
		return "<=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 <= e2;
	}
};

template<>
struct eq<double, false> {
	static std::string get_name() {
		return "==";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 == e2;
	}
};

template<>
struct neq<double, false> {
	static std::string get_name() {
		return "!=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 != e2;
	}
};

template<>
struct gt<double, false> {
	static std::string get_name() {
		return ">";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 > e2;
	}
};

template<>
struct ge<double, false> {
	static std::string get_name() {
		return ">=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 >= e2;
	}
};

template<>
struct lt<double, false> {
	static std::string get_name() {
		return "<";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 < e2;
	}
};

template<>
struct le<double, false> {
	static std::string get_name() {
		return "<=";
	}
	static int get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const double &e1, const double &e2) const {
		if (std::isnan(e1) || std::isnan(e2))
			return R_get_na<int, true>();
		else
			return e1 <= e2;
	}
};

template<class Type, bool is_logical>
struct logic_or {
	static std::string get_name() {
		return "||";
	}
	static int get_agg_init() {
		return false;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 || e2;
	}
};

template<class Type, bool is_logical>
struct logic_and {
	static std::string get_name() {
		return "&&";
	}
	static int get_agg_init() {
		return true;
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e1, const Type &e2) const {
		return e1 && e2;
	}
};

//////////////////////// unary operators /////////////////////////

template<class Type, bool is_logical>
struct uop_neg {
	static std::string get_name() {
		return "-";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return -e;
	}
};

template<class Type, bool is_logical>
struct uop_sqrt {
	static std::string get_name() {
		return "sqrt";
	}
	static R_type get_output_type() {
		return get_Rtype<double, false>();
	}
	double operator()(const Type &e) const {
		return std::sqrt(e);
	}
};

template<class Type, bool is_logical>
struct uop_abs {
	static std::string get_name() {
		return "abs";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return std::abs(e);
	}
};

template<class Type, bool is_logical>
struct uop_not {
	static std::string get_name() {
		return "not";
	}
	static R_type get_output_type() {
		return get_Rtype<int, true>();
	}
	int operator()(const Type &e) const {
		return !e;
	}
};

template<class Type, bool is_logical>
struct sq {
	static std::string get_name() {
		return "sq";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return e * e;
	}
};

template<class Type, bool is_logical>
struct ceil {
	static std::string get_name() {
		return "ceil";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return std::ceil(e);
	}
};

template<class Type, bool is_logical>
struct floor {
	static std::string get_name() {
		return "floor";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return std::floor(e);
	}
};

template<class Type, bool is_logical>
struct round {
	static std::string get_name() {
		return "round";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	Type operator()(const Type &e) const {
		return std::round(e);
	}
};

template<class Type, bool is_logical>
struct log {
	static std::string get_name() {
		return "log";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	double operator()(const Type &e) const {
		return std::log((double) e);
	}
};

template<class Type, bool is_logical>
struct log2 {
	static std::string get_name() {
		return "log2";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	double operator()(const Type &e) const {
		return std::log2((double) e);
	}
};

template<class Type, bool is_logical>
struct log10 {
	static std::string get_name() {
		return "log10";
	}
	static R_type get_output_type() {
		return get_Rtype<Type, false>();
	}
	double operator()(const Type &e) const {
		return std::log10((double) e);
	}
};

////////////////////// Define basic Ops //////////////////////

class basic_Ruops: public basic_uops
{
public:
	typedef std::shared_ptr<basic_Ruops> ptr;
	virtual R_type get_output_type(op_idx idx) const = 0;
};

class basic_Rops: public basic_ops
{
public:
	typedef std::shared_ptr<basic_Rops> ptr;
	virtual R_type get_output_type(op_idx idx) const = 0;
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class Type, bool is_logical>
class basic_Ruops_impl: public basic_Ruops
{
	bulk_uoperate_impl<uop_neg<Type, is_logical>, Type, Type> neg_op;
	bulk_uoperate_impl<uop_sqrt<Type, is_logical>, Type, double> sqrt_op;
	bulk_uoperate_impl<uop_abs<Type, is_logical>, Type, Type> abs_op;
	bulk_uoperate_impl<uop_not<Type, is_logical>, Type, bool> not_op;
	bulk_uoperate_impl<sq<Type, is_logical>, Type, Type> sq_op;
	bulk_uoperate_impl<ceil<Type, is_logical>, Type, Type> ceil_op;
	bulk_uoperate_impl<floor<Type, is_logical>, Type, Type> floor_op;
	bulk_uoperate_impl<round<Type, is_logical>, Type, Type> round_op;
	bulk_uoperate_impl<log<Type, is_logical>, Type, double> log_op;
	bulk_uoperate_impl<log2<Type, is_logical>, Type, double> log2_op;
	bulk_uoperate_impl<log10<Type, is_logical>, Type, double> log10_op;

	std::vector<bulk_uoperate *> ops;
	std::vector<R_type> R_output_types;
public:
	basic_Ruops_impl() {
		ops.resize(NUM_OPS);
		ops[NEG] = &neg_op;
		ops[SQRT] = &sqrt_op;
		ops[ABS] = &abs_op;
		ops[NOT] = &not_op;
		ops[SQ] = &sq_op;
		ops[CEIL] = &ceil_op;
		ops[FLOOR] = &floor_op;
		ops[ROUND] = &round_op;
		ops[LOG] = &log_op;
		ops[LOG2] = &log2_op;
		ops[LOG10] = &log10_op;

		R_output_types.resize(NUM_OPS);
		R_output_types[NEG] = uop_neg<Type, is_logical>::get_output_type();
		R_output_types[SQRT] = uop_sqrt<Type, is_logical>::get_output_type();
		R_output_types[ABS] = uop_abs<Type, is_logical>::get_output_type();
		R_output_types[NOT] = uop_not<Type, is_logical>::get_output_type();
		R_output_types[SQ] = sq<Type, is_logical>::get_output_type();
		R_output_types[CEIL] = ceil<Type, is_logical>::get_output_type();
		R_output_types[FLOOR] = floor<Type, is_logical>::get_output_type();
		R_output_types[ROUND] = round<Type, is_logical>::get_output_type();
		R_output_types[LOG] = log<Type, is_logical>::get_output_type();
		R_output_types[LOG2] = log2<Type, is_logical>::get_output_type();
		R_output_types[LOG10] = log10<Type, is_logical>::get_output_type();
	}

	virtual const bulk_uoperate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}

	virtual R_type get_output_type(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return R_type::R_NTYPES;
		return R_output_types[idx];
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class Type, bool is_logical>
class basic_Rops_impl: public basic_Rops
{
	bulk_operate_impl<add<Type, is_logical>, Type, Type, Type> add_op;
	bulk_operate_impl<sub<Type, is_logical>, Type, Type, Type> sub_op;
	bulk_operate_impl<multiply<Type, is_logical>, Type, Type, Type> mul_op;
	bulk_operate_impl<divide<Type, is_logical>, Type, Type, double> div_op;
	bulk_operate_impl<mod<Type, is_logical>, Type, Type, Type> mod_op;
	bulk_operate_impl<idiv<Type, is_logical>, Type, Type, Type> idiv_op;
	bulk_operate_impl<min<Type, is_logical>, Type, Type, Type> min_op;
	bulk_operate_impl<max<Type, is_logical>, Type, Type, Type> max_op;
	bulk_operate_impl<pow<Type, is_logical>, Type, Type, Type> pow_op;
	bulk_operate_impl<eq<Type, is_logical>, Type, Type, int> eq_op;
	bulk_operate_impl<neq<Type, is_logical>, Type, Type, int> neq_op;
	bulk_operate_impl<gt<Type, is_logical>, Type, Type, int> gt_op;
	bulk_operate_impl<ge<Type, is_logical>, Type, Type, int> ge_op;
	bulk_operate_impl<lt<Type, is_logical>, Type, Type, int> lt_op;
	bulk_operate_impl<le<Type, is_logical>, Type, Type, int> le_op;
	bulk_operate_impl<logic_or<Type, is_logical>, Type, Type, int> or_op;
	bulk_operate_impl<logic_and<Type, is_logical>, Type, Type, int> and_op;

	std::vector<R_type> R_output_types;
	std::vector<bulk_operate *> ops;
public:
	basic_Rops_impl() {
		ops.resize(NUM_OPS);
		ops[ADD] = &add_op;
		ops[SUB] = &sub_op;
		ops[MUL] = &mul_op;
		ops[DIV] = &div_op;
		ops[MIN] = &min_op;
		ops[MAX] = &max_op;
		ops[POW] = &pow_op;
		ops[EQ] = &eq_op;
		ops[NEQ] = &neq_op;
		ops[GT] = &gt_op;
		ops[GE] = &ge_op;
		ops[LT] = &lt_op;
		ops[LE] = &le_op;
		ops[OR] = &or_op;
		ops[AND] = &and_op;
		ops[MOD] = &mod_op;
		ops[IDIV] = &idiv_op;

		R_output_types.resize(NUM_OPS);
		R_output_types[ADD] = add<Type, is_logical>::get_output_type();
		R_output_types[SUB] = sub<Type, is_logical>::get_output_type();
		R_output_types[MUL] = multiply<Type, is_logical>::get_output_type();
		R_output_types[DIV] = divide<Type, is_logical>::get_output_type();
		R_output_types[MIN] = min<Type, is_logical>::get_output_type();
		R_output_types[MAX] = max<Type, is_logical>::get_output_type();
		R_output_types[POW] = pow<Type, is_logical>::get_output_type();
		R_output_types[EQ] = eq<Type, is_logical>::get_output_type();
		R_output_types[NEQ] = neq<Type, is_logical>::get_output_type();
		R_output_types[GT] = gt<Type, is_logical>::get_output_type();
		R_output_types[GE] = ge<Type, is_logical>::get_output_type();
		R_output_types[LT] = lt<Type, is_logical>::get_output_type();
		R_output_types[LE] = le<Type, is_logical>::get_output_type();
		R_output_types[OR] = logic_or<Type, is_logical>::get_output_type();
		R_output_types[AND] = logic_and<Type, is_logical>::get_output_type();
		R_output_types[MOD] = mod<Type, is_logical>::get_output_type();
		R_output_types[IDIV] = idiv<Type, is_logical>::get_output_type();
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}

	virtual R_type get_output_type(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return R_type::R_NTYPES;
		return R_output_types[idx];
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class Type, bool is_logical>
class basic_Ruops_NA_impl: public basic_Ruops
{
	struct uop_neg_na: public uop_neg<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : -e;
		}
	};

	struct uop_sqrt_na: public uop_sqrt<Type, is_logical> {
		double operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<double, false>() : std::sqrt(e);
		}
	};

	struct uop_abs_na: public uop_abs<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : std::abs(e);
		}
	};

	struct uop_not_na: public uop_not<Type, is_logical> {
		int operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<int, true>() : !e;
		}
	};

	struct sq_na: public sq<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : e * e;
		}
	};

	struct ceil_na: public ceil<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : std::ceil(e);
		}
	};

	struct floor_na: public floor<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : std::floor(e);
		}
	};

	struct round_na: public round<Type, is_logical> {
		Type operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<Type, false>() : std::round(e);
		}
	};

	struct log_na: public log<Type, is_logical> {
		double operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<double, false>() : std::log((double) e);
		}
	};

	struct log2_na: public log2<Type, is_logical> {
		double operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<double, false>() : std::log2((double) e);
		}
	};

	struct log10_na: public log10<Type, is_logical> {
		double operator()(const Type &e) const {
			return R_is_na<Type, is_logical>(e)
				? R_get_na<double, false>() : std::log10((double) e);
		}
	};

	bulk_uoperate_impl<uop_neg_na, Type, Type> neg_op;
	bulk_uoperate_impl<uop_sqrt_na, Type, double> sqrt_op;
	bulk_uoperate_impl<uop_abs_na, Type, Type> abs_op;
	bulk_uoperate_impl<uop_not_na, Type, int> not_op;
	bulk_uoperate_impl<sq_na, Type, Type> sq_op;
	bulk_uoperate_impl<ceil_na, Type, Type> ceil_op;
	bulk_uoperate_impl<floor_na, Type, Type> floor_op;
	bulk_uoperate_impl<round_na, Type, Type> round_op;
	bulk_uoperate_impl<log_na, Type, double> log_op;
	bulk_uoperate_impl<log2_na, Type, double> log2_op;
	bulk_uoperate_impl<log10_na, Type, double> log10_op;

	std::vector<bulk_uoperate *> ops;
	std::vector<R_type> R_output_types;
public:
	basic_Ruops_NA_impl() {
		ops.resize(NUM_OPS);
		ops[NEG] = &neg_op;
		ops[SQRT] = &sqrt_op;
		ops[ABS] = &abs_op;
		ops[NOT] = &not_op;
		ops[SQ] = &sq_op;
		ops[CEIL] = &ceil_op;
		ops[FLOOR] = &floor_op;
		ops[ROUND] = &round_op;
		ops[LOG] = &log_op;
		ops[LOG2] = &log2_op;
		ops[LOG10] = &log10_op;

		R_output_types.resize(NUM_OPS);
		R_output_types[NEG] = uop_neg<Type, is_logical>::get_output_type();
		R_output_types[SQRT] = uop_sqrt<Type, is_logical>::get_output_type();
		R_output_types[ABS] = uop_abs<Type, is_logical>::get_output_type();
		R_output_types[NOT] = uop_not<Type, is_logical>::get_output_type();
		R_output_types[SQ] = sq<Type, is_logical>::get_output_type();
		R_output_types[CEIL] = ceil<Type, is_logical>::get_output_type();
		R_output_types[FLOOR] = floor<Type, is_logical>::get_output_type();
		R_output_types[ROUND] = round<Type, is_logical>::get_output_type();
		R_output_types[LOG] = log<Type, is_logical>::get_output_type();
		R_output_types[LOG2] = log2<Type, is_logical>::get_output_type();
		R_output_types[LOG10] = log10<Type, is_logical>::get_output_type();
	}

	virtual const bulk_uoperate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}

	virtual R_type get_output_type(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return R_type::R_NTYPES;
		return R_output_types[idx];
	}
};

template<class Type, bool is_logical>
struct min_na: public min<Type, is_logical>
{
	Type operator()(const Type &e1, const Type &e2) const {
		return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
			? R_get_na<Type, is_logical>() : std::min(e1, e2);
	}
};

template<>
struct min_na<double, false>: public min<double, false>
{
	double operator()(const double &e1, const double &e2) const {
		if (R_is_na<double, false>(e1) || R_is_na<double, false>(e2))
			return R_get_na<double, false>();
		else if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::min(e1, e2);
	}
};

template<class Type, bool is_logical>
struct max_na: public max<Type, is_logical>
{
	Type operator()(const Type &e1, const Type &e2) const {
		return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
			? R_get_na<Type, is_logical>() : std::max(e1, e2);
	}
};

template<>
struct max_na<double, false>: public max<double, false>
{
	double operator()(const double &e1, const double &e2) const {
		if (R_is_na<double, false>(e1) || R_is_na<double, false>(e2))
			return R_get_na<double, false>();
		else if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::max(e1, e2);
	}
};

template<class Type, bool is_logical>
struct mod_na: public mod<Type, is_logical>
{
	Type operator()(const Type &e1, const Type &e2) const {
		return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
			? R_get_na<Type, is_logical>() : e1 % e2;
	}
};

template<>
struct mod_na<double, false>: public mod<double, false>
{
	double operator()(const double &e1, const double &e2) const {
		return R_is_na<double, false>(e1) || R_is_na<double, false>(e2)
			? R_get_na<double, false>() : std::fmod(e1, e2);
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class Type, bool is_logical>
class basic_Rops_NA_impl: public basic_Rops
{
	struct add_na: public add<Type, is_logical> {
		Type operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				// This operation on a logical outputs an integer
				? R_get_na<Type, false>() : e1 + e2;
		}
	};

	struct sub_na: public sub<Type, is_logical> {
		Type operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<Type, false>() : e1 - e2;
		}
	};

	struct multiply_na: public multiply<Type, is_logical> {
		Type operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<Type, false>() : e1 * e2;
		}
	};

	// Division is special. Its output should be float point.
	// Therefore, we convert both input values to float point.
	struct divide_na: public divide<Type, is_logical> {
		double operator()(const Type &e1, const Type &e2) const {
			double d1 = e1;
			double d2 = e2;
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<double, false>() : d1 / d2;
		}
	};
	struct idiv_na: public idiv<Type, is_logical> {
		divide<Type, is_logical> div;
		Type operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<Type, false>() : std::floor(div(e1, e2));
		}
	};

	struct pow_na: public pow<Type, is_logical> {
		Type operator()(const Type &e1, const Type &e2) const {
			if (e1 == 1 || e2 == 0)
				return 1;
			else if (R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2))
				return R_get_na<Type, false>();
			else
				return pow<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct eq_na: public eq<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : eq<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct neq_na: public neq<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : neq<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct gt_na: public gt<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : gt<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct ge_na: public ge<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : ge<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct lt_na: public lt<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : lt<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct le_na: public le<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			return R_is_na<Type, is_logical>(e1) || R_is_na<Type, is_logical>(e2)
				? R_get_na<int, true>() : le<Type, is_logical>::operator()(e1, e2);
		}
	};

	struct logic_or_na: public logic_or<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			bool is_na1 = R_is_na<Type, is_logical>(e1);
			bool is_na2 = R_is_na<Type, is_logical>(e2);
			if (is_na1 && is_na2)
				return R_get_na<int, true>();
			else if ((is_na1 && e2) || (is_na2 && e1))
				return true;
			else if (is_na1 || is_na2)
				return R_get_na<int, true>();
			else
				return e1 || e2;
		}
	};

	struct logic_and_na: public logic_and<Type, is_logical> {
		int operator()(const Type &e1, const Type &e2) const {
			bool is_na1 = R_is_na<Type, is_logical>(e1);
			bool is_na2 = R_is_na<Type, is_logical>(e2);
			if (is_na1 && is_na2)
				return R_get_na<int, true>();
			else if ((is_na1 && !e2) || (is_na2 && !e1))
				return false;
			else if (is_na1 || is_na2)
				return R_get_na<int, true>();
			else
				return e1 && e2;
		}
	};

	bulk_operate_impl<add_na, Type, Type, Type> add_op;
	bulk_operate_impl<sub_na, Type, Type, Type> sub_op;
	bulk_operate_impl<multiply_na, Type, Type, Type> mul_op;
	bulk_operate_impl<divide_na, Type, Type, double> div_op;
	bulk_operate_impl<mod_na<Type, is_logical>, Type, Type, Type> mod_op;
	bulk_operate_impl<idiv_na, Type, Type, Type> idiv_op;
	bulk_operate_impl<min_na<Type, is_logical>, Type, Type, Type> min_op;
	bulk_operate_impl<max_na<Type, is_logical>, Type, Type, Type> max_op;
	bulk_operate_impl<pow_na, Type, Type, Type> pow_op;
	bulk_operate_impl<eq_na, Type, Type, int> eq_op;
	bulk_operate_impl<neq_na, Type, Type, int> neq_op;
	bulk_operate_impl<gt_na, Type, Type, int> gt_op;
	bulk_operate_impl<ge_na, Type, Type, int> ge_op;
	bulk_operate_impl<lt_na, Type, Type, int> lt_op;
	bulk_operate_impl<le_na, Type, Type, int> le_op;
	bulk_operate_impl<logic_or_na, Type, Type, int> or_op;
	bulk_operate_impl<logic_and_na, Type, Type, int> and_op;

	std::vector<bulk_operate *> ops;
	std::vector<R_type> R_output_types;
public:
	basic_Rops_NA_impl() {
		ops.resize(NUM_OPS);
		ops[ADD] = &add_op;
		ops[SUB] = &sub_op;
		ops[MUL] = &mul_op;
		ops[DIV] = &div_op;
		ops[MIN] = &min_op;
		ops[MAX] = &max_op;
		ops[POW] = &pow_op;
		ops[EQ] = &eq_op;
		ops[NEQ] = &neq_op;
		ops[GT] = &gt_op;
		ops[GE] = &ge_op;
		ops[LT] = &lt_op;
		ops[LE] = &le_op;
		ops[OR] = &or_op;
		ops[AND] = &and_op;
		ops[MOD] = &mod_op;
		ops[IDIV] = &idiv_op;

		R_output_types.resize(NUM_OPS);
		R_output_types[ADD] = add<Type, is_logical>::get_output_type();
		R_output_types[SUB] = sub<Type, is_logical>::get_output_type();
		R_output_types[MUL] = multiply<Type, is_logical>::get_output_type();
		R_output_types[DIV] = divide<Type, is_logical>::get_output_type();
		R_output_types[MIN] = min<Type, is_logical>::get_output_type();
		R_output_types[MAX] = max<Type, is_logical>::get_output_type();
		R_output_types[POW] = pow<Type, is_logical>::get_output_type();
		R_output_types[EQ] = eq<Type, is_logical>::get_output_type();
		R_output_types[NEQ] = neq<Type, is_logical>::get_output_type();
		R_output_types[GT] = gt<Type, is_logical>::get_output_type();
		R_output_types[GE] = ge<Type, is_logical>::get_output_type();
		R_output_types[LT] = lt<Type, is_logical>::get_output_type();
		R_output_types[LE] = le<Type, is_logical>::get_output_type();
		R_output_types[OR] = logic_or<Type, is_logical>::get_output_type();
		R_output_types[AND] = logic_and<Type, is_logical>::get_output_type();
		R_output_types[MOD] = mod<Type, is_logical>::get_output_type();
		R_output_types[IDIV] = idiv<Type, is_logical>::get_output_type();
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}

	virtual R_type get_output_type(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return R_type::R_NTYPES;
		return R_output_types[idx];
	}
};

class generic_bulk_operate
{
	std::string name;
	// a bulk_operate for a different type.
	std::vector<bulk_operate::const_ptr> ops;
public:
	generic_bulk_operate(const std::string &name,
			const std::vector<bulk_operate::const_ptr> &ops) {
		this->name = name;
		this->ops = ops;
	}

	bulk_operate::const_ptr get_op(R_type type) const {
		size_t idx = type;
		if (ops.size() <= idx)
			return bulk_operate::const_ptr();
		else
			return ops[idx];
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
	generic_bulk_uoperate(const std::string &name,
			const std::vector<bulk_uoperate::const_ptr> &ops) {
		this->name = name;
		this->ops = ops;
	}

	bulk_uoperate::const_ptr get_op(R_type type) const {
		size_t idx = type;
		if (ops.size() <= idx)
			return bulk_uoperate::const_ptr();
		else
			return ops[idx];
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
	bulk_ops.emplace_back(name, ops);
}

/*
 * Register a unary UDF.
 * A user has to provide UDFs for all different types.
 */
void register_udf(const std::vector<bulk_uoperate::const_ptr> &ops,
		const std::string &name)
{
	bulk_uops.emplace_back(name, ops);
}

static bool use_na_op = true;

void set_use_na_op(bool val)
{
	use_na_op = val;
}

std::vector<basic_Rops::ptr> bops((int) R_type::R_NTYPES);
std::vector<basic_Rops::ptr> bops_na((int) R_type::R_NTYPES);
std::vector<basic_Ruops::ptr> buops((int) R_type::R_NTYPES);
std::vector<basic_Ruops::ptr> buops_na((int) R_type::R_NTYPES);

static R_type get_op_output_type(basic_ops::op_idx bo_idx, R_type in_type)
{
	if (bo_idx < basic_ops::op_idx::NUM_OPS) {
		if (bops.size() <= (size_t) in_type) {
			fprintf(stderr, "wrong input type\n");
			return R_type::R_NTYPES;
		}

		basic_Rops::ptr ops = bops[(int) in_type];
		return ops->get_output_type(bo_idx);
	}
	else if ((size_t) (bo_idx - basic_ops::op_idx::NUM_OPS) < bulk_ops.size()) {
		size_t off = bo_idx - basic_ops::op_idx::NUM_OPS;
		bulk_operate::const_ptr op = bulk_ops[off].get_op(in_type);
		if (op == NULL) {
			fprintf(stderr,
					"can't find the specified operation with the right type\n");
			return R_type::R_NTYPES;
		}
		return trans_FM2R(op->get_output_type());
	}
	else {
		fprintf(stderr, "can't find the specified operation\n");
		return R_type::R_NTYPES;
	}
}

static bulk_operate::const_ptr _get_op(basic_ops::op_idx bo_idx, int noperands,
		R_type type)
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
		if (bops.size() <= (size_t) type) {
			fprintf(stderr, "wrong type\n");
			return bulk_operate::const_ptr();
		}

		basic_Rops::ptr ops;
		if (use_na_op)
			ops = bops_na[(int) type];
		else
			ops = bops[(int) type];
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
std::pair<bulk_operate::const_ptr, R_type> get_op(SEXP pfun, R_type type)
{
	Rcpp::S4 fun_obj(pfun);
	Rcpp::IntegerVector info = fun_obj.slot("info");
	auto op = _get_op((basic_ops::op_idx) info[0], info[1], type);
	if (op == NULL)
		return std::pair<bulk_operate::const_ptr, R_type>(NULL,
				R_type::R_NTYPES);
	else
		return std::pair<bulk_operate::const_ptr, R_type>(op,
				get_op_output_type((basic_ops::op_idx) info[0], type));
}

/*
 * This construct an aggregation operator from binary operators.
 */
std::pair<agg_operate::const_ptr, R_type> get_agg_op(SEXP pfun, R_type type)
{
	Rcpp::S4 sym_op(pfun);
	Rcpp::IntegerVector agg_info = sym_op.slot("agg");
	basic_ops::op_idx agg_op_idx = (basic_ops::op_idx) agg_info[0];
	bulk_operate::const_ptr agg_op = _get_op(agg_op_idx, agg_info[1], type);
	if (agg_op == NULL)
		return std::pair<agg_operate::const_ptr, R_type>(NULL, R_type::R_NTYPES);
	R_type out_type = get_op_output_type(agg_op_idx, type);

	Rcpp::IntegerVector combine_info = sym_op.slot("combine");
	bulk_operate::const_ptr combine_op;
	if (combine_info[0] >= 0) {
		basic_ops::op_idx comb_op_idx = (basic_ops::op_idx) combine_info[0];
		combine_op = _get_op(comb_op_idx, combine_info[1], out_type);
		if (combine_op == NULL)
			return std::pair<agg_operate::const_ptr, R_type>(NULL,
					R_type::R_NTYPES);
		out_type = get_op_output_type(comb_op_idx, out_type);
	}
	auto ret = agg_operate::create(agg_op, combine_op);
	return std::pair<agg_operate::const_ptr, R_type>(ret, out_type);
}

/*
 * Get a unary operator.
 */
std::pair<bulk_uoperate::const_ptr, R_type> get_uop(SEXP pfun, R_type type)
{
	Rcpp::S4 fun_obj(pfun);
	Rcpp::IntegerVector info = fun_obj.slot("info");
	basic_uops::op_idx bo_idx = (basic_uops::op_idx) info[0];
	int noperands = info[1];
	if (noperands != 1) {
		fprintf(stderr, "This isn't a unary operator\n");
		return std::pair<bulk_uoperate::const_ptr, R_type>(NULL,
				R_type::R_NTYPES);
	}

	R_type output_type = R_type::R_NTYPES;
	bulk_uoperate::const_ptr op;
	if (bo_idx < basic_uops::op_idx::NUM_OPS) {
		if (bops.size() <= (size_t) type) {
			fprintf(stderr, "wrong unary operation type\n");
			return std::pair<bulk_uoperate::const_ptr, R_type>(NULL,
					R_type::R_NTYPES);
		}

		basic_Ruops::ptr ops;
		if (use_na_op)
			ops = buops_na[(int) type];
		else
			ops = buops[(int) type];

		output_type = ops->get_output_type(bo_idx);
		op = bulk_uoperate::conv2ptr(*ops->get_op(bo_idx));
		if (op == NULL) {
			fprintf(stderr, "invalid basic unary operator\n");
			return std::pair<bulk_uoperate::const_ptr, R_type>(NULL,
					R_type::R_NTYPES);
		}
	}
	else if ((size_t) (bo_idx - basic_uops::op_idx::NUM_OPS) < bulk_uops.size()) {
		size_t off = bo_idx - basic_uops::op_idx::NUM_OPS;
		op = bulk_uops[off].get_op(type);
		if (op == NULL) {
			fprintf(stderr,
					"can't find the specified unary operation with right type\n");
			return std::pair<bulk_uoperate::const_ptr, R_type>(NULL,
					R_type::R_NTYPES);
		}
		output_type = trans_FM2R(op->get_output_type());
	}
	else {
		fprintf(stderr, "can't find the specified unary operation\n");
		return std::pair<bulk_uoperate::const_ptr, R_type>(NULL,
				R_type::R_NTYPES);
	}
	return std::pair<bulk_uoperate::const_ptr, R_type>(op, output_type);
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
	else if (name == "mod")
		return basic_ops::op_idx::MOD;
	else if (name == "%%")
		return basic_ops::op_idx::MOD;
	else if (name == "%/%")
		return basic_ops::op_idx::IDIV;
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

	virtual void runAgg(size_t num_eles, const void *in, void *output) const {
		double *t_out = (double *) output;
		t_out[0] = num_eles;
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		// We need to use floating-points. "int" may overflow for large vectors.
		return get_scalar_type<double>();
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

	virtual void runAgg(size_t num_eles, const void *in, void *output) const {
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

	virtual void runAgg(size_t num_eles, const void *in, void *output) const {
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

	virtual void runAgg(size_t num_eles, const void *in, void *output) const {
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

template<class InT, class OutT, bool in_logical, bool out_logical>
class cast_ele
{
public:
	OutT operator()(InT v) const {
		if (R_is_na<InT, in_logical>(v))
			return R_get_na<OutT, out_logical>();
		else
			return v;
	}
};

template<>
class cast_ele<double, int, false, false>
{
public:
	int operator()(double v) const {
		if (R_is_na<double, false>(v) || std::isnan(v))
			return R_get_na<int, false>();
		else
			return v;
	}
};

template<>
class cast_ele<double, int, false, true>
{
public:
	int operator()(double v) const {
		if (R_is_na<double, false>(v) || std::isnan(v))
			return R_get_na<int, true>();
		else
			return v != 0;
	}
};

template<class InT, class OutT, bool in_logical, bool out_logical>
class ele_type_cast: public bulk_uoperate
{
	cast_ele<InT, OutT, in_logical, out_logical> cast;
public:
	virtual void runA(size_t num_eles, const void *in_arr,
			void *out_arr) const {
		OutT *t_out = reinterpret_cast<OutT *>(out_arr);
		const InT *t_in = reinterpret_cast<const InT *>(in_arr);
		for (size_t i = 0; i < num_eles; i++)
			t_out[i] = cast(t_in[i]);
	}
	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<InT>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<OutT>();
	}
	virtual std::string get_name() const {
		return std::string("cast(") + get_scalar_type<InT>().get_name() + ", "
			+ get_scalar_type<OutT>().get_name() + ")";
	}
};

void init_udf_ext()
{
	bops[(int) R_type::R_LOGICAL]
		= basic_Rops::ptr(new basic_Rops_impl<int, true>());
	bops[(int) R_type::R_INT]
		= basic_Rops::ptr(new basic_Rops_impl<int, false>());
	bops[(int) R_type::R_REAL]
		= basic_Rops::ptr(new basic_Rops_impl<double, false>());

	bops_na[(int) R_type::R_LOGICAL]
		= basic_Rops::ptr(new basic_Rops_NA_impl<int, true>());
	bops_na[(int) R_type::R_INT]
		= basic_Rops::ptr(new basic_Rops_NA_impl<int, false>());
	bops_na[(int) R_type::R_REAL]
		= basic_Rops::ptr(new basic_Rops_NA_impl<double, false>());

	buops[(int) R_type::R_LOGICAL]
		= basic_Ruops::ptr(new basic_Ruops_impl<int, true>());
	buops[(int) R_type::R_INT]
		= basic_Ruops::ptr(new basic_Ruops_impl<int, false>());
	buops[(int) R_type::R_REAL]
		= basic_Ruops::ptr(new basic_Ruops_impl<double, false>());

	buops_na[(int) R_type::R_LOGICAL]
		= basic_Ruops::ptr(new basic_Ruops_NA_impl<int, true>());
	buops_na[(int) R_type::R_INT]
		= basic_Ruops::ptr(new basic_Ruops_NA_impl<int, false>());
	buops_na[(int) R_type::R_REAL]
		= basic_Ruops::ptr(new basic_Ruops_NA_impl<double, false>());

	std::vector<bulk_operate::const_ptr> ops(R_type::R_NTYPES);
	// Add count.
	ops[R_type::R_LOGICAL]
		= bulk_operate::const_ptr(new r_count_operate<int>());
	ops[R_type::R_INT]
		= bulk_operate::const_ptr(new r_count_operate<int>());
	ops[R_type::R_REAL]
		= bulk_operate::const_ptr(new r_count_operate<double>());
	register_udf(ops, "count");

	ops[R_type::R_LOGICAL]
		= bulk_operate::const_ptr(new r_which_max_operate<int>());
	ops[R_type::R_INT]
		= bulk_operate::const_ptr(new r_which_max_operate<int>());
	ops[R_type::R_REAL]
		= bulk_operate::const_ptr(new r_which_max_operate<double>());
	register_udf(ops, "which.max");

	ops[R_type::R_LOGICAL]
		= bulk_operate::const_ptr(new r_which_min_operate<int>());
	ops[R_type::R_INT]
		= bulk_operate::const_ptr(new r_which_min_operate<int>());
	ops[R_type::R_REAL]
		= bulk_operate::const_ptr(new r_which_min_operate<double>());
	register_udf(ops, "which.min");

	// TODO does this make sense for logicals?
	ops[R_type::R_LOGICAL]
		= bulk_operate::const_ptr(new r_euclidean_operate<int>());
	ops[R_type::R_INT]
		= bulk_operate::const_ptr(new r_euclidean_operate<int>());
	ops[R_type::R_REAL]
		= bulk_operate::const_ptr(new r_euclidean_operate<double>());
	register_udf(ops, "euclidean");

	std::vector<bulk_uoperate::const_ptr> uops(R_type::R_NTYPES);
	uops[R_type::R_LOGICAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, int, true, false>());
	uops[R_type::R_INT] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, int, false, false>());
	uops[R_type::R_REAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<double, int, false, false>());
	register_udf(uops, "as.int");

	uops[R_type::R_LOGICAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, double, true, false>());
	uops[R_type::R_INT] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, double, false, false>());
	uops[R_type::R_REAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<double, double, false, false>());
	register_udf(uops, "as.numeric");

	uops[R_type::R_LOGICAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, int, true, true>());
	uops[R_type::R_INT] = bulk_uoperate::const_ptr(
			new ele_type_cast<int, int, false, true>());
	uops[R_type::R_REAL] = bulk_uoperate::const_ptr(
			new ele_type_cast<double, int, false, true>());
	register_udf(uops, "as.logical");
}

typedef std::vector<arr_apply_operate::const_ptr> app_op_vec;

static std::unordered_map<std::string, app_op_vec> apply_ops;

template<class T>
class rank_apply_operate: public arr_apply_operate
{
	typedef std::pair<T, int> indexed_entry;
	std::vector<std::vector<indexed_entry> > bufs;
	struct {
		bool operator()(const indexed_entry &e1, const indexed_entry &e2) const {
			return e1.first < e2.first;
		}
	} entry_less;
public:
	rank_apply_operate() {
		bufs.resize(detail::mem_thread_pool::get_global_num_threads());
	}

	virtual void run(const local_vec_store &in,
			local_vec_store &out) const {
		assert(out.get_length() == in.get_length());
		const T *in_arr = reinterpret_cast<const T *>(in.get_raw_arr());
		int *out_arr = reinterpret_cast<int *>(out.get_raw_arr());
		int thread_id = detail::mem_thread_pool::get_curr_thread_id();
		std::vector<std::pair<T, int> > &buf
			= const_cast<rank_apply_operate *>(this)->bufs[thread_id];
		buf.resize(in.get_length());
		for (size_t i = 0; i < in.get_length(); i++) {
			buf[i].first = in_arr[i];
			buf[i].second = i;
		}
		std::sort(buf.begin(), buf.end(), entry_less);
		for (size_t i = 0; i < out.get_length(); i++)
			out_arr[i] = buf[i].second;
	}
	virtual size_t get_num_out_eles(size_t num_input) const {
		return num_input;
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
};

template<class T>
class sort_apply_operate: public arr_apply_operate
{
public:
	virtual void run(const local_vec_store &in,
			local_vec_store &out) const {
		assert(out.get_length() == in.get_length());
		memcpy(out.get_raw_arr(), in.get_raw_arr(),
				in.get_entry_size() * in.get_length());
		T *out_arr = reinterpret_cast<T *>(out.get_raw_arr());
		std::sort(out_arr, out_arr + out.get_length());
	}
	virtual size_t get_num_out_eles(size_t num_input) const {
		return num_input;
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<T>();
	}
};

bool register_apply_op(const std::string &name, const app_op_vec &ops)
{
	auto ret = apply_ops.insert(std::pair<std::string, app_op_vec>(name, ops));
	return ret.second;
}

void init_apply_ops()
{
	app_op_vec ops(R_type::R_NTYPES);

	// For logicals
	ops[R_type::R_LOGICAL]
		= arr_apply_operate::const_ptr(new rank_apply_operate<int>());
	// For integers
	ops[R_type::R_INT]
		= arr_apply_operate::const_ptr(new rank_apply_operate<int>());
	// For floating-points
	ops[R_type::R_REAL]
		= arr_apply_operate::const_ptr(new rank_apply_operate<double>());
	register_apply_op("rank", ops);

	// For logicals
	ops[R_type::R_LOGICAL]
		= arr_apply_operate::const_ptr(new sort_apply_operate<int>());
	// For integers
	ops[R_type::R_INT]
		= arr_apply_operate::const_ptr(new sort_apply_operate<int>());
	// For floating-points
	ops[R_type::R_REAL]
		= arr_apply_operate::const_ptr(new sort_apply_operate<double>());
	register_apply_op("sort", ops);
}

std::pair<arr_apply_operate::const_ptr, R_type> get_apply_op(SEXP pfun, R_type type)
{
	Rcpp::S4 sym_op(pfun);
	std::string name = sym_op.slot("name");

	auto it = apply_ops.find(name);
	if (it == apply_ops.end()) {
		fprintf(stderr, "apply function %s doesn't exist\n", name.c_str());
		return std::pair<arr_apply_operate::const_ptr, R_type>(NULL,
				R_type::R_NTYPES);
	}

	const app_op_vec &vec = it->second;
	if (vec.size() <= (size_t) type) {
		fprintf(stderr, "can't find the right type for %s\n", name.c_str());
		return std::pair<arr_apply_operate::const_ptr, R_type>(NULL,
				R_type::R_NTYPES);
	}
	else {
		auto op = vec[(int) type];
		return std::pair<arr_apply_operate::const_ptr, R_type>(op,
				trans_FM2R(op->get_output_type()));
	}
}

dense_matrix::ptr cast_Rtype(dense_matrix::ptr mat, R_type in_type,
		R_type out_type)
{
	if (in_type == out_type)
		return mat;
	else if (out_type == R_type::R_INT) {
		int op_idx = _get_uop_id("as.int");
		size_t off = op_idx - basic_uops::op_idx::NUM_OPS;
		if (off >= bulk_uops.size()) {
			fprintf(stderr, "Can't cast to int\n");
			return dense_matrix::ptr();
		}
		auto op = bulk_uops[off].get_op(in_type);
		return mat->sapply(op);
	}
	else if (out_type == R_type::R_REAL) {
		int op_idx = _get_uop_id("as.numeric");
		size_t off = op_idx - basic_uops::op_idx::NUM_OPS;
		if (off >= bulk_uops.size()) {
			fprintf(stderr, "Can't cast to floating-points\n");
			return dense_matrix::ptr();
		}
		auto op = bulk_uops[off].get_op(in_type);
		return mat->sapply(op);
	}
	else {
		fprintf(stderr, "can't cast to other types.\n");
		return dense_matrix::ptr();
	}
}

}
