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
#include "log.h"

#include "bulk_operate.h"
#include "bulk_operate_impl.h"
#include "bulk_operate_ext.h"
#include "generic_type.h"
#include "rand_gen.h"

namespace fm
{

agg_operate::const_ptr agg_operate::create(bulk_operate::const_ptr agg,
			bulk_operate::const_ptr combine)
{
	if (agg == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< "the agg operator is required in agg_operate";
		return const_ptr();
	}
	if (combine && agg->get_output_type() != combine->get_left_type()
			&& combine->get_left_type() != combine->get_output_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the agg and combine operators have incompatible types";
		return const_ptr();
	}
	return const_ptr(new agg_operate(agg, combine));
}

agg_operate::const_ptr agg_operate::create(bulk_operate::const_ptr agg)
{
	if (agg->get_left_type() != agg->get_right_type()
			|| agg->get_left_type() != agg->get_output_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The agg operator needs to have the same input and output type";
		return const_ptr();
	}
	return const_ptr(new agg_operate(agg, agg));
}

/*
 * This template implements all basic binary operators for different types.
 */
template<class InType, class OutType>
class basic_uops_impl: public basic_uops
{
	struct uop_neg {
		static std::string get_name() {
			return "-";
		}
		OutType operator()(const InType &e) const {
			return -e;
		}
	};

	struct uop_sqrt {
		static std::string get_name() {
			return "sqrt";
		}
		double operator()(const InType &e) const {
			return std::sqrt(e);
		}
	};

	struct uop_abs {
		static std::string get_name() {
			return "abs";
		}
		OutType operator()(const InType &e) const {
			return (OutType) std::abs(e);
		}
	};

	struct uop_not {
		static std::string get_name() {
			return "not";
		}
		bool operator()(const InType &e) const {
			return !e;
		}
	};

	struct sq {
		static std::string get_name() {
			return "sq";
		}
		OutType operator()(const InType &e) const {
			return e * e;
		}
	};

	struct ceil {
		static std::string get_name() {
			return "ceil";
		}
		OutType operator()(const InType &e) const {
			return std::ceil(e);
		}
	};

	struct floor {
		static std::string get_name() {
			return "floor";
		}
		OutType operator()(const InType &e) const {
			return std::floor(e);
		}
	};

	struct round {
		static std::string get_name() {
			return "round";
		}
		OutType operator()(const InType &e) const {
			return std::round(e);
		}
	};

	struct log {
		static std::string get_name() {
			return "log";
		}
		OutType operator()(const InType &e) const {
			return std::log(e);
		}
	};

	struct log2 {
		static std::string get_name() {
			return "log2";
		}
		OutType operator()(const InType &e) const {
			return std::log2(e);
		}
	};

	struct log10 {
		static std::string get_name() {
			return "log10";
		}
		OutType operator()(const InType &e) const {
			return std::log10(e);
		}
	};

	struct sign {
		static std::string get_name() {
			return "sign";
		}
		OutType operator()(const InType &e) const {
			if (e > 0)
				return 1;
			else if (e < 0)
				return -1;
			else
				return 0;
		}
	};

	bulk_uoperate_impl<uop_neg, InType, OutType> neg_op;
	bulk_uoperate_impl<uop_sqrt, InType, double> sqrt_op;
	bulk_uoperate_impl<uop_abs, InType, OutType> abs_op;
	bulk_uoperate_impl<uop_not, InType, bool> not_op;
	bulk_uoperate_impl<sq, InType, OutType> sq_op;
	bulk_uoperate_impl<ceil, InType, OutType> ceil_op;
	bulk_uoperate_impl<floor, InType, OutType> floor_op;
	bulk_uoperate_impl<round, InType, OutType> round_op;
	bulk_uoperate_impl<log, InType, OutType> log_op;
	bulk_uoperate_impl<log2, InType, OutType> log2_op;
	bulk_uoperate_impl<log10, InType, OutType> log10_op;
	bulk_uoperate_impl<sign, InType, OutType> sign_op;

	std::vector<bulk_uoperate *> ops;
public:
	basic_uops_impl() {
		ops.push_back(&neg_op);
		ops.push_back(&sqrt_op);
		ops.push_back(&abs_op);
		ops.push_back(&not_op);
		ops.push_back(&sq_op);
		ops.push_back(&ceil_op);
		ops.push_back(&floor_op);
		ops.push_back(&round_op);
		ops.push_back(&log_op);
		ops.push_back(&log2_op);
		ops.push_back(&log10_op);
		ops.push_back(&sign_op);
	}

	virtual const bulk_uoperate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}
};

template<class LeftType1, class RightType1, class ResType1>
struct multiply
{
	static std::string get_name() {
		return "*";
	}
	static ResType1 get_agg_init() {
		return 1;
	}
	ResType1 operator()(const LeftType1 &e1, const RightType1 &e2) const {
		return e1 * e2;
	}
};

template<>
struct multiply<double, double, double>
{
	static std::string get_name() {
		return "*";
	}
	static double get_agg_init() {
		return 1;
	}
	double operator()(const double &e1, const double &e2) const {
		long double first = e1;
		long double second = e2;
		return first * second;
	}
};

template<class LeftType, class RightType, class ResType>
struct mod
{
	static std::string get_name() {
		return "%";
	}
	static ResType get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	ResType operator()(const LeftType &e1, const RightType &e2) const {
		return e1 % e2;
	}
};

template<>
struct mod<float, float, float>
{
	static std::string get_name() {
		return "%";
	}
	static float get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	float operator()(const float &e1, const float &e2) const {
		return std::fmod(e1, e2);
	}
};

template<>
struct mod<double, double, double>
{
	static std::string get_name() {
		return "%";
	}
	static double get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	double operator()(const double &e1, const double &e2) const {
		return std::fmod(e1, e2);
	}
};

template<>
struct mod<long double, long double, long double>
{
	static std::string get_name() {
		return "%";
	}
	static long double get_agg_init() {
		// This operation isn't used in aggregation, so we
		// don't care this agg init.
		return 0;
	}
	long double operator()(const long double &e1,
			const long double &e2) const {
		return std::fmod(e1, e2);
	}
};

template<class LeftType, class RightType, class ResType>
struct min
{
	static std::string get_name() {
		return "min";
	}
	static ResType get_agg_init() {
		return std::numeric_limits<ResType>::max();
	}
	ResType operator()(const LeftType &e1, const RightType &e2) const {
		return std::min(e1, e2);
	}
};

template<>
struct min<double, double, double>
{
	static std::string get_name() {
		return "min";
	}
	static double get_agg_init() {
		return std::numeric_limits<double>::max();
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

template<class LeftType, class RightType, class ResType>
struct max
{
	static std::string get_name() {
		return "max";
	}
	static ResType get_agg_init() {
		return std::numeric_limits<ResType>::min();
	}
	ResType operator()(const LeftType &e1, const RightType &e2) const {
		return std::max(e1, e2);
	}
};

template<>
struct max<double, double, double>
{
	static std::string get_name() {
		return "max";
	}
	static double get_agg_init() {
		// We need to define the minimum float-point differently.
		return -std::numeric_limits<double>::max();
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

template<>
struct max<float, float, float>
{
	static std::string get_name() {
		return "max";
	}
	static float get_agg_init() {
		// We need to define the minimum float-point differently.
		return -std::numeric_limits<float>::max();
	}
	float operator()(const float &e1, const float &e2) const {
		if (std::isnan(e1))
			return e1;
		else if (std::isnan(e2))
			return e2;
		else
			return std::max(e1, e2);
	}
};

/*
 * This template implements all basic binary operators for different types.
 */
template<class LeftType, class RightType, class ResType>
class basic_ops_impl: public basic_ops
{
	struct add {
		static std::string get_name() {
			return "+";
		}
		static ResType get_agg_init() {
			return 0;
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 + e2;
		}
	};

	struct sub {
		static std::string get_name() {
			return "-";
		}
		static ResType get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return 0;
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return e1 - e2;
		}
	};

	// Division is special. Its output should be float point.
	// Therefore, we convert both input values to float point.
	struct divide {
		static std::string get_name() {
			return "/";
		}
		static ResType get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return 0;
		}
		double operator()(const LeftType &e1, const RightType &e2) const {
			double d1 = e1;
			double d2 = e2;
			return d1 / d2;
		}
	};
	struct divide_float {
		static std::string get_name() {
			return "/";
		}
		static float get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return 0;
		}
		float operator()(const float &e1, const float &e2) const {
			return e1 / e2;
		}
	};
	struct idiv {
		divide div;
		static std::string get_name() {
			return "%/%";
		}
		static ResType get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return 0;
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return floor(div(e1, e2));
		}
	};

	struct pow {
		static std::string get_name() {
			return "pow";
		}
		static ResType get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return 0;
		}
		ResType operator()(const LeftType &e1, const RightType &e2) const {
			return (ResType) std::pow(e1, e2);
		}
	};

	struct eq {
		static std::string get_name() {
			return "==";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 == e2;
		}
	};

	struct neq {
		static std::string get_name() {
			return "!=";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 != e2;
		}
	};

	struct gt {
		static std::string get_name() {
			return ">";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 > e2;
		}
	};

	struct ge {
		static std::string get_name() {
			return ">=";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 >= e2;
		}
	};

	struct lt {
		static std::string get_name() {
			return "<";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 < e2;
		}
	};

	struct le {
		static std::string get_name() {
			return "<=";
		}
		static bool get_agg_init() {
			// This operation isn't used in aggregation, so we
			// don't care this agg init.
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 <= e2;
		}
	};

	struct logic_or {
		static std::string get_name() {
			return "||";
		}
		static bool get_agg_init() {
			return false;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 || e2;
		}
	};

	struct logic_and {
		static std::string get_name() {
			return "&&";
		}
		static bool get_agg_init() {
			return true;
		}
		bool operator()(const LeftType &e1, const RightType &e2) const {
			return e1 && e2;
		}
	};

	bulk_operate_impl<add, LeftType, RightType, ResType> add_op;
	bulk_operate_impl<sub, LeftType, RightType, ResType> sub_op;
	bulk_operate_impl<multiply<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> mul_op;
	bulk_operate_impl<divide, LeftType, RightType, double> div_op;
	bulk_operate_impl<divide_float, float, float, float> div_float_op;
	bulk_operate_impl<mod<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> mod_op;
	bulk_operate_impl<idiv, LeftType, RightType, ResType> idiv_op;
	bulk_operate_impl<min<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> min_op;
	bulk_operate_impl<max<LeftType, RightType, ResType>,
		LeftType, RightType, ResType> max_op;
	bulk_operate_impl<pow, LeftType, RightType, ResType> pow_op;
	bulk_operate_impl<eq, LeftType, RightType, bool> eq_op;
	bulk_operate_impl<neq, LeftType, RightType, bool> neq_op;
	bulk_operate_impl<gt, LeftType, RightType, bool> gt_op;
	bulk_operate_impl<ge, LeftType, RightType, bool> ge_op;
	bulk_operate_impl<lt, LeftType, RightType, bool> lt_op;
	bulk_operate_impl<le, LeftType, RightType, bool> le_op;
	bulk_operate_impl<logic_or, LeftType, RightType, bool> or_op;
	bulk_operate_impl<logic_and, LeftType, RightType, bool> and_op;

	std::vector<bulk_operate *> ops;
public:
	basic_ops_impl() {
		ops.resize(NUM_OPS);
		ops[ADD] = &add_op;
		ops[SUB] = &sub_op;
		ops[MUL] = &mul_op;
		// We should treat float differently, because we want to output float.
		if (get_type<LeftType>() == prim_type::P_FLOAT
				&& get_type<RightType>() == prim_type::P_FLOAT)
			ops[DIV] = &div_float_op;
		else
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
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}
};

/*
 * Find the first location where its element is different from the first
 * element in the array.
 */
template<class T>
class find_next_impl: public bulk_operate
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

	/*
	 * The first element in the array is pointed by `left_arr' and there are
	 * `num_eles' in the array.
	 * If all elements are the same, it returns the number of elements
	 * in the array.
	 */
	virtual void runAgg(size_t num_eles, const void *left_arr,
			void *output) const {
		const T *curr = (const T *) left_arr;
		T val = *curr;
		size_t loc = 1;
		for (; loc < num_eles && curr[loc] == val; loc++);
		*(size_t *) output = loc;
	}

	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
	virtual std::string get_name() const {
		return "find_next";
	}
};

/*
 * Search backwards and find the first location where its element is different
 * from the last element in the array.
 */
template<class T>
class find_prev_impl: public bulk_operate
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

	/*
	 * The end of the array is indicated by `arr_end'. The last element
	 * is right before `arr_end'. There are `num_eles' elements in the array.
	 * If all elements are the same, it returns the number of elements
	 * in the array.
	 */
	virtual void runAgg(size_t num_eles, const void *arr_end,
			void *output) const {
		const T *curr = ((const T *) arr_end) - 1;
		T val = *curr;
		const T *first = ((const T *) arr_end) - num_eles;
		for (; curr > first && *curr == val; curr--);
		if (*curr != val)
			curr++;
		assert(*curr == val);
		*(size_t *) output = ((const T *) arr_end) - curr;
	}

	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
	virtual std::string get_name() const {
		return "find_prev";
	}
};

template<class T>
class count_operate: public bulk_operate
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
		size_t *t_out = (size_t *) output;
		t_out[0] = num_eles;
	}

	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		throw unsupported_exception();
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
	virtual std::string get_name() const {
		return "count";
	}
};

template<class T>
class argmax_operate: public bulk_operate
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
		t_out[0] = idx;
	}

	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		throw unsupported_exception();
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
		return "argmax";
	}
};

template<class T>
class argmin_operate: public bulk_operate
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
		t_out[0] = idx;
	}

	virtual void runCum(size_t num_eles, const void *left_arr,
			const void *prev, void *output) const {
		throw unsupported_exception();
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
		return "argmin";
	}
};

template<class InType, class OutType>
class agg_ops_impl: public agg_ops
{
public:
	agg_ops_impl() {
		bulk_operate::const_ptr count_agg
			= bulk_operate::const_ptr(new count_operate<InType>());
		ops[COUNT] = agg_operate::create(count_agg,
				bulk_operate::conv2ptr(
					count_agg->get_output_type().get_basic_ops().get_add()));
		ops[FIND_NEXT] = agg_operate::create(
				bulk_operate::const_ptr(new find_next_impl<InType>()),
				bulk_operate::const_ptr());
		ops[FIND_PREV] = agg_operate::create(
				bulk_operate::const_ptr(new find_prev_impl<InType>()),
				bulk_operate::const_ptr());
		ops[ARGMIN] = agg_operate::create(
				bulk_operate::const_ptr(new argmin_operate<InType>()),
				bulk_operate::const_ptr());
		ops[ARGMAX] = agg_operate::create(
				bulk_operate::const_ptr(new argmax_operate<InType>()),
				bulk_operate::const_ptr());

		auto min = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::MIN));
		ops[MIN] = agg_operate::create(min, min);
		auto max = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::MAX));
		ops[MAX] = agg_operate::create(max, max);
		auto add = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::ADD));
		ops[SUM] = agg_operate::create(add, add);
		auto mul = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::MUL));
		ops[PROD] = agg_operate::create(mul, mul);
		auto agg_and = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::AND));
		auto combine_and = bulk_operate::conv2ptr(
				*get_scalar_type<bool>().get_basic_ops().get_op(
					basic_ops::op_idx::AND));
		ops[AND] = agg_operate::create(agg_and, combine_and);
		auto agg_or = bulk_operate::conv2ptr(
				*get_scalar_type<InType>().get_basic_ops().get_op(
					basic_ops::op_idx::OR));
		auto combine_or = bulk_operate::conv2ptr(
				*get_scalar_type<bool>().get_basic_ops().get_op(
					basic_ops::op_idx::OR));
		ops[OR] = agg_operate::create(agg_or, combine_or);
	}
};

#define set_basic_ops_impl(T) do {\
	int idx = (int) fm::get_type<T>(); \
	basic_ops_impls[idx] = basic_ops::ptr(new basic_ops_impl<T, T, T>()); \
} while (0)

#define set_basic_uops_impl(T) do {\
	int idx = (int) fm::get_type<T>(); \
	basic_uops_impls[idx] = basic_uops::ptr(new basic_uops_impl<T, T>()); \
} while(0)

#define set_agg_ops_impl(T) do {   \
	int idx = (int) fm::get_type<T>(); \
	agg_ops_impls[idx] = agg_ops::ptr(new agg_ops_impl<T, T>()); \
} while(0)

void scalar_type::init_ops()
{
	basic_ops_impls.resize((int) prim_type::NUM_TYPES);
	basic_uops_impls.resize((int) prim_type::NUM_TYPES);
	agg_ops_impls.resize((int) prim_type::NUM_TYPES);
	printf("Initialize ops. There are %ld types\n", agg_ops_impls.size());

	set_basic_ops_impl(char);
	set_basic_ops_impl(short);
	set_basic_ops_impl(int);
	set_basic_ops_impl(long);
	set_basic_ops_impl(float);
	set_basic_ops_impl(double);
	set_basic_ops_impl(long double);
	set_basic_ops_impl(bool);
	set_basic_ops_impl(unsigned short);
	set_basic_ops_impl(unsigned int);
	set_basic_ops_impl(unsigned long);
	for (size_t i = 0; i < basic_ops_impls.size(); i++)
		if (basic_ops_impls[i] == NULL)
			throw unsupported_exception("find an unsupported type");

	set_basic_uops_impl(char);
	set_basic_uops_impl(short);
	set_basic_uops_impl(int);
	set_basic_uops_impl(long);
	set_basic_uops_impl(float);
	set_basic_uops_impl(double);
	set_basic_uops_impl(long double);
	set_basic_uops_impl(bool);
	set_basic_uops_impl(unsigned short);
	set_basic_uops_impl(unsigned int);
	set_basic_uops_impl(unsigned long);
	for (size_t i = 0; i < basic_uops_impls.size(); i++)
		if (basic_uops_impls[i] == NULL)
			throw unsupported_exception("find an unsupported type");

	set_agg_ops_impl(char);
	set_agg_ops_impl(short);
	set_agg_ops_impl(int);
	set_agg_ops_impl(long);
	set_agg_ops_impl(float);
	set_agg_ops_impl(double);
	set_agg_ops_impl(long double);
	set_agg_ops_impl(bool);
	set_agg_ops_impl(unsigned short);
	set_agg_ops_impl(unsigned int);
	set_agg_ops_impl(unsigned long);
	for (size_t i = 0; i < agg_ops_impls.size(); i++)
		if (agg_ops_impls[i] == NULL)
			throw unsupported_exception("find an unsupported type");
}

namespace
{

/*
 * This class set elements in a container randomly.
 * set_operate can't change its own state and has to be thread-safe when
 * running on multiple threads. However, random generators aren't
 * thread-safe, so we have to create a random generator for each thread.
 */
class rand_init: public set_operate
{
public:
	enum rand_dist_type {
		NORM,
		UNIF,
		MAX_NUM,
	};
private:
	class rand_gen_wrapper {
		rand_gen::ptr gen;
	public:
		rand_gen_wrapper(rand_gen::ptr gen) {
			this->gen = gen;
		}

		rand_gen &get_gen() {
			return *gen;
		}
	};

	pthread_key_t gen_key;
	const scalar_type &type;
	scalar_variable::const_ptr var1;
	scalar_variable::const_ptr var2;
	rand_dist_type rand_dist;

	rand_gen &get_rand_gen() const {
		void *addr = pthread_getspecific(gen_key);
		if (addr == NULL) {
			if (rand_dist == rand_dist_type::NORM)
				addr = new rand_gen_wrapper(type.create_randn_gen(*var1, *var2));
			else if (rand_dist == rand_dist_type::UNIF)
				addr = new rand_gen_wrapper(type.create_randu_gen(*var1, *var2));
			else
				assert(0);
			int ret = pthread_setspecific(gen_key, addr);
			assert(ret == 0);
		}
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) addr;
		return wrapper->get_gen();
	}

	static void destroy_rand_gen(void *gen) {
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) gen;
		delete wrapper;
		printf("destroy rand gen\n");
	}
public:
	rand_init(scalar_variable::const_ptr var1, scalar_variable::const_ptr var2,
			rand_dist_type rand_dist): type(var1->get_type()) {
		this->var1 = var1;
		this->var2 = var2;
		int ret = pthread_key_create(&gen_key, destroy_rand_gen);
		this->rand_dist = rand_dist;
		assert(ret == 0);
	}

	~rand_init() {
		pthread_key_delete(gen_key);
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		get_rand_gen().gen(arr, num_eles);
	}
	virtual const scalar_type &get_type() const {
		return get_rand_gen().get_type();
	}
	virtual set_operate::const_ptr transpose() const {
		return set_operate::const_ptr();
	}
};

}

set_operate::const_ptr create_urand_init(scalar_variable::const_ptr min,
		scalar_variable::const_ptr max)
{
	assert(min->get_type() == max->get_type());
	return set_operate::const_ptr(new rand_init(min, max,
				rand_init::rand_dist_type::UNIF));
}

set_operate::const_ptr create_nrand_init(
		scalar_variable::const_ptr mean, scalar_variable::const_ptr var)
{
	assert(mean->get_type() == var->get_type());
	return set_operate::const_ptr(new rand_init(mean, var,
				rand_init::rand_dist_type::NORM));
}

}
