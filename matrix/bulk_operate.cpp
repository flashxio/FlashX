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
 * This template implements the interface of a bulk unary operator.
 */
template<class OpType, class InType, class OutType>
class bulk_uoperate_impl: public bulk_uoperate
{
	OpType op;
public:
	bulk_uoperate_impl() {
	}

	bulk_uoperate_impl(const OpType &_op): op(_op) {
	}

	virtual void runA(size_t num_eles, const void *in_arr1,
			void *output_arr1) const {
		const InType *in_arr = (const InType *) in_arr1;
		OutType *output_arr = (OutType *) output_arr1;
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(in_arr[i]);
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<InType>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<OutType>();
	}
	virtual std::string get_name() const {
		return OpType::get_name();
	}
};

/*
 * This template implements the interface of the bulk binary operator.
 */
template<class OpType, class LeftType, class RightType, class ResType>
class bulk_operate_impl: public bulk_operate
{
	static const size_t BULK_LEN = 128;
	OpType op;

	void runAA_short(size_t num_eles, const LeftType *left_arr,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}
	void runAA_fixed(const LeftType *left_arr,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(left_arr[i], right_arr[i]);
	}

	void runAE_short(size_t num_eles, const LeftType *left_arr,
			const RightType &entry, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(left_arr[i], entry);
	}
	void runAE_fixed(const LeftType *left_arr, const RightType &entry,
			ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(left_arr[i], entry);
	}

	void runEA_short(size_t num_eles, const LeftType &entry,
			const RightType *right_arr, ResType *output_arr) const {
		for (size_t i = 0; i < num_eles; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}
	void runEA_fixed(const LeftType &entry, const RightType *right_arr,
			ResType *output_arr) const {
		for (size_t i = 0; i < BULK_LEN; i++)
			output_arr[i] = op(entry, right_arr[i]);
	}

	ResType runAgg_short(size_t num_eles, const LeftType *left_arr,
			const ResType &orig) const {
		ResType res = orig;
		for (size_t i = 0; i < num_eles; i++)
			res = op(left_arr[i], res);
		return res;
	}
	ResType runAgg_fixed(const LeftType *left_arr, const ResType &orig) const {
		// GCC cannot optimize reduction well for many operations.
		// Because some optimizations may change the result of the reduction.
		// For example, float-point operations aren't strictly commutative or
		// associative; even though min/max is commutative and associative
		// semantically, GCC cannot discover this property by default.
		// Arranging the operation as below helps GCC apply more optimizations
		// and vectorize operations.
		// NOTE: such arranging will change the result of float-point operations
		// slightly. We may not care because parallelization changes the result
		// anyway.
		const size_t tmp_len = 8;
		ResType tmp[tmp_len];
		for (size_t i = 0; i < tmp_len; i++)
			tmp[i] = OpType::get_agg_init();
		for (size_t i = 0; i < BULK_LEN; i += tmp_len) {
			for (size_t j = 0; j < tmp_len; j++)
				tmp[j] = op(left_arr[i + j], tmp[j]);
		}
		ResType res = orig;
		for (size_t i = 0; i < tmp_len; i++)
			res = op(tmp[i], res);
		return res;
	}
public:
	virtual void runAA(size_t num_eles, const void *left_arr1,
			const void *right_arr1, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		if (num_eles < BULK_LEN)
			runAA_short(num_eles, left_arr, right_arr, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runAA_fixed(left_arr + i, right_arr + i, output_arr + i);
			runAA_short(num_eles - num_eles_align, left_arr + num_eles_align,
					right_arr + num_eles_align, output_arr + num_eles_align);
		}
	}

	virtual void runAE(size_t num_eles, const void *left_arr1,
			const void *right, void *output_arr1) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		RightType entry = *(const RightType *) right;
		if (num_eles < BULK_LEN)
			runAE_short(num_eles, left_arr, entry, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runAE_fixed(left_arr + i, entry, output_arr + i);
			runAE_short(num_eles - num_eles_align, left_arr + num_eles_align,
					entry, output_arr + num_eles_align);
		}
	}

	virtual void runEA(size_t num_eles, const void *left,
			const void *right_arr1, void *output_arr1) const {
		LeftType entry = *(const LeftType *) left;
		const RightType *right_arr = (const RightType *) right_arr1;
		ResType *output_arr = (ResType *) output_arr1;
		if (num_eles < BULK_LEN)
			runEA_short(num_eles, entry, right_arr, output_arr);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				runEA_fixed(entry, right_arr + i, output_arr + i);
			runEA_short(num_eles - num_eles_align, entry,
					right_arr + num_eles_align, output_arr + num_eles_align);
		}
	}

	virtual void runAgg(size_t num_eles, const void *left_arr1,
			void *output) const {
		const LeftType *left_arr = (const LeftType *) left_arr1;
		ResType res = OpType::get_agg_init();
		if (num_eles < BULK_LEN)
			*(ResType *) output = runAgg_short(num_eles, left_arr, res);
		else {
			size_t num_eles_align = ROUND(num_eles, BULK_LEN);
			for (size_t i = 0; i < num_eles_align; i += BULK_LEN)
				res = runAgg_fixed(left_arr + i, res);
			*(ResType *) output = runAgg_short(num_eles - num_eles_align,
					left_arr + num_eles_align, res);
		}
	}

	virtual const scalar_type &get_left_type() const {
		return get_scalar_type<LeftType>();
	}
	virtual const scalar_type &get_right_type() const {
		return get_scalar_type<RightType>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<ResType>();
	}
	virtual std::string get_name() const {
		return OpType::get_name();
	}
};

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
	}

	virtual const bulk_operate *get_op(op_idx idx) const {
		if (idx >= op_idx::NUM_OPS)
			return NULL;
		return ops[idx];
	}
};

template class basic_ops_impl<char, char, char>;
template class basic_ops_impl<short, short, short>;
template class basic_ops_impl<int, int, int>;
template class basic_ops_impl<long, long, long>;
template class basic_ops_impl<float, float, float>;
template class basic_ops_impl<double, double, double>;
template class basic_ops_impl<long double, long double, long double>;
template class basic_ops_impl<bool, bool, bool>;
template class basic_ops_impl<unsigned short, unsigned short, unsigned short>;
template class basic_ops_impl<unsigned int, unsigned int, unsigned int>;
template class basic_ops_impl<unsigned long, unsigned long, unsigned long>;

template class basic_uops_impl<char, char>;
template class basic_uops_impl<short, short>;
template class basic_uops_impl<int, int>;
template class basic_uops_impl<long, long>;
template class basic_uops_impl<float, float>;
template class basic_uops_impl<double, double>;
template class basic_uops_impl<long double, long double>;
template class basic_uops_impl<bool, bool>;
template class basic_uops_impl<unsigned short, unsigned short>;
template class basic_uops_impl<unsigned int, unsigned int>;
template class basic_uops_impl<unsigned long, unsigned long>;

template<class T>
const basic_uops &scalar_type_impl<T>::get_basic_uops() const
{
	static basic_uops_impl<T, T> uops;
	return uops;
}

template<class T>
const basic_ops &scalar_type_impl<T>::get_basic_ops() const
{
	static basic_ops_impl<T, T, T> ops;
	return ops;
}

template<class T>
const agg_ops &scalar_type_impl<T>::get_agg_ops() const
{
	static agg_ops_impl<T, T> aops;
	return aops;
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
