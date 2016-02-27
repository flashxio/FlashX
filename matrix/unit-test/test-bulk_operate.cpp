#include "bulk_operate.h"

using namespace fm;

template<class T>
void test_add(size_t num_eles)
{
	std::vector<T> vec(num_eles);
	for (size_t i = 0; i < vec.size(); i++)
		vec[i] = i;

	const scalar_type &type = get_scalar_type<T>();
	const bulk_operate &op = type.get_basic_ops().get_add();

	// Test agg
	T res;
	op.runAgg(vec.size(), vec.data(), &res);
	assert(vec.size() * (vec.size() - 1) / 2 == res);

	// Test runAA
	std::vector<T> res_vec(num_eles);
	op.runAA(vec.size(), vec.data(), vec.data(), res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(vec[i] * 2 == res_vec[i]);

	// Test runAE
	T val = 2;
	op.runAE(vec.size(), vec.data(), &val, res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(vec[i] + val == res_vec[i]);

	// Test runEA
	op.runEA(vec.size(), &val, vec.data(), res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(vec[i] + val == res_vec[i]);
}

template<class T>
void test_sub(size_t num_eles)
{
	std::vector<T> vec(num_eles);
	for (size_t i = 0; i < vec.size(); i++)
		vec[i] = i;

	const scalar_type &type = get_scalar_type<T>();
	const bulk_operate &op = type.get_basic_ops().get_sub();

	// Test runAA
	std::vector<T> res_vec(num_eles);
	op.runAA(vec.size(), vec.data(), vec.data(), res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(0 == res_vec[i]);

	// Test runAE
	T val = 2;
	op.runAE(vec.size(), vec.data(), &val, res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(vec[i] - val == res_vec[i]);

	// Test runEA
	op.runEA(vec.size(), &val, vec.data(), res_vec.data());
	for (size_t i = 0; i < res_vec.size(); i++)
		assert(val - vec[i] == res_vec[i]);
}

int main()
{
	test_add<int>(100);
	test_add<int>(10000);
	test_sub<int>(100);
	test_sub<int>(10000);

	test_add<double>(100);
	test_add<double>(10000);
	test_sub<double>(100);
	test_sub<double>(10000);
}
