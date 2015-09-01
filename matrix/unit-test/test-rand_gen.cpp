#include <assert.h>

#include <iostream>
#include <memory>

#include "rand_gen.h"

using namespace fm;

int num_tries = 1000000;

template<class T>
void test_rand(T min, T max)
{
	rand_gen::ptr gen = rand_gen::create_randu<T>(min, max);
	assert(gen->get_type() == get_scalar_type<T>());
	std::vector<T> randv(num_tries);
	gen->gen((char *) randv.data(), randv.size());
	for (size_t i = 0; i < randv.size(); i++)
		assert(randv[i] >= min && randv[i] <= max);
}

void test_rand_bool()
{
	rand_gen::ptr gen = rand_gen::create_randu<bool>(0, 1);
	assert(gen->get_type() == get_scalar_type<bool>());
	std::unique_ptr<bool[]> randv(new bool[num_tries]);
	gen->gen((char *) randv.get(), num_tries);
	size_t num_trues = 0;
	for (int i = 0; i < num_tries; i++)
		num_trues += randv[i];
	printf("There are %ld true values out of %d values\n", num_trues, num_tries);
}

int main()
{
	test_rand<int>(1, 10);
	test_rand<long>(1L << 40, 10L << 40);
	test_rand<float>(1, 10);
	test_rand<double>(1, 10);
	test_rand_bool();
}
