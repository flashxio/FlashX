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

#include "rand_gen.h"

namespace fm
{

unsigned int rand_gen::gen_rand_seed()
{
	static std::random_device rd;
	return rd();
}

////////////////////////// Random distribution /////////////////////////////

void urand_gen_impl<int64_t>::gen(void *buf, size_t len)
{
	int64_t *t_buf = (int64_t *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

void urand_gen_impl<float>::gen(void *buf, size_t len)
{
	float *t_buf = (float *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

void urand_gen_impl<double>::gen(void *buf, size_t len)
{
	double *t_buf = (double *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

void urand_gen_impl<long double>::gen(void *buf, size_t len)
{
	long double *t_buf = (long double *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

void urand_gen_impl<bool>::gen(void *buf, size_t len)
{
	bool *t_buf = (bool *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

////////////////////////// Normal distribution /////////////////////////////

template<>
class nrand_gen_impl<double>: public rand_gen
{
	typedef std::mt19937_64 base_generator_type;
	typedef std::normal_distribution<double> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	nrand_gen_impl(double mean, double stddev,
			double seed): generator((int64_t) seed), dist(mean, stddev) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<double>();
	}
};

template<>
class nrand_gen_impl<long double>: public rand_gen
{
	typedef std::mt19937_64 base_generator_type;
	typedef std::normal_distribution<double> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	nrand_gen_impl(long double mean, long double stddev,
			long double seed): generator((int64_t) seed), dist(mean, stddev) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<long double>();
	}
};

void nrand_gen_impl<double>::gen(void *buf, size_t len)
{
	double *t_buf = (double *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

void nrand_gen_impl<long double>::gen(void *buf, size_t len)
{
	long double *t_buf = (long double *) buf;
	for (size_t i = 0; i < len; i++)
		t_buf[i] = dist(generator);
}

}
