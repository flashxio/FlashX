#ifndef __RAND_GEN_H__
#define __RAND_GEN_H__

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

#include <memory>
#include <random>

#include "generic_type.h"

namespace fm
{

class rand_gen
{
protected:
	static unsigned int gen_rand_seed();
public:
	typedef std::shared_ptr<rand_gen> ptr;
	template<class T>
	static ptr create(T min, T max);
	template<class T>
	static ptr create(T min, T max, T seed);

	virtual void gen(void *buf, size_t len) = 0;
	virtual const scalar_type &get_type() const = 0;
};

template<class T>
class rand_gen_impl: public rand_gen
{
	// By default, it generates 32-bit integers.
	typedef std::mt19937 base_generator_type;
	// By default, I assume it generates integers.
	typedef std::uniform_int_distribution<T> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(T min, T max, T seed): generator(seed), dist(min, max) {
	}

	void gen(void *buf, size_t len) {
		T *t_buf = (T *) buf;
		for (size_t i = 0; i < len; i++)
			t_buf[i] = dist(generator);
	}

	virtual const scalar_type &get_type() const {
		return get_scalar_type<T>();
	}
};

template<class T>
rand_gen::ptr rand_gen::create(T min, T max)
{
	return rand_gen::ptr(new rand_gen_impl<T>(min, max, gen_rand_seed()));
}

template<class T>
rand_gen::ptr rand_gen::create(T min, T max, T seed)
{
	return rand_gen::ptr(new rand_gen_impl<T>(min, max, seed));
}

template<>
class rand_gen_impl<int64_t>: public rand_gen
{
	typedef std::mt19937_64 base_generator_type;
	typedef std::uniform_int_distribution<int64_t> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(int64_t min, int64_t max,
			int64_t seed): generator(seed), dist(min, max) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<int64_t>();
	}
};

template<>
class rand_gen_impl<float>: public rand_gen
{
	typedef std::mt19937 base_generator_type;
	typedef std::uniform_real_distribution<float> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(float min, float max,
			float seed): generator((int32_t) seed), dist(min, max) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<float>();
	}
};

template<>
class rand_gen_impl<double>: public rand_gen
{
	typedef std::mt19937_64 base_generator_type;
	typedef std::uniform_real_distribution<double> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(double min, double max,
			double seed): generator((int64_t) seed), dist(min, max) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<double>();
	}
};

template<>
class rand_gen_impl<long double>: public rand_gen
{
	typedef std::mt19937_64 base_generator_type;
	typedef std::uniform_real_distribution<double> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(long double min, long double max,
			long double seed): generator((int64_t) seed), dist(min, max) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<long double>();
	}
};

template<>
class rand_gen_impl<bool>: public rand_gen
{
	typedef std::mt19937 base_generator_type;
	typedef std::uniform_int_distribution<int> distribution_type;

	base_generator_type generator;
	distribution_type dist;
public:
	rand_gen_impl(bool min, bool max,
			// Use a better seed for the generator.
			bool seed): generator(gen_rand_seed()), dist(0, 1) {
	}

	void gen(void *buf, size_t len);

	virtual const scalar_type &get_type() const {
		return get_scalar_type<bool>();
	}
};

}

#endif
