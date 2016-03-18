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
