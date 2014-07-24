/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
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

#define _BSD_SOURCE
#include <endian.h>

#include <stdio.h>

#include "FG_vector.h"
#include "matrix/FG_dense_matrix.h"

uint32_t read32(FILE *f)
{
	uint32_t v;
	size_t ret = fread(&v, sizeof(v), 1, f);
	assert(ret == 1);
	return be32toh(v);
}

FG_eigen_matrix<unsigned char>::ptr read_mnist_data(const std::string &file,
		int num_samples)
{
	FILE *f = fopen(file.c_str(), "r");
	assert(f);

	uint32_t magic_number = read32(f);
	assert(magic_number == 2051);
	uint32_t num_imgs = read32(f);
	printf("There are %d images\n", num_imgs);
	uint32_t num_rows = read32(f);
	uint32_t num_cols = read32(f);
	printf("An image has %d rows and %d columns\n", num_rows, num_cols);


	int num;
	if (num_samples < 0)
		num = num_imgs;
	else
		num = min(num_samples, num_imgs);
	FG_eigen_matrix<unsigned char>::ptr m = FG_eigen_matrix<unsigned char>::create(
			num, num_rows * num_cols);
	m->resize(num, num_rows * num_cols);
	FG_vector<unsigned char>::ptr v = FG_vector<unsigned char>::create(
			num_rows * num_cols);
	for (int i = 0; i < num; i++) {
		size_t ret = fread(v->get_data(), v->get_size(), 1, f);
		m->set_row(i, *v);
		assert(ret == 1);
	}

	fclose(f);
	return m;
}

FG_vector<unsigned char>::ptr read_mnist_label(const std::string file, int num_samples)
{
	FILE *f = fopen(file.c_str(), "r");
	assert(f);

	uint32_t magic_number = read32(f);
	assert(magic_number == 2049);
	uint32_t num_imgs = read32(f);
	printf("There are %d images\n", num_imgs);

	int num;
	if (num_samples < 0)
		num = num_imgs;
	else
		num = min(num_samples, num_imgs);
	FG_vector<unsigned char>::ptr v = FG_vector<unsigned char>::create(num);
	size_t ret = fread(v->get_data(), v->get_size(), 1, f);
	assert(ret == 1);

	fclose(f);
	return v;
}
