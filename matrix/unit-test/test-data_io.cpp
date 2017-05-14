/*
 * Copyright 2017 Open Connectome Project (http://openconnecto.me)
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

#include "native_file.h"

#include "data_io.h"
#include "dense_matrix.h"

using namespace fm;

std::unique_ptr<char []> create_text_file(const std::string file, size_t size)
{
	std::unique_ptr<char []> buf(new char[size]);
	for (size_t i = 0; i < size; i++) {
		buf[i] = 'a' + random() % 26;
		if (random() % 1000 == 0)
			buf[i] = '\n';
	}
	buf[size - 1] = '\n';
	FILE *f = fopen(file.c_str(), "w");
	fwrite(buf.get(), size, 1, f);
	fclose(f);
	return buf;
}

void test_read_file()
{
	printf("test text_io on a regular text file\n");
	std::string file = "/tmp/tmp.txt";
	size_t size = 9999999;
	std::unique_ptr<char []> data = create_text_file(file, size);
	text_io::ptr io = text_io::create(file);
	size_t off = 0;
	while (off < size) {
		size_t wanted_size = random() % 100000;
		size_t read_size = 0;
		if (random() % 2 == 0) {
			auto read_buf = io->read_lines(wanted_size, read_size);
			assert(read_size <= size - off);
			assert(memcmp(read_buf.get(), data.get() + off, read_size) == 0);
			off += read_size;
		}
		else {
			auto read_buf = io->peek(wanted_size, read_size);
			assert(read_size <= size - off);
			assert(memcmp(read_buf.get(), data.get() + off, read_size) == 0);
		}
	}
}

void test_read_matrix()
{
	printf("test reading a dense matrix\n");
	std::string file = "/tmp/tmp.txt";
	dense_matrix::ptr mat = dense_matrix::create_randu<int>(0, 1000,
			999999, 10, matrix_layout_t::L_COL);
	mat = mat->cast_ele_type(get_scalar_type<double>());
	mat = mat->multiply_scalar<double>(0.001);
	FILE *f = fopen(file.c_str(), "w");
	mat->print(f);
	fclose(f);

	std::vector<std::string> files(1, file);
	dense_matrix::ptr read_mat = read_matrix(files, true, true, "D", "auto");
	assert(read_mat->get_num_rows() == mat->get_num_rows());
	assert(read_mat->get_num_cols() == mat->get_num_cols());
	auto diff_sum = read_mat->minus(*mat)->sum();
	printf("%g\n", scalar_variable::get_val<double>(*diff_sum));
	assert(std::abs(scalar_variable::get_val<double>(*diff_sum)) < 1e-6);
}

int main()
{
	test_read_file();
	test_read_matrix();
}
