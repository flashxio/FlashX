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

#include <stdio.h>

#include "matrix/FG_dense_matrix.h"
#include "FG_vector.h"
#include "mnist_io.h"

int main(int argc, char *argv[])
{
	if (argc < 5) {
		fprintf(stderr, "print_mnist in_data in_label out_data out_label\n");
		exit(1);
	}

	const std::string in_data = argv[1];
	const std::string in_label = argv[2];
	const std::string out_data = argv[3];
	const std::string out_label = argv[4];

	FILE *f_data = fopen(out_data.c_str(), "w");
	assert(f_data);
	FG_eigen_matrix<unsigned char>::ptr m = read_mnist_data(in_data);
	printf("There are %ld rows and %ld columns\n", m->get_num_rows(), m->get_num_cols());
	for (size_t i = 0; i < m->get_num_rows(); i++) {
		for (size_t j = 0; j < m->get_num_cols(); j++) {
			fprintf(f_data, "%d,", m->get(i, j));
		}
		fprintf(f_data, "\n");
	}
	fclose(f_data);

	FILE *f_label = fopen(out_label.c_str(), "w");
	assert(f_label);
	FG_vector<unsigned char>::ptr labels = read_mnist_label(in_label);
	for (size_t i = 0; i < labels->get_size(); i++)
		fprintf(f_label, "%d\n", labels->get(i));
	fclose(f_label);
}
