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

#include <string>

#include "matrix/FG_dense_matrix.h"
#include "FG_vector.h"

const int NUM_SAMPLES = 100;

FG_row_wise_matrix<unsigned char>::ptr read_mnist_data(const std::string &file,
		int num_samples);

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr, "m2g input_file output_file\n");
		exit(1);
	}

	const std::string input_file = argv[1];
	const std::string output_file = argv[2];

	FG_row_wise_matrix<unsigned char>::ptr m = read_mnist_data(input_file, NUM_SAMPLES);
	printf("There are %ld rows and %ld columns\n", m->get_num_rows(), m->get_num_cols());
	for (size_t i = 0; i < m->get_num_rows(); i++) {
		for (size_t j = 0; j < m->get_num_cols(); j++) {
			fprintf(stderr, "%d,", m->get(i, j));
		}
		fprintf(stderr, "\n");
	}
	exit(1);

	graph::ptr g = m->conv2graph();
	printf("There are %ld vertices and %ld edges in the graph\n",
			g->get_num_vertices(), g->get_num_edges());
	g->dump(output_file + ".index", output_file + ".adj");
}
