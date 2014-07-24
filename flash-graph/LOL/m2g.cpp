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
#include "mnist_io.h"

void print_usage()
{
	fprintf(stderr, "m2g input_file output_file #samples\n");
	fprintf(stderr, "-r num: the number of rows used to generate the graph\n");
	fprintf(stderr, "-t: transpose the matrix\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	int opt;
	int num_opts = 0;
	int num_rows = -1;
	bool transpose = false;
	while ((opt = getopt(argc, argv, "r:t")) != -1) {
		num_opts++;
		switch (opt) {
			case 'r':
				num_rows = atoi(optarg);
				num_opts++;
				break;
			case 't':
				transpose = true;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 2) {
		print_usage();
	}

	const std::string input_file = argv[0];
	const std::string output_file = argv[1];

	FG_eigen_matrix<unsigned char>::ptr m = read_mnist_data(input_file, num_rows);
	if (transpose)
		m = m->transpose();
	printf("There are %ld rows and %ld columns\n", m->get_num_rows(), m->get_num_cols());
	graph::ptr g = m->conv2graph();
	printf("There are %ld vertices and %ld edges in the graph\n",
			g->get_num_vertices(), g->get_num_edges());
	g->dump(output_file + ".index", output_file + ".adj");
}
