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

#include "matrix/FG_sparse_matrix.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"

const int NUM_SAMPLES = 100;

FG_vector<unsigned char>::ptr read_mnist_label(const std::string file, int num_samples);

void train(FG_general_sparse_matrix<unsigned char>::ptr input,
		FG_vector<int>::ptr labels, int k)
{
	typedef std::map<int, FG_vector<double>::ptr> mean_map_t;
	mean_map_t mean;
	input->group_by_mean(*labels, true, mean);
	std::vector<eigen_pair_t> eigens;
	compute_SVD<FG_general_sparse_matrix<unsigned char> >(input, 2 * k, k,
			"LA", "RS", eigens);
	BOOST_FOREACH(eigen_pair_t p, eigens) {
		printf("value: %f, ||vector||: %f\n", p.first, p.second->norm1());
	}

	size_t nrow = input->get_num_cols();
	size_t ncol = mean.size() + k;
	FG_eigen_matrix<double>::ptr qr_matrix
		= FG_eigen_matrix<double>::create(nrow, ncol);
	qr_matrix->resize(nrow, ncol);
	printf("construct QR matrix: %ld, %ld\n", qr_matrix->get_num_rows(),
			qr_matrix->get_num_cols());
	int col_idx = 0;
	BOOST_FOREACH(mean_map_t::value_type v, mean) {
		qr_matrix->set_col(col_idx++, *v.second);
	}
	for (int i = 0; i < k; i++)
		qr_matrix->set_col(col_idx++, *eigens[i].second);
}

void print_usage()
{
	fprintf(stderr, "LOL conf-file train-data train-index train-label\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	graph_conf.print_help();
	params.print_help();
	exit(1);
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	while ((opt = getopt(argc, argv, "c:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 4) {
		print_usage();
		exit(-1);
	}

	const std::string conf_file = argv[0];
	const std::string train_data_file = argv[1];
	const std::string train_index_file = argv[2];
	const std::string train_label_file = argv[3];

	config_map configs(conf_file);
	configs.add_options(confs);

	printf("LOL starts\n");
	FG_graph::ptr graph = FG_graph::create(train_data_file, train_index_file,
			configs);
	FG_general_sparse_matrix<unsigned char>::ptr train_matrix
		= FG_general_sparse_matrix<unsigned char>::create(graph);
	train_matrix->resize(NUM_SAMPLES, train_matrix->get_num_cols());

#if 0
	class print_apply
	{
	public:
		void operator()(vertex_id_t vid, general_get_edge_iter<unsigned char>::iterator &it,
				size_t ncol) const {
			vsize_t idx = 0;
			while (it.has_next()) {
				vertex_id_t nid = it.get_curr_id();
				assert(idx == nid);
				if (nid >= ncol)
					break;
				fprintf(stderr, "%d %d %d\n", vid, nid, it.get_curr_value());
				it.next();
				idx++;
			}
		}
	} apply;
	train_matrix->apply(true, apply);
#endif

	FG_vector<unsigned char>::ptr train_label_tmp = read_mnist_label(train_label_file,
			NUM_SAMPLES);
	FG_vector<int>::ptr train_labels
		= FG_vector<int>::create(train_label_tmp->get_size());
	for (size_t i = 0; i < train_label_tmp->get_size(); i++)
		train_labels->set(i, train_label_tmp->get(i));
	train(train_matrix, train_labels, 10);
}
