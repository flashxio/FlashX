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

#include <shogun/labels/MulticlassLabels.h>
#include <shogun/multiclass/MCLDA.h>
#include <shogun/features/DenseFeatures.h>
#include <shogun/labels/Labels.h>
#include <shogun/base/init.h>
#include <shogun/multiclass/QDA.h>

#include "matrix/FG_sparse_matrix.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/matrix_eigensolver.h"
#include "mnist_io.h"

const int NUM_SAMPLES = 100;
typedef FG_sparse_matrix<general_get_edge_iter<unsigned char> > FG_general_sparse_char_matrix;

/**
 * This procedure generates a dense projection matrix of Dxk.
 */
FG_eigen_matrix<double>::ptr LOL(FG_general_sparse_char_matrix::ptr input,
		FG_vector<int>::ptr labels, int ndim)
{
	typedef std::map<int, FG_vector<double>::ptr> mean_map_t;
	mean_map_t mean;
	input->group_by_mean(*labels, true, mean);
	std::vector<eigen_pair_t> eigens;
	int num_eigen = ndim - mean.size();
	if (num_eigen > 0)
		compute_SVD<FG_general_sparse_char_matrix>(input,
				2 * num_eigen, num_eigen, "LA", "RS", eigens);

	size_t nrow = input->get_num_cols();
	size_t ncol = ndim;
	FG_eigen_matrix<double>::ptr qr_matrix
		= FG_eigen_matrix<double>::create(nrow, ncol);
	qr_matrix->resize(nrow, ncol);
	printf("construct QR matrix: %ld, %ld\n", qr_matrix->get_num_rows(),
			qr_matrix->get_num_cols());
	size_t col_idx = 0;
	BOOST_FOREACH(mean_map_t::value_type v, mean) {
		if (col_idx >= ncol)
			break;
		qr_matrix->set_col(col_idx++, *v.second);
	}
	for (int i = 0; i < num_eigen; i++) {
		if (col_idx >= ncol)
			break;
		qr_matrix->set_col(col_idx++, *eigens[i].second);
	}
	return qr_matrix->householderQ();
}

FG_col_wise_matrix<double>::ptr multiply(
		FG_general_sparse_char_matrix &input1,
		FG_eigen_matrix<double> &input2)
{
	assert(input1.get_num_cols() == input2.get_num_rows());
	FG_col_wise_matrix<double>::ptr ret = FG_col_wise_matrix<double>::create(
			input1.get_num_rows(), input2.get_num_cols());
	ret->resize(input1.get_num_rows(), input2.get_num_cols());
	for (size_t i = 0; i < input2.get_num_cols(); i++) {
		FG_vector<double>::ptr vec = input2.get_col(i);
		input1.multiply(*vec, *ret->get_col_ref(i));
	}
	return ret;
}

void print_usage()
{
	fprintf(stderr, "LOL conf-file train-data train-index train-label classify_type\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-d dimensions: the number of dimensions\n");
	fprintf(stderr, "-s num: the number of samples\n");
	graph_conf.print_help();
	params.print_help();
	exit(1);
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	int ndim = 10;
	while ((opt = getopt(argc, argv, "c:d:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'd':
				ndim = atoi(optarg);
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 5) {
		print_usage();
		exit(-1);
	}

	const std::string conf_file = argv[0];
	const std::string train_data_file = argv[1];
	const std::string train_index_file = argv[2];
	const std::string train_label_file = argv[3];
	const std::string classify_type = argv[4];

	config_map configs(conf_file);
	configs.add_options(confs);

	printf("LOL starts\n");
	FG_general_sparse_char_matrix::ptr train_matrix
		= FG_general_sparse_char_matrix::create(
				FG_graph::create(train_data_file, train_index_file, configs));
	// TODO a matrix needs to detect the number of rows or the number of cols.
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

	// Create the projection matrix
	FG_eigen_matrix<double>::ptr q = LOL(train_matrix, train_labels, ndim);
	printf("q size: %ld, %ld\n", q->get_num_rows(), q->get_num_cols());
	for (size_t i = 0; i < q->get_num_cols(); i++)
		printf("||col%ld||: %f\n", i, q->get_col(i)->norm1());

	shogun::init_shogun_with_defaults();

	// Prepare training data.
	FG_col_wise_matrix<double>::ptr proj_train = multiply(*train_matrix, *q);
	printf("There are %ld rows and %ld cols in the projected train matrix\n",
			proj_train->get_num_rows(), proj_train->get_num_cols());
	shogun::SGMatrix<double> sg_train(proj_train->get_num_cols(),
			proj_train->get_num_rows());
	for (size_t i = 0; i < proj_train->get_num_cols(); i++) {
		for (size_t j = 0; j < proj_train->get_num_rows(); j++) {
			sg_train(i, j) = proj_train->get(j, i);
		}
	}
	shogun::CDenseFeatures<double>* sg_train_features
		= new shogun::CDenseFeatures<double>();
	sg_train_features->set_feature_matrix(sg_train);
	printf("There are %d feature vectors and %d features in the training data\n",
			sg_train_features->get_num_vectors(),
			sg_train_features->get_num_features());
	shogun::CMulticlassLabels* sg_train_labels=new shogun::CMulticlassLabels(
			train_labels->get_size());
	for (size_t i = 0; i < train_labels->get_size(); i++)
		sg_train_labels->set_label(i, train_labels->get(i));
	shogun::SGVector<double>::display_vector(sg_train_labels->get_labels(),
			sg_train_labels->get_num_labels());

	// Prepare testing data.
	// For now, we just test on the training data.
	FG_col_wise_matrix<double>::ptr proj_test = multiply(*train_matrix, *q);
	shogun::SGMatrix<double> sg_test(proj_test->get_num_cols(),
			proj_test->get_num_rows());
	for (size_t i = 0; i < proj_test->get_num_cols(); i++) {
		for (size_t j = 0; j < proj_test->get_num_rows(); j++) {
			sg_test(i, j) = proj_test->get(j, i);
		}
	}
	shogun::CDenseFeatures<double>* sg_test_features
		= new shogun::CDenseFeatures<double>();
	sg_test_features->set_feature_matrix(sg_test);
	printf("There are %d feature vectors and %d features in the testing data\n",
			sg_test_features->get_num_vectors(),
			sg_test_features->get_num_features());

	if (classify_type == "LDA") {
		shogun::CMCLDA* lda = new shogun::CMCLDA(sg_train_features, sg_train_labels);
		SG_REF(lda);
		lda->train();
		shogun::CMulticlassLabels* output = shogun::CLabelsFactory::to_multiclass(
				lda->apply(sg_test_features));
		SG_REF(output);
		shogun::SGVector< float64_t > output_labels = output->get_labels();
		shogun::SGVector<double>::display_vector(output_labels.vector,
				output->get_num_labels());

		size_t num_same = 0;
		for (size_t i = 0; i < train_labels->get_size(); i++)
			if (train_labels->get(i) == output_labels[i])
				num_same++;
		printf("accuracy rate: %f\n", ((double) num_same) / train_labels->get_size());
		// Free memory
		SG_UNREF(output);
		SG_UNREF(lda);
	}
	else if (classify_type == "QDA") {
		shogun::CQDA* qda = new shogun::CQDA(sg_train_features, sg_train_labels);
		SG_REF(qda);
		qda->train();
		shogun::CMulticlassLabels* output = shogun::CLabelsFactory::to_multiclass(
				qda->apply(sg_test_features));
		SG_REF(output);
		shogun::SGVector< float64_t > output_labels = output->get_labels();
		shogun::SGVector<double>::display_vector(output_labels.vector,
				output->get_num_labels());

		size_t num_same = 0;
		for (size_t i = 0; i < train_labels->get_size(); i++)
			if (train_labels->get(i) == output_labels[i])
				num_same++;
		printf("accuracy rate: %f\n", ((double) num_same) / train_labels->get_size());
		// Free memory
		SG_UNREF(output);
		SG_UNREF(qda);
	}

	shogun::exit_shogun();
}
