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

#include <signal.h>
#include <google/profiler.h>

#include <vector>

#include <Eigen/Eigenvalues>


#include "thread.h"
#include "io_interface.h"
#include "container.h"
#include "concurrency.h"

#include "vertex_index.h"
#include "graph_engine.h"
#include "graph_config.h"
#include "FG_vector.h"
#include "matrix/FG_dense_matrix.h"
#include "matrix/FG_sparse_matrix.h"

typedef double ev_float_t;

const int RHO = 1;
const ev_float_t TOL = 1e-8;

class SPMV
{
public:
	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) = 0;
	virtual size_t get_vector_size() = 0;
};

typedef std::pair<ev_float_t, FG_vector<ev_float_t>::ptr> eigen_pair_t;

class substract_store
{
	FG_vector<ev_float_t>::ptr r;
public:
	substract_store(FG_vector<ev_float_t>::ptr r) {
		this->r = r;
	}

	void operator()(size_t idx, ev_float_t v) {
		r->set(idx, r->get(idx) - v);
	}
};

void orthogonalization(FG_col_wise_matrix<ev_float_t>::ptr V,
		FG_vector<ev_float_t>::ptr r, ev_float_t &alpha, ev_float_t &beta)
{
	// h = transpose(V) * r
	FG_vector<ev_float_t>::ptr h = V->transpose_ref()->multiply(*r);
	assert(V->get_num_cols() >= 2);
	alpha += h->get(V->get_num_cols() - 1);
	beta += h->get(V->get_num_cols() - 2);

	// r = r - V * h
	V->multiply(*h, substract_store(r));
}

class divide_apply
{
	ev_float_t beta;
public:
	divide_apply(ev_float_t beta) {
		this->beta = beta;
	}

	ev_float_t operator()(ev_float_t v) {
		return v / beta;
	}
};

class r_apply
{
	ev_float_t alpha;
	ev_float_t beta;
public:
	r_apply(ev_float_t alpha, ev_float_t beta) {
		this->alpha = alpha;
		this->beta = beta;
	}

	ev_float_t operator()(size_t idx,
			const std::vector<FG_vector<ev_float_t>::ptr> &arrs) {
		// r - alpha * v_i - beta * v_i-1
		ev_float_t ret = arrs[0]->get(idx) - alpha * arrs[1]->get(idx);
		if (arrs.size() > 2)
			ret -= beta * arrs[2]->get(idx);
		return ret;
	}
};

void lanczos_factorization(SPMV &spmv,
		FG_col_wise_matrix<ev_float_t>::ptr V, FG_vector<ev_float_t>::ptr r,
		int k, int m, Eigen::VectorXd &alphas, Eigen::VectorXd &betas,
		Eigen::MatrixXd &T)
{
	ev_float_t beta = r->norm2();
	if (k > 0) {
		T(k, k - 1) = beta;
		T(k - 1, k) = beta;
	}
	size_t num_rows = spmv.get_vector_size();
	printf("first beta: %f\n", beta);
	for (int i = k; i < m; i++) {
		struct timeval start, end;
		struct timeval iter_start;

		// v_i = r / beta
		// r = A * v_i
		gettimeofday(&start, NULL);
		iter_start = start;
		V->resize(num_rows, i + 1);
		FG_vector<ev_float_t>::ptr vi = V->get_col_ref(i);
		r->apply(divide_apply(beta), *vi);
		spmv.compute(*vi, *r);
		gettimeofday(&end, NULL);
		printf("SPMV takes %f seconds\n", time_diff(start, end));

		// Compute alpha_i = r . v_i
		start = end;
		ev_float_t alpha = vi->dot_product(*r);
		ev_float_t orth_threshold = sqrt(alpha * alpha + beta * beta) * RHO;
		gettimeofday(&end, NULL);
		printf("dot product takes %f seconds\n", time_diff(start, end));

		// Compute r = r - alpha_i * v_i - beta_i * v_i-1
		// beta_i+1 = || r ||
		start = end;
		r_apply apply(alpha, beta);
		std::vector<FG_vector<ev_float_t>::ptr> inputs;
		inputs.push_back(r);
		inputs.push_back(vi);
		if (i > 0)
			inputs.push_back(V->get_col_ref(i - 1));
		multi_vec_apply<ev_float_t, r_apply>(inputs, r, apply);
		beta = r->norm2();
		gettimeofday(&end, NULL);
		printf("adjusting w takes %f seconds\n", time_diff(start, end));

		if (beta < orth_threshold && i > 0) {
			start = end;
			assert(V->get_num_cols() == (size_t) i + 1);
			orthogonalization(V, r, alpha, beta);
			gettimeofday(&end, NULL);
			printf("orthogonalization takes %f seconds\n", time_diff(start, end));
		}

		alphas(i) = alpha;
		betas(i) = beta;
		T(i, i) = alpha;
		if (i < m - 1) {
			T(i, i + 1) = beta;
			T(i + 1, i) = beta;
		}
		printf("a%d: %f, b%d: %f\n", i, alpha, i + 1, beta);

		gettimeofday(&end, NULL);
		printf("Iteration %d takes %f seconds\n", i, time_diff(iter_start, end));
	}
}

void reset_matrix(Eigen::MatrixXd &T, std::pair<int, int> &size)
{
	for (int i = 0; i < size.first; i++)
		for (int j = 0; j < size.second; j++)
			T(i, j) = 0;
}

void reset_matrix_remain(Eigen::MatrixXd &T, std::pair<int, int> &size,
		std::pair<int, int> &keep_region_size)
{
	for (int i = keep_region_size.first; i < size.first; i++)
		for (int j = 0; j < size.second; j++)
			T(i, j) = 0;
	for (int j = keep_region_size.second; j < size.second; j++)
		for (int i = 0; i < size.first; i++)
			T(i, j) = 0;
}

// eigen values, index
typedef std::pair<ev_float_t, int> ev_pair_t;

class LA_comp
{
public:
	bool operator()(const ev_pair_t &v1, const ev_pair_t &v2) {
		return v1.first >= v2.first;
	}
};

class SA_comp
{
public:
	bool operator()(const ev_pair_t &v1, const ev_pair_t &v2) {
		return v1.first < v2.first;
	}
};

class LM_comp
{
public:
	bool operator()(const ev_pair_t &v1, const ev_pair_t &v2) {
		return abs(v1.first) >= abs(v2.first);
	}
};

class SM_comp
{
public:
	bool operator()(const ev_pair_t &v1, const ev_pair_t &v2) {
		return abs(v1.first) < abs(v2.first);
	}
};

int get_converged_eigen(Eigen::MatrixXd &T, const std::string &which,
		ev_float_t last_beta, int k, int m,
		std::vector<ev_float_t> &wanted, std::vector<ev_float_t> &unwanted,
		std::vector<FG_vector<ev_float_t>::ptr> &wanted_eigen_vectors)
{
	Eigen::EigenSolver<Eigen::MatrixXd> es(T);

	Eigen::MatrixXcd eigen_vectors = es.eigenvectors();
	Eigen::VectorXcd eigen_values = es.eigenvalues();

	std::vector<ev_pair_t> eigen_val_vec(m);
	for (int i = 0; i < m; i++) {
		eigen_val_vec[i].first = eigen_values(i).real();
		eigen_val_vec[i].second = i;
	}

	// sort the vector of eigen values so that the first k are wanted eigenvalues.
	if (which == "LA") {
		std::sort(eigen_val_vec.begin(), eigen_val_vec.end(), LA_comp());
	}
	else if (which == "SA") {
		std::sort(eigen_val_vec.begin(), eigen_val_vec.end(), SA_comp());
	}
	else if (which == "LM") {
		std::sort(eigen_val_vec.begin(), eigen_val_vec.end(), LM_comp());
	}
	else if (which == "SM") {
		std::sort(eigen_val_vec.begin(), eigen_val_vec.end(), SM_comp());
	}

	int num_converged = 0;
	for (int i = 0; i < k; i++) {
		int idx = eigen_val_vec[i].second;
		ev_float_t bound = abs(last_beta * eigen_vectors(m - 1, i).real());
		if (bound < TOL * abs(eigen_val_vec[i].first)) {
			wanted.push_back(eigen_val_vec[i].first);
			num_converged++;

			// Get the eigen vectors corresponding to the wanted eigen values.
			FG_vector<ev_float_t>::ptr eigen_vector
				= FG_vector<ev_float_t>::create(m);
			for (int j = 0; j < m; j++)
				eigen_vector->set(j, eigen_vectors(j, idx).real());
			wanted_eigen_vectors.push_back(eigen_vector);
		}
	}

	for (int i = k; i < m; i++)
		unwanted.push_back(eigen_val_vec[i].first);

	return num_converged;
}

class post_QR_apply
{
	ev_float_t beta;
	ev_float_t sigma;
public:
	post_QR_apply(ev_float_t beta, ev_float_t sigma) {
		this->beta = beta;
		this->sigma = sigma;
	}

	ev_float_t operator()(size_t idx,
			const std::vector<FG_vector<ev_float_t>::ptr> &arrs) {
		// v_k+1 * beta_k + w_m * sigma_k,
		return arrs[0]->get(idx) * beta + arrs[1]->get(idx) * sigma;
	}
};

class eigen_SPMV: public SPMV
{
	FG_adj_matrix::ptr A;
public:
	eigen_SPMV(FG_adj_matrix::ptr A) {
		this->A = A;
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->multiply(input, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

class LS_SPMV: public SPMV
{
	FG_adj_matrix::ptr A;
	FG_vector<ev_float_t>::ptr tmp;
public:
	LS_SPMV(FG_adj_matrix::ptr A) {
		this->A = A;
		tmp = FG_vector<ev_float_t>::create(A->get_num_rows());
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->transpose()->multiply(input, *tmp);
		A->multiply(*tmp, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

class RS_SPMV: public SPMV
{
	FG_adj_matrix::ptr A;
	FG_vector<ev_float_t>::ptr tmp;
public:
	RS_SPMV(FG_adj_matrix::ptr A) {
		this->A = A;
		tmp = FG_vector<ev_float_t>::create(A->get_num_rows());
	}

	virtual void compute(const FG_vector<ev_float_t> &input,
			FG_vector<ev_float_t> &output) {
		A->multiply(input, *tmp);
		A->transpose()->multiply(*tmp, output);
	}

	virtual size_t get_vector_size() {
		return A->get_num_rows();
	}
};

void eigen_solver(SPMV &spmv, int m, int nv, const std::string &which,
		std::vector<eigen_pair_t> &eigen_pairs)
{
	FG_col_wise_matrix<ev_float_t>::ptr V
		= FG_col_wise_matrix<ev_float_t>::create(spmv.get_vector_size(), m);
	FG_vector<ev_float_t>::ptr r = FG_vector<ev_float_t>::create(
			spmv.get_vector_size());
	r->init(1);

	struct timeval start, end;
	gettimeofday(&start, NULL);
	Eigen::MatrixXd T;
	Eigen::VectorXd betas;
	Eigen::VectorXd alphas;
	T.conservativeResize(m, m);
	std::pair<int, int> matrix_size(m, m);
	reset_matrix(T, matrix_size);
	betas.conservativeResize(m);
	alphas.conservativeResize(m);

	Eigen::VectorXd I_vec;
	I_vec.conservativeResize(m);
	for (int i = 0; i < m; i++)
		I_vec(i) = 1;
	Eigen::MatrixXd I = I_vec.asDiagonal();

	std::vector<FG_vector<ev_float_t>::ptr> wanted_eigen_vectors;
	std::vector<ev_float_t> wanted_eigen_values;
	lanczos_factorization(spmv, V, r, 0, m, alphas, betas, T);
	while (true) {
		struct timeval start, end;
		gettimeofday(&start, NULL);

		std::vector<ev_float_t> unwanted;
		wanted_eigen_vectors.clear();
		wanted_eigen_values.clear();
		int num_converged = get_converged_eigen(T, which, betas(m - 1),
				nv, m, wanted_eigen_values, unwanted, wanted_eigen_vectors);
		if (num_converged >= nv)
			break;

		Eigen::MatrixXd Q = I;
		assert(unwanted.size() == (size_t) (m - nv));
		for (int i = 0; i < m - nv; i++) {
			ev_float_t mu = unwanted[i];
			Eigen::MatrixXd tmp = T - (I * mu);
			Eigen::HouseholderQR<Eigen::MatrixXd> qr = tmp.householderQr();
			Eigen::MatrixXd Qj = qr.householderQ();
			T = Qj.transpose() * T;
			T = T * Qj;
			Q = Q * Qj;
		}

		// w_k = v_k+1 * beta_k + w_m * sigma_k,
		// where beta_k = T_m[k + 1, k] and sigma_k = Q[m, k]
		ev_float_t beta_k = T(nv, nv - 1);
		std::cout << "beta: " << beta_k
			<< ", sigma: " << Q(m - 1, nv - 1) << std::endl;
		std::vector<FG_vector<ev_float_t>::ptr> inputs(2);
		inputs[0] = V->get_col_ref(nv);
		inputs[1] = r;
		post_QR_apply apply(beta_k, Q(m - 1, nv - 1));
		multi_vec_apply<ev_float_t, post_QR_apply>(inputs, r, apply);
		// V_k = V_m * Q[:, 1:k]
		FG_eigen_matrix<ev_float_t> subQ(Q, m, nv);
		printf("subQ: %ld, %ld\n", subQ.get_num_rows(), subQ.get_num_cols());
		printf("V: %ld, %ld\n", V->get_num_rows(), V->get_num_cols());
		V->multiply_in_place(subQ);
		// T_k = T_m[1:k, 1:k]
		std::pair<int, int> keep_region_size(nv, nv);
		reset_matrix_remain(T, matrix_size, keep_region_size);
		gettimeofday(&end, NULL);
		printf("Eigen lib takes %f seconds\n", time_diff(start, end));

		lanczos_factorization(spmv, V, r, nv, m, alphas, betas, T);
	}
	gettimeofday(&end, NULL);
	printf("The total running time is %f seconds\n", time_diff(start, end));

	assert((size_t) nv == wanted_eigen_vectors.size());
	std::vector<FG_vector<ev_float_t>::ptr> orig_eigen_vectors(nv);
	for (int i = 0; i < nv; i++)
		orig_eigen_vectors[i] = V->multiply(*wanted_eigen_vectors[i]);

	assert(wanted_eigen_values.size() == orig_eigen_vectors.size());
	for (size_t i = 0; i < wanted_eigen_values.size(); i++) {
		eigen_pairs.push_back(eigen_pair_t(wanted_eigen_values[i],
					orig_eigen_vectors[i]));
	}
}

void int_handler(int sig_num)
{
	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	exit(0);
}

void print_usage()
{
	fprintf(stderr,
			"eigensolver [options] conf_file graph_file index_file\n");
	fprintf(stderr, "-c confs: add more configurations to the system\n");
	fprintf(stderr, "-m: the dimension of the tridiagonal matrix\n");
	fprintf(stderr, "-k: the number of required eigenvalues\n");
	fprintf(stderr, "-w: which type of eigenvalues\n");
	fprintf(stderr, "-t: the type of eigenvalues and eigenvectors.\n");
	fprintf(stderr, "\tEV (eigen vectors of a symmetric matrix),\n");
	fprintf(stderr, "\tLS (left-singular vectors of SVD),\n");
	fprintf(stderr, "\tRS (right-singular vectors of SVD)\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	int m = 0;
	int nv = 0;
	std::string which = "LA";
	std::string type = "EV";
	while ((opt = getopt(argc, argv, "c:m:k:w:t:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'm':
				m = atoi(optarg);
				num_opts++;
				break;
			case 'k':
				nv = atoi(optarg);
				num_opts++;
				break;
			case 'w':
				which = optarg;
				num_opts++;
				break;
			case 't':
				type = optarg;
				num_opts++;
				break;
			default:
				print_usage();
		}
	}
	argv += 1 + num_opts;
	argc -= 1 + num_opts;

	if (argc < 3) {
		print_usage();
		exit(-1);
	}

	std::string conf_file = argv[0];
	std::string graph_file = argv[1];
	std::string index_file = argv[2];

	config_map configs(conf_file);
	configs.add_options(confs);

	signal(SIGINT, int_handler);

	assert(nv < m);
	printf("Eigensolver starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	FG_graph::ptr graph = FG_graph::create(graph_file, index_file, configs);

	std::vector<eigen_pair_t> eigen_pairs;
	if (type == "EV") {
		eigen_SPMV spmv(FG_adj_matrix::create(graph));
		eigen_solver(spmv, m, nv, which, eigen_pairs);
	}
	else if (type == "LS") {
		LS_SPMV spmv(FG_adj_matrix::create(graph));
		eigen_solver(spmv, m, nv, which, eigen_pairs);
	}
	else if (type == "RS") {
		RS_SPMV spmv(FG_adj_matrix::create(graph));
		eigen_solver(spmv, m, nv, which, eigen_pairs);
	}
	else
		assert(0);
	assert(eigen_pairs.size() == (size_t) nv);

	for (int i = 0; i < nv; i++) {
		ev_float_t v = eigen_pairs[i].first;
		if (type == "LS" || type == "RS")
			v = sqrt(v);
		printf("eigen value: %f, norm1(vector): %f\n",
				v, eigen_pairs[i].second->norm1());
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
