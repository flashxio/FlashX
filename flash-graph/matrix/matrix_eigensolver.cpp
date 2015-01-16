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

#ifdef USE_EIGEN

#include <vector>

#include <eigen3/Eigen/Eigenvalues>

#include "graph_engine.h"
#include "graph_config.h"
#include "FG_vector.h"
#include "FG_dense_matrix.h"
#include "FG_sparse_matrix.h"
#include "FGlib.h"
#include "matrix_eigensolver.h"

namespace fg
{

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
		FG_vector<ev_float_t>::ptr r, ev_float_t &alpha)
{
	// h = transpose(V) * r
	FG_vector<ev_float_t>::ptr h = V->transpose_ref()->multiply(*r);
	assert(V->get_num_cols() >= 2);
	alpha += h->get(V->get_num_cols() - 1);

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

ev_float_t lanczos_factorization(SPMV &spmv, int k, int m,
		// Input/output
		FG_col_wise_matrix<ev_float_t>::ptr V, FG_vector<ev_float_t>::ptr r,
		// output
		Eigen::MatrixXd &T)
{
	ev_float_t beta = r->norm2();
	size_t num_rows = spmv.get_vector_size();
	BOOST_LOG_TRIVIAL(info) << "first beta: " << beta;
	for (int i = k; i < m; i++) {
		struct timeval start, end;
		struct timeval iter_start;
		bool restart = false;

		gettimeofday(&iter_start, NULL);
		// Step 1: we need to reinitialize r if necessary
		if (beta == 0) {
			// TODO
			assert(0);
		}

		// Step 2
		// v_i = r / beta
		gettimeofday(&start, NULL);
		V->resize(num_rows, i + 1);
		FG_vector<ev_float_t>::ptr vi = V->get_col_ref(i);
		r->apply(divide_apply(beta), *vi);

		// Step 3
		// r = A * v_i
		spmv.compute(*vi, *r);
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("SPMV takes %1% seconds")
			% time_diff(start, end);

		// Step 4
		ev_float_t wnorm = r->norm2();
		start = end;
		FG_vector<ev_float_t>::ptr w = V->transpose_ref()->multiply(*r);
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("t(V)*r takes %1% seconds, V has %2% cols")
			% time_diff(start, end) % V->get_num_cols();
		// r = r - V * h
		start = end;
		V->multiply(*w, substract_store(r));
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("V*w takes %1% seconds")
			% time_diff(start, end);
		assert((size_t) i == w->get_size() - 1);
		ev_float_t alpha = w->get(i);
		T(i, i) = alpha;
		if (i > 0) {
			if (restart) {
				T(i, i - 1) = 0;
				T(i - 1, i) = 0;
			}
			else {
				T(i, i - 1) = beta;
				T(i - 1, i) = beta;
			}
		}
		start = end;
		beta = r->norm2();
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("norm2 takes %1% seconds")
			% time_diff(start, end);

		// Step 5
		int orth_iter = 0;
		start = end;
		while (beta <= 0.717 * wnorm && orth_iter < 2) {
			orth_iter++;
			assert(V->get_num_cols() == (size_t) i + 1);
			orthogonalization(V, r, alpha);
			ev_float_t new_beta = r->norm2();
			if (new_beta > 0.717 * beta) {
				beta = new_beta;
				break;
			}
			else if (orth_iter < 2) {
				beta = new_beta;
				continue;
			}
			else {
				r->init(0);
				beta = 0;
				break;
			}
		}
		BOOST_LOG_TRIVIAL(info) << boost::format("a%1%: %2%, b%3%: %4%")
			% i % alpha % (i + 1) % beta;
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("reorthogonalization takes %1% seconds")
			% time_diff(start, end);

		if (i > 0 && T(i, i - 1) < 0) {
			T(i, i - 1) = -T(i, i - 1);
			T(i - 1, i) = -T(i - 1, i);
			if (i < m - 1)
				V->get_col_ref(i)->neg_in_place();
			else
				r->neg_in_place();
		}

		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("Iteration %1% takes %2% seconds")
			% i % time_diff(iter_start, end);
	}
	return beta;
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
		return std::abs(v1.first) >= std::abs(v2.first);
	}
};

class SM_comp
{
public:
	bool operator()(const ev_pair_t &v1, const ev_pair_t &v2) {
		return std::abs(v1.first) < std::abs(v2.first);
	}
};

int get_converged_eigen(Eigen::MatrixXd &T, const std::string &which,
		ev_float_t last_beta, int k, int m, double tol,
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
		ev_float_t bound = std::abs(last_beta * eigen_vectors(m - 1, idx).real());
		if (bound < tol * std::abs(eigen_val_vec[i].first)) {
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

void eigen_solver(SPMV &spmv, int m, int nv, const std::string &which,
		double tol, std::vector<eigen_pair_t> &eigen_pairs)
{
	FG_col_wise_matrix<ev_float_t>::ptr V
		= FG_col_wise_matrix<ev_float_t>::create(spmv.get_vector_size(), m);
	FG_vector<ev_float_t>::ptr r = FG_vector<ev_float_t>::create(
			spmv.get_vector_size());
	r->init_rand();

	struct timeval start, end;
	gettimeofday(&start, NULL);
	Eigen::MatrixXd T;
	T.conservativeResize(m, m);
	std::pair<int, int> matrix_size(m, m);
	reset_matrix(T, matrix_size);

	Eigen::VectorXd I_vec;
	I_vec.conservativeResize(m);
	for (int i = 0; i < m; i++)
		I_vec(i) = 1;
	Eigen::MatrixXd I = I_vec.asDiagonal();

	std::vector<FG_vector<ev_float_t>::ptr> wanted_eigen_vectors;
	std::vector<ev_float_t> wanted_eigen_values;
	ev_float_t beta = lanczos_factorization(spmv, 0, m, V, r, T);
	while (true) {
		struct timeval start, end;
		gettimeofday(&start, NULL);

		std::vector<ev_float_t> unwanted;
		wanted_eigen_vectors.clear();
		wanted_eigen_values.clear();
		int num_converged = get_converged_eigen(T, which, beta,
				nv, m, tol, wanted_eigen_values, unwanted, wanted_eigen_vectors);
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
		BOOST_LOG_TRIVIAL(info) << boost::format("beta: %1%, sigma: %2%")
			% beta_k % Q(m - 1, nv - 1);
#if 0
		// TODO we should fix this.
		std::vector<FG_vector<ev_float_t>::ptr> inputs(2);
		inputs[0] = V->get_col_ref(nv);
		inputs[1] = r;
		post_QR_apply apply(beta_k, Q(m - 1, nv - 1));
		multi_vec_apply<ev_float_t, post_QR_apply>(inputs, r, apply);
#endif
		r->multiply_in_place(Q(m - 1, nv - 1));
		// V_k = V_m * Q[:, 1:k]
		FG_eigen_matrix<ev_float_t> subQ(Q, m, nv);
		BOOST_LOG_TRIVIAL(info) << boost::format("subQ: %1%, %2%")
			% subQ.get_num_rows() % subQ.get_num_cols();
		BOOST_LOG_TRIVIAL(info) << boost::format("V: %1%, %2%")
			% V->get_num_rows() % V->get_num_cols();
		V->multiply_in_place(subQ);
		// T_k = T_m[1:k, 1:k]
		std::pair<int, int> keep_region_size(nv, nv);
		reset_matrix_remain(T, matrix_size, keep_region_size);
		gettimeofday(&end, NULL);
		BOOST_LOG_TRIVIAL(info) << boost::format("Eigen lib takes %1% seconds")
			% time_diff(start, end);

		lanczos_factorization(spmv, nv, m, V, r, T);
	}
	gettimeofday(&end, NULL);
	BOOST_LOG_TRIVIAL(info) << boost::format("The total running time is %1% seconds")
		% time_diff(start, end);

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

}

#endif
