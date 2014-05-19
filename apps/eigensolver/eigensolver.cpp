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

const int VS_SIZE = 10;
const int RHO = 1;

typedef double ev_float_t;

/**
 * This eigensolver implements the Lanzcos algorithm.
 */

class eigen_vertex: public compute_vertex
{
	ev_float_t vs[VS_SIZE];
	ev_float_t w;
public:
	eigen_vertex() {
		memset(vs, 0, sizeof(vs));
		w = 0;
	}

	eigen_vertex(vertex_id_t id,
			const vertex_index &index): compute_vertex(id, index) {
		memset(vs, 0, sizeof(vs));
		w = 1;
//		w = random() % 1000;
	}

	void init_eigen(ev_float_t beta, int curr_v) {
		vs[curr_v] = w / beta;
	}

	ev_float_t adjust_w(ev_float_t alpha, ev_float_t beta, int curr_v) {
		w = w - alpha * vs[curr_v];
		if (curr_v > 0)
			w = w - beta * vs[curr_v - 1];
		return w;
	}

	void adjust_w2(int k, ev_float_t beta, ev_float_t sigma) {
		this->w = vs[k] * beta + w * sigma;
	}

	void set_w(ev_float_t w) {
		this->w = w;
	}

	ev_float_t get_w() const {
		return w;
	}

	void set_v(int idx, ev_float_t v) {
		this->vs[idx] = v;
	}

	ev_float_t get_v(int idx) const {
		return vs[idx];
	}

	void run(vertex_program &prog) {
		vertex_id_t id = get_id();
		request_vertices(&id, 1);
	}

	void run(vertex_program &prog, const page_vertex &vertex);

	void run_on_message(vertex_program &prog, const vertex_message &msg) {
	}
};

class eigen_vertex_program: public vertex_program_impl<eigen_vertex>
{
	ev_float_t alpha;
	int curr_v;
public:
	typedef std::shared_ptr<eigen_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<eigen_vertex_program, vertex_program>(
				prog);
	}

	eigen_vertex_program(int curr_v) {
		alpha = 0;
		this->curr_v = curr_v;
	}

	void add_vertex(const eigen_vertex &v) {
		alpha += v.get_w() * v.get_v(curr_v);
	}

	int get_curr_vidx() const {
		return curr_v;
	}

	ev_float_t get_alpha() const {
		return alpha;
	}
};

void eigen_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	eigen_vertex_program &eigen_vprog = (eigen_vertex_program &) prog;
	int curr_v = eigen_vprog.get_curr_vidx();
	edge_seq_iterator it = vertex.get_neigh_seq_it(edge_type::BOTH_EDGES, 0,
			vertex.get_num_edges(edge_type::BOTH_EDGES));
	w = 0;
	PAGE_FOREACH(vertex_id_t, id, it) {
		const eigen_vertex &eigen_v = (const eigen_vertex &) prog.get_graph().get_vertex(id);
		w += eigen_v.get_v(curr_v);
	} PAGE_FOREACH_END
	eigen_vprog.add_vertex(*this);
}

class eigen_vertex_initiator: public vertex_initiator
{
	ev_float_t beta;
	int curr_v;
public:
	eigen_vertex_initiator(ev_float_t beta, int curr_v) {
		this->beta = beta;
		this->curr_v = curr_v;
	}

	void init(compute_vertex &v) {
		eigen_vertex &ev = (eigen_vertex &) v;
		ev.init_eigen(beta, curr_v);
	}
};

class eigen_vertex_program_creater: public vertex_program_creater
{
	int curr_v;
public:
	eigen_vertex_program_creater(int curr_v) {
		this->curr_v = curr_v;
	}

	vertex_program::ptr create() const {
		return vertex_program::ptr(new eigen_vertex_program(curr_v));
	}
};

class norm2w_query: public vertex_query
{
	ev_float_t v_sq_sum;
public:
	norm2w_query() {
		v_sq_sum = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		v_sq_sum += eigen_v.get_w() * eigen_v.get_w();
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		v_sq_sum += ((norm2w_query *) q.get())->v_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new norm2w_query());
	}

	ev_float_t get_norm2() const {
		return sqrt(v_sq_sum);
	}
};

class w_query: public vertex_query
{
	const ev_float_t alpha;
	const ev_float_t beta;
	ev_float_t w_sq_sum;
	int curr_v;
public:
	w_query(ev_float_t _alpha, ev_float_t _beta, int curr_v): alpha(_alpha), beta(_beta) {
		w_sq_sum = 0;
		this->curr_v = curr_v;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		ev_float_t w = eigen_v.adjust_w(alpha, beta, curr_v);
		w_sq_sum += w * w;
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		w_sq_sum += ((w_query *) q.get())->w_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new w_query(alpha, beta, curr_v));
	}

	ev_float_t get_new_beta() {
		return sqrt(w_sq_sum);
	}
};

/**
 * This class multiplies the transpose of a matrix (n x k) with a vector of size n.
 * A row of the matrix is stored in a vertex.
 */
template<class GetLeft, class GetRight>
class matrixT_vector_multiply: public vertex_query
{
	std::vector<ev_float_t> res;
	GetLeft get_left;
	GetRight get_right;
public:
	matrixT_vector_multiply(GetLeft get_left, GetRight get_right,
			int matrix_width) {
		this->get_left = get_left;
		this->get_right = get_right;
		res.resize(matrix_width);
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		for (size_t i = 0; i < res.size(); i++) {
			res[i] += get_left(eigen_v, i) * get_right(eigen_v);
		}
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		std::vector<ev_float_t> &other_res = ((matrixT_vector_multiply *) q.get())->res;
		if (res.empty())
			res.resize(other_res.size());
		if (res.size() != other_res.size())
			printf("%ld, %ld\n", res.size(), other_res.size());
		assert(res.size() == other_res.size());
		for (size_t i = 0; i < res.size(); i++)
			res[i] += other_res[i];
	}

	virtual ptr clone() {
		return vertex_query::ptr(new matrixT_vector_multiply<GetLeft, GetRight>(
					get_left, get_right, res.size()));
	}

	size_t get_result(std::vector<ev_float_t> &res) {
		res = this->res;
		return res.size();
	}
};

/**
 * This multiplies a large matrix (n x k) with a vector of size k.
 * A row of the matrix is stored in a vertex.
 */
template<class GetLeft, class Store>
class matrix_vector_multiply: public vertex_query
{
	std::vector<ev_float_t> vec;
	GetLeft get_left;
	Store store;
public:
	matrix_vector_multiply(std::vector<ev_float_t> &vec, GetLeft get_left,
			Store store) {
		this->vec = vec;
		this->get_left = get_left;
		this->store = store;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		ev_float_t res = 0;
		for (size_t i = 0; i < vec.size(); i++)
			res += get_left(eigen_v, i) * vec[i];
		store(eigen_v, res);
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return vertex_query::ptr(new matrix_vector_multiply<GetLeft, Store>(
					vec, get_left, store));
	}
};

template<class GetLeft, class Store>
class matrix_small_matrix_multiply: public vertex_query
{
	Eigen::MatrixXd Q;
	int n_rows;
	int n_cols;
	GetLeft get_left;
	Store store;
	std::vector<ev_float_t> buf;
public:
	matrix_small_matrix_multiply(Eigen::MatrixXd &Q, int n_rows, int n_cols,
			GetLeft get_left, Store store) {
		this->Q = Q;
		this->n_rows = n_rows;
		this->n_cols = n_cols;
		this->get_left = get_left;
		this->store = store;
		buf.resize(n_cols);
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		for (int j = 0; j < n_cols; j++) {
			ev_float_t res = 0;
			for (int i = 0; i < n_rows; i++)
				res += get_left(eigen_v, i) * Q(i, j);
			buf[j] = res;
		}
		for (int j = 0; j < n_cols; j++)
			store(eigen_v, j, buf[j]);
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return vertex_query::ptr(new matrix_small_matrix_multiply<GetLeft, Store>(
					Q, n_rows, n_cols, get_left, store));
	}
};

class post_QR_adjust_w: public vertex_initiator
{
	int nv;
	ev_float_t beta;
	ev_float_t sigma;
public:
	post_QR_adjust_w(int nv, ev_float_t beta, ev_float_t sigma) {
		this->nv = nv;
		this->beta = beta;
		this->sigma = sigma;
	}

	void init(compute_vertex &v) {
		eigen_vertex &ev = (eigen_vertex &) v;
		ev.adjust_w2(nv, beta, sigma);
	}
};

class get_matrix_row
{
public:
	ev_float_t operator()(const eigen_vertex &v, int idx) {
		return v.get_v(idx);
	}
};

class get_vector_element
{
public:
	ev_float_t operator()(const eigen_vertex &v) {
		return v.get_w();
	}
};

class store_vector_element
{
public:
	void operator()(eigen_vertex &v, ev_float_t res) {
		v.set_w(v.get_w() - res);
	}
};

class store_matrix_row
{
public:
	void operator()(eigen_vertex &v, int idx, ev_float_t res) {
		v.set_v(idx, res);
	}
};

void orthogonalization(graph_engine::ptr graph, int matrix_width,
		ev_float_t &alpha, ev_float_t &beta)
{
	// transpose(V) * w
	vertex_query::ptr q(new matrixT_vector_multiply<get_matrix_row,
			get_vector_element>(get_matrix_row(), get_vector_element(), matrix_width));
	graph->query_on_all(q);
	std::vector<ev_float_t> res;
	((matrixT_vector_multiply<get_matrix_row, get_vector_element> *) q.get(
		))->get_result(res);
	assert(matrix_width >= 2);
	alpha += res[matrix_width - 1];
	beta += res[matrix_width - 2];

	// V * (transpose(V) * w)
	q = vertex_query::ptr(new matrix_vector_multiply<get_matrix_row,
			store_vector_element>(res, get_matrix_row(), store_vector_element()));
	graph->query_on_all(q);
}

/**
 * This computes norm2 of the w vector stored in the vertex state of the graph.
 */
ev_float_t norm2w(graph_engine::ptr graph)
{
	vertex_query::ptr norm2_q(new norm2w_query());
	graph->query_on_all(norm2_q);
	return ((norm2w_query *) norm2_q.get())->get_norm2();
}

void lanczos_factorization(graph_engine::ptr graph, int k, int m,
		Eigen::VectorXd &alphas, Eigen::VectorXd &betas, Eigen::MatrixXd &T)
{
	ev_float_t beta = norm2w(graph);
	if (k > 0) {
		T(k, k - 1) = beta;
		T(k - 1, k) = beta;
	}
	printf("first beta: %f\n", beta);
	for (int i = k; i < m; i++) {
		struct timeval start, end;
		struct timeval iter_start;
		gettimeofday(&start, NULL);
		iter_start = start;
		graph->start_all(vertex_initiator::ptr(new eigen_vertex_initiator(beta, i)),
				vertex_program_creater::ptr(new eigen_vertex_program_creater(i)));
		graph->wait4complete();
		gettimeofday(&end, NULL);
		printf("SPMV takes %f seconds\n", time_diff(start, end));

		// Compute alpha_j = w_j . v_j
		start = end;
		std::vector<vertex_program::ptr> vprogs;
		graph->get_vertex_programs(vprogs);
		ev_float_t alpha = 0;
		BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
			eigen_vertex_program::ptr eigen_vprog = eigen_vertex_program::cast2(vprog);
			alpha += eigen_vprog->get_alpha();
		}
		ev_float_t orth_threshold = sqrt(alpha * alpha + beta * beta) * RHO;
		gettimeofday(&end, NULL);
		printf("dot product takes %f seconds\n", time_diff(start, end));

		// Compute w_j = w_j - alpha_j * v_j - beta_j * v_j-1
		// beta_j+1 = || w_j ||
		start = end;
		vertex_query::ptr wq(new w_query(alpha, beta, i));
		graph->query_on_all(wq);
		beta = ((w_query *) wq.get())->get_new_beta();
		gettimeofday(&end, NULL);
		printf("adjusting w takes %f seconds\n", time_diff(start, end));

		if (beta < orth_threshold && i > 0) {
			start = end;
			orthogonalization(graph, i + 1, alpha, beta);
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
	fprintf(stderr, "-p: preload the graph\n");
	fprintf(stderr, "-m: the dimension of the tridiagonal matrix\n");
	fprintf(stderr, "-k: the number of required eigenvalues\n");
	graph_conf.print_help();
	params.print_help();
}

int main(int argc, char *argv[])
{
	int opt;
	std::string confs;
	int num_opts = 0;
	bool preload = false;
	int m = 0;
	int nv = 0;
	while ((opt = getopt(argc, argv, "c:pm:k:")) != -1) {
		num_opts++;
		switch (opt) {
			case 'c':
				confs = optarg;
				num_opts++;
				break;
			case 'p':
				preload = true;
				break;
			case 'm':
				m = atoi(optarg);
				num_opts++;
				break;
			case 'k':
				nv = atoi(optarg);
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

	assert(m <= VS_SIZE);
	assert(nv < m);
	graph_index::ptr index = NUMA_graph_index<eigen_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("Eigensolver starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

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
	lanczos_factorization(graph, 0, m, alphas, betas, T);
	Eigen::VectorXd I_vec;
	I_vec.conservativeResize(m);
	for (int i = 0; i < m; i++)
		I_vec(i) = 1;
	Eigen::MatrixXd I = I_vec.asDiagonal();

	while (true) {
		struct timeval start, end;
		gettimeofday(&start, NULL);
		Eigen::EigenSolver<Eigen::MatrixXd> es(T);

		Eigen::MatrixXcd eigen_vectors = es.eigenvectors();
		Eigen::VectorXcd eigen_values = es.eigenvalues();
		ev_float_t tol = 1e-8;
		int kk = 0;
		ev_float_t last_beta = betas(m - 1);
		for (int i = 0; i < m; i++) {
			if (abs(last_beta * eigen_vectors(m - 1, i).real()) < tol) {
				kk++;
				printf("eigen value[%d]: %f\n", i, eigen_values(i).real());
			}
		}
		printf("There are %d eigenvalues\n", kk);
		if (kk >= nv)
			break;

		// The p largest eigenvalues are not wanted.
		std::vector<ev_float_t> eigen_val_vec(m);
		for (int i = 0; i < m; i++)
			eigen_val_vec[i] = eigen_values(i).real();
		std::sort(eigen_val_vec.begin(), eigen_val_vec.end(),
				std::greater<ev_float_t>());

		Eigen::MatrixXd Q = I;
		for (int i = nv; i < m; i++) {
			ev_float_t mu = eigen_val_vec[i];
			Eigen::MatrixXd tmp = T - (I * mu);
			Eigen::HouseholderQR<Eigen::MatrixXd> qr = tmp.householderQr();
			Eigen::MatrixXd Qj = qr.householderQ();
			T = Qj.transpose() * T;
			T = T * Qj;
			Q = Q * Qj;
		}

		// w_k = v_k+1 + beta_k + w_m * sigma_k,
		// where beta_k = T_m[k + 1, k] and sigma_k = Q[m, k]
		ev_float_t beta_k = T(nv, nv - 1);
		std::cout << "beta: " << beta_k
			<< ", sigma: " << Q(m - 1, nv - 1) << std::endl;
		graph->init_all_vertices(vertex_initiator::ptr(new post_QR_adjust_w(
						nv, 0, Q(m - 1, nv - 1))));
		// V_k = V_m * Q[:, 1:k]
		graph->query_on_all(vertex_query::ptr(
					new matrix_small_matrix_multiply<get_matrix_row,
					store_matrix_row>(Q, m, nv, get_matrix_row(),
						store_matrix_row())));
		// T_k = T_m[1:k, 1:k]
		std::pair<int, int> keep_region_size(nv, nv);
		reset_matrix_remain(T, matrix_size, keep_region_size);
		gettimeofday(&end, NULL);
		printf("Eigen lib takes %f seconds\n", time_diff(start, end));

		lanczos_factorization(graph, nv, m, alphas, betas, T);
	}
	gettimeofday(&end, NULL);
	printf("The total running time is %f seconds\n", time_diff(start, end));

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
