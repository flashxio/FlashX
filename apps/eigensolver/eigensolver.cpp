/**
 * Copyright 2014 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
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

float beta = 0;
const int VS_SIZE = 11;

/**
 * This eigensolver implements the Lanzcos algorithm.
 */

class eigen_vertex: public compute_vertex
{
	float vs[VS_SIZE];
	float w;
	int curr_v;
public:
	eigen_vertex() {
		memset(vs, 0, sizeof(vs));
		w = 0;
		curr_v = 1;
	}

	eigen_vertex(vertex_id_t id,
			const vertex_index &index): compute_vertex(id, index) {
		memset(vs, 0, sizeof(vs));
		vs[1] = random() % 1000;
		w = 0;
		curr_v = 1;
	}

	void first_init(float normalize) {
		curr_v = 1;
		vs[1] /= normalize;
	}

	void init_eigen() {
		curr_v++;
		vs[curr_v] = w / beta;
	}

	float adjust_w(float alpha) {
		w = w - alpha * vs[curr_v] - beta * vs[curr_v - 1];
		return w;
	}

	void set_w(float w) {
		this->w = w;
	}

	float get_w() const {
		return w;
	}

	float get_v() const {
		return vs[curr_v];
	}

	float get_v(int idx) const {
		assert(idx <= curr_v);
		return vs[idx];
	}

	int get_num_vs() const {
		return curr_v + 1;
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
	float alpha;
public:
	typedef std::shared_ptr<eigen_vertex_program> ptr;

	static ptr cast2(vertex_program::ptr prog) {
		return std::static_pointer_cast<eigen_vertex_program, vertex_program>(
				prog);
	}

	eigen_vertex_program() {
		alpha = 0;
	}

	void add_vertex(const eigen_vertex &v) {
		alpha += v.get_w() * v.get_v();
	}

	float get_alpha() const {
		return alpha;
	}
};

void eigen_vertex::run(vertex_program &prog, const page_vertex &vertex)
{
	edge_seq_iterator it = vertex.get_neigh_seq_it(edge_type::BOTH_EDGES, 0,
			vertex.get_num_edges(edge_type::BOTH_EDGES));
	w = 0;
	PAGE_FOREACH(vertex_id_t, id, it) {
		const eigen_vertex &eigen_v = (const eigen_vertex &) prog.get_graph().get_vertex(id);
		w += eigen_v.get_v();
	} PAGE_FOREACH_END

	((eigen_vertex_program &) prog).add_vertex(*this);
}

class first_initiator: public vertex_initiator
{
	float normalize;
public:
	first_initiator(float normalize) {
		this->normalize = normalize;
	}

	void init(compute_vertex &v) {
		eigen_vertex &ev = (eigen_vertex &) v;
		ev.first_init(normalize);
	}
};

class eigen_vertex_initiator: public vertex_initiator
{
public:
	void init(compute_vertex &v) {
		eigen_vertex &ev = (eigen_vertex &) v;
		ev.init_eigen();
	}
};

class eigen_vertex_program_creater: public vertex_program_creater
{
public:
	vertex_program::ptr create() const {
		return vertex_program::ptr(new eigen_vertex_program());
	}
};

class norm2_query: public vertex_query
{
	float v_sq_sum;
public:
	norm2_query() {
		v_sq_sum = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		v_sq_sum += eigen_v.get_v(1) * eigen_v.get_v(1);
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		v_sq_sum += ((norm2_query *) q.get())->v_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new norm2_query());
	}

	float get_norm2() const {
		return sqrt(v_sq_sum);
	}
};

class w_query: public vertex_query
{
	const float alpha;
	float w_sq_sum;
public:
	w_query(float _alpha): alpha(_alpha) {
		w_sq_sum = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		float w = eigen_v.adjust_w(alpha);
		w_sq_sum += w * w;
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		w_sq_sum += ((w_query *) q.get())->w_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new w_query(alpha));
	}

	float get_new_beta() {
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
	std::vector<float> res;
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
		std::vector<float> &other_res = ((matrixT_vector_multiply *) q.get())->res;
		if (res.empty())
			res.resize(other_res.size());
		if (res.size() != other_res.size())
			printf("%ld, %ld\n", res.size(), other_res.size());
		assert(res.size() == other_res.size());
		for (size_t i = 0; i < res.size(); i++)
			res[i] += other_res[i];
	}

	virtual ptr clone() {
		return vertex_query::ptr(new matrixT_vector_multiply(get_left,
					get_right, res.size()));
	}

	size_t get_result(std::vector<float> &res) {
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
	std::vector<float> vec;
	float w_sq_sum;
	GetLeft get_left;
	Store store;
public:
	matrix_vector_multiply(std::vector<float> &vec, GetLeft get_left,
			Store store) {
		this->vec = vec;
		this->get_left = get_left;
		this->store = store;
		w_sq_sum = 0;
	}

	virtual void run(graph_engine &graph, compute_vertex &v) {
		eigen_vertex &eigen_v = (eigen_vertex &) v;
		assert((size_t) eigen_v.get_num_vs() == vec.size());
		float res = 0;
		for (size_t i = 0; i < vec.size(); i++)
			res += get_left(eigen_v, i) * vec[i];
		store(eigen_v, res);
		w_sq_sum += eigen_v.get_w() * eigen_v.get_w();
	}

	virtual void merge(graph_engine &graph, vertex_query::ptr q) {
		w_sq_sum += ((matrix_vector_multiply *) q.get())->w_sq_sum;
	}

	virtual ptr clone() {
		return vertex_query::ptr(new matrix_vector_multiply(vec, get_left,
					store));
	}

	float get_new_beta() {
		return sqrt(w_sq_sum);
	}
};

class get_matrix_row
{
public:
	float operator()(const eigen_vertex &v, int idx) {
		return v.get_v(idx);
	}
};

class get_vector_element
{
public:
	float operator()(const eigen_vertex &v) {
		return v.get_w();
	}
};

class store_vector_element
{
public:
	void operator()(eigen_vertex &v, float res) {
		v.set_w(v.get_w() - res);
	}
};

float orthogonalization(graph_engine::ptr graph, int matrix_width)
{
	// transpose(V) * w
	vertex_query::ptr q(new matrixT_vector_multiply<get_matrix_row,
			get_vector_element>(get_matrix_row(), get_vector_element(), matrix_width));
	graph->query_on_all(q);
	std::vector<float> res;
	((matrixT_vector_multiply<get_matrix_row, get_vector_element> *) q.get(
		))->get_result(res);

	// V * (transpose(V) * w)
	q = vertex_query::ptr(new matrix_vector_multiply<get_matrix_row,
			store_vector_element>(res, get_matrix_row(), store_vector_element()));
	graph->query_on_all(q);
	return ((matrix_vector_multiply<get_matrix_row, store_vector_element> *) q.get(
				))->get_new_beta();
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

	graph_index::ptr index = NUMA_graph_index<eigen_vertex>::create(index_file);
	graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);
	if (preload)
		graph->preload_graph();
	printf("Eigensolver starts\n");
	printf("prof_file: %s\n", graph_conf.get_prof_file().c_str());
	if (!graph_conf.get_prof_file().empty())
		ProfilerStart(graph_conf.get_prof_file().c_str());

	vertex_query::ptr norm2_q(new norm2_query());
	graph->query_on_all(norm2_q);
	float norm2_v = ((norm2_query *) norm2_q.get())->get_norm2();

	int num_curr = 0;
	int dimT = 0;
	int k = 0;
	Eigen::VectorXf alphas;
	std::vector<float> betas;
	while (num_curr < nv) {
		printf("There are %d iterations\n", m);
		dimT += m;
		alphas.conservativeResize(dimT);

		assert(VS_SIZE > dimT);
		for (int i = 0; i < m; i++) {
			struct timeval start, end;
			gettimeofday(&start, NULL);
			if (beta)
				graph->start_all(vertex_initiator::ptr(new eigen_vertex_initiator()),
						vertex_program_creater::ptr(new eigen_vertex_program_creater()));
			else
				graph->start_all(vertex_initiator::ptr(new first_initiator(norm2_v)),
						vertex_program_creater::ptr(new eigen_vertex_program_creater()));
			graph->wait4complete();

			// Compute alpha_j = w_j . v_j
			std::vector<vertex_program::ptr> vprogs;
			graph->get_vertex_programs(vprogs);
			float alpha = 0;
			BOOST_FOREACH(vertex_program::ptr vprog, vprogs) {
				eigen_vertex_program::ptr eigen_vprog = eigen_vertex_program::cast2(vprog);
				alpha += eigen_vprog->get_alpha();
			}
			alphas(k * m + i) = alpha;

			// Compute w_j = w_j - alpha_j * v_j - beta_j * v_j-1
			// beta_j+1 = || w_j ||
			vertex_query::ptr wq(new w_query(alpha));
			graph->query_on_all(wq);
			beta = ((w_query *) wq.get())->get_new_beta();

			beta = orthogonalization(graph, k * m + i + 2);
			betas.push_back(beta);
			printf("a%d: %f, b%d: %f\n", i, alpha, i + 1, beta);

			gettimeofday(&end, NULL);
			printf("Iteration %d takes %f seconds\n", i, time_diff(start, end));
		}
		k++;

		Eigen::MatrixXf T = alphas.asDiagonal();
		for (int i = 0; i < k * m - 1; i++) {
			T(i, i + 1) = betas[i];
			T(i + 1, i) = betas[i];
		}
		Eigen::EigenSolver<Eigen::MatrixXf> es(T);
		std::cout << "The eigenvalues: " << std::endl << es.eigenvalues() << std::endl;

		Eigen::MatrixXcf eigen_vectors = es.eigenvectors();
		Eigen::VectorXcf eigen_values = es.eigenvalues();
		float tol = 1e-8;
		int kk = 0;
		for (int i = 0; i < k * m; i++) {
			if (abs(beta * eigen_vectors(k * m - 1, i).real()) < tol) {
				kk++;
				printf("eigen value[%d]: %f\n", i, eigen_values(i).real());
			}
		}
		printf("There are %d eigenvalues\n", kk);
		if (kk >= nv)
			break;
	}

	if (!graph_conf.get_prof_file().empty())
		ProfilerStop();
	if (graph_conf.get_print_io_stat())
		print_io_thread_stat();
}
