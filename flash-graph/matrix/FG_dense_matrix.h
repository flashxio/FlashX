#ifndef __FG_DENSE_MATRIX_H__
#define __FG_DENSE_MATRIX_H__

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

#include <vector>
#include <memory>

#include <Eigen/Dense>

#include "FG_vector.h"
#include "graph.h"

template<class T>
class col_wise_matrix_store
{
	std::vector<typename FG_vector<T>::ptr> cols;
public:
	col_wise_matrix_store(size_t nrow, size_t ncol) {
		cols.resize(ncol);
		for (size_t i = 0; i < ncol; i++)
			cols[i] = FG_vector<T>::create(nrow);
	}

	void set(size_t row, size_t col, const T &v) {
		return cols[col]->set(row, v);
	}

	const T &get(size_t row, size_t col) const {
		return cols[col]->get(row);
	}

	typename FG_vector<T>::ptr get(size_t col) {
		return cols[col];
	}

	const std::vector<typename FG_vector<T>::ptr> &get_cols() const {
		return cols;
	}

	void set_cols(const std::vector<typename FG_vector<T>::ptr> &cols) {
		this->cols = cols;
	}

	size_t get_num_rows() const {
		return cols[0]->get_size();
	}

	size_t get_num_cols() const {
		return cols.size();
	}
};

template<class T>
class row_wise_matrix_store
{
	std::vector<typename FG_vector<T>::ptr> rows;
public:
	row_wise_matrix_store(size_t nrow, size_t ncol) {
		rows.resize(nrow);
		for (size_t i = 0; i < nrow; i++)
			rows[i] = FG_vector<T>::create(ncol);
	}

	void set(size_t row, size_t col, const T &v) {
		return rows[row]->set(col, v);
	}

	const T &get(size_t row, size_t col) const {
		return rows[row]->get(col);
	}

	typename FG_vector<T>::ptr get(size_t row) {
		return rows[row];
	}

	const std::vector<typename FG_vector<T>::ptr> &get_rows() const {
		return rows;
	}

	void set_rows(const std::vector<typename FG_vector<T>::ptr> &rows) {
		this->rows = rows;
	}

	size_t get_num_rows() const {
		return rows.size();
	}

	size_t get_num_cols() const {
		return rows[0]->get_size();
	}
};

template<class T>
class Eigen_matrix_store
{
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;
public:
	typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> matrix_type;

	Eigen_matrix_store(size_t nrow, size_t ncol) {
		mat.resize(nrow, ncol);
		assert((size_t) mat.rows() == nrow);
		assert((size_t) mat.cols() == ncol);
	}

	Eigen_matrix_store(
			const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat) {
		this->mat = mat;
	}

	void set(size_t row, size_t col, const T &v) {
		mat(row, col) = v;
	}

	const T &get(size_t row, size_t col) const {
		return mat(row, col);
	}

	size_t get_num_rows() const {
		return mat.rows();
	}

	size_t get_num_cols() const {
		return mat.cols();
	}

	const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &get_matrix() const {
		return mat;
	}

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &get_matrix() {
		return mat;
	}
};

template<class T, class MatrixStore>
class FG_dense_matrix
{
protected:
	// The number of rows and columns used by the matrix.
	size_t nrow;
	size_t ncol;
	// The data structure storing the matrix. Its space needs to be allocated
	// in advance.
	MatrixStore matrix_store;

	FG_dense_matrix(): matrix_store(0, 0) {
		nrow = 0;
		ncol = 0;
	}

	FG_dense_matrix(size_t nrow, size_t ncol): matrix_store(nrow, ncol) {
		this->nrow = 0;
		this->ncol = 0;
	}
public:
	typedef std::shared_ptr<FG_dense_matrix<T, MatrixStore> > ptr;

	/**
	 * Create a dense matrix.
	 * `nrow' and `ncol' specify the reserved space for the matrix.
	 */
	static ptr create(size_t nrow, size_t ncol) {
		return ptr(new FG_dense_matrix<T, MatrixStore>(nrow, ncol));
	}

	/**
	 * Resize the matrix.
	 * `nrow' and `ncol' defined the size of the matrix. They must be smaller
	 * than or equal to the space reserved for the matrix.
	 */
	void resize(size_t nrow, size_t ncol) {
		assert(matrix_store.get_num_rows() >= nrow);
		assert(matrix_store.get_num_cols() >= ncol);
		this->nrow = nrow;
		this->ncol = ncol;
	}

	/**
	 * Multiply the matrix by a vector and return the result vector.
	 */
	typename FG_vector<T>::ptr multiply(FG_vector<T> &vec) const {
		typename FG_vector<T>::ptr ret = FG_vector<T>::create(nrow);

		struct {
			typename FG_vector<T>::ptr fg_vec;
			void operator()(size_t i, const T &v) {
				fg_vec->set(i, v);
			}
		} identity_store;
		identity_store.fg_vec = ret;
		multiply(vec, identity_store);
		return ret;
	}

	/**
	 * Multiply the matrix by a vector and the user can specify how the result
	 * is stored.
	 */
	template<class Store>
	void multiply(FG_vector<T> &vec, Store store) const {
		assert(vec.get_size() == ncol);
#pragma omp parallel for
		for (size_t i = 0; i < nrow; i++) {
			T res = 0;
			for (size_t j = 0; j < ncol; j++)
				res += get(i, j) * vec.get(j);
			store(i, res);
		}
	}

	/**
	 * Multiply the matrix by another matrix.
	 * The other matrix needs to have fewer columns so that the result can be
	 * stored in the same matrix. The multiplication also changes the number
	 * of columns of the current matrix.
	 */
	template<class MatrixStore1>
	void multiply_in_place(FG_dense_matrix<T, MatrixStore1> &matrix) {
		assert(ncol == matrix.get_num_rows());
		assert(ncol >= matrix.get_num_cols());
		std::vector<T> buf(matrix.get_num_cols());
#pragma omp parallel for private(buf)
		for (size_t i = 0; i < nrow; i++) {
			buf.resize(matrix.get_num_cols());
			for (size_t j = 0; j < matrix.get_num_cols(); j++) {
				T res = 0;
				for (size_t k = 0; k < ncol; k++)
					res += get(i, k) * matrix.get(k, j);
				buf[j] = res;
			}
			for (size_t j = 0; j < matrix.get_num_cols(); j++)
				matrix_store.set(i, j, buf[j]);
		}
		ncol = matrix.get_num_cols();
	}

	const T &get(size_t row, size_t col) const {
		return matrix_store.get(row, col);
	}

	size_t get_num_rows() const {
		return nrow;
	}

	size_t get_num_cols() const {
		return ncol;
	}
};

template<class T> class FG_col_wise_matrix;
template<class T> class FG_row_wise_matrix;

template<class T>
class FG_row_wise_matrix: public FG_dense_matrix<T, row_wise_matrix_store<T> >
{
	FG_row_wise_matrix(size_t nrow,
			size_t ncol): FG_dense_matrix<T, row_wise_matrix_store<T> >(
				nrow, ncol) {
		// We assume the row-wise matrix has more columns than rows.
		assert(nrow < ncol);
	}

	FG_row_wise_matrix(const FG_col_wise_matrix<T> &mat, bool transpose) {
		if (transpose) {
			this->nrow = mat.get_num_cols();
			this->ncol = mat.get_num_rows();
			this->matrix_store.set_rows(mat.matrix_store.get_cols());
		}
		else {
			// TODO
			assert(0);
		}
		// We assume the row-wise matrix has more columns than rows.
		assert(this->nrow < this->ncol);
	}
public:
	typedef std::shared_ptr<FG_row_wise_matrix<T> > ptr;

	static ptr create(size_t nrow, size_t ncol) {
		return ptr(new FG_row_wise_matrix<T>(nrow, ncol));
	}

	typename FG_vector<T>::ptr get_row_ref(size_t row) {
		assert(row < this->get_num_rows());
		return this->matrix_store.get(row);
	}

	graph::ptr conv2graph() const {
		// TODO we can generate an undirected graph for a symmetric matrix.
		typename directed_graph<T>::ptr g = directed_graph<T>::create(true);
		for (size_t i = 0; i < this->get_num_rows(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_cols(); j++) {
				edge<T> e(i, j, this->get(i, j));
				v.add_out_edge(e);
			}
			for (size_t j = 0; j < this->get_num_rows(); j++) {
				edge<T> e(j, i, this->get(j, i));
				v.add_in_edge(e);
			}
			g->add_vertex(v);
		}
		// There are more columns than rows in the matrix.
		for (size_t i = this->get_num_rows(); i < this->get_num_cols(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_rows(); j++) {
				edge<T> e(j, i, this->get(j, i));
				v.add_in_edge(e);
			}
			g->add_vertex(v);
		}

		return g;
	}

	template<class T1>
	friend class FG_col_wise_matrix;
};

template<class T>
class FG_col_wise_matrix: public FG_dense_matrix<T, col_wise_matrix_store<T> >
{
	FG_col_wise_matrix(size_t nrow,
			size_t ncol): FG_dense_matrix<T, col_wise_matrix_store<T> >(
				nrow, ncol) {
	}
public:
	typedef std::shared_ptr<FG_col_wise_matrix<T> > ptr;

	static ptr create(size_t nrow, size_t ncol) {
		return ptr(new FG_col_wise_matrix<T>(nrow, ncol));
	}

	typename FG_vector<T>::ptr get_col_ref(size_t col) {
		assert(col < this->get_num_cols());
		return this->matrix_store.get(col);
	}

	typename FG_row_wise_matrix<T>::ptr transpose_ref() {
		return typename FG_row_wise_matrix<T>::ptr(
				new FG_row_wise_matrix<T>(*this, true));
	}

	template<class T1>
	friend class FG_row_wise_matrix;
};

template<class T>
class FG_eigen_matrix: public FG_dense_matrix<T, Eigen_matrix_store<T> >
{
	FG_eigen_matrix(size_t nrow,
			size_t ncol): FG_dense_matrix<T, Eigen_matrix_store<T> >(nrow, ncol) {
	}
public:
	typedef std::shared_ptr<FG_eigen_matrix<T> > ptr;

	FG_eigen_matrix(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat,
			size_t nrow, size_t ncol): FG_dense_matrix<T, Eigen_matrix_store<T> >(
				nrow, ncol) {
		this->matrix_store = Eigen_matrix_store<T>(mat);
		assert(nrow <= this->matrix_store.get_num_rows());
		assert(ncol <= this->matrix_store.get_num_cols());
		this->nrow = nrow;
		this->ncol = ncol;
	}

	static ptr create(size_t nrow, size_t ncol) {
		return ptr(new FG_eigen_matrix<T>(nrow, ncol));
	}

	void set_col(size_t idx, const FG_vector<T> &vec) {
		assert(vec.get_size() == this->get_num_rows());
		for (size_t i = 0; i < vec.get_size(); i++)
			this->matrix_store.set(i, idx, vec.get(i));
	}

	void set_row(size_t idx, const FG_vector<T> &vec) {
		assert(vec.get_size() == this->get_num_cols());
		for (size_t i = 0; i < vec.get_size(); i++)
			this->matrix_store.set(idx, i, vec.get(i));
	}

	typename FG_eigen_matrix<T>::ptr householderQ() const {
		size_t nrows = this->matrix_store.get_num_rows();
		size_t ncols = this->matrix_store.get_num_cols();
		typename FG_eigen_matrix<T>::ptr ret = FG_eigen_matrix<T>::ptr(
				new FG_eigen_matrix(nrows, ncols));
		ret->matrix_store.get_matrix()
			= this->matrix_store.get_matrix().householderQr().householderQ(
					) * Eigen::MatrixXd::Identity(nrows, ncols);
		ret->nrow = ret->matrix_store.get_num_rows();
		ret->ncol = ret->matrix_store.get_num_cols();
		return ret;
	}

	typename FG_vector<T>::ptr get_col(size_t col) const {
		assert(col < this->ncol);
		typename FG_vector<T>::ptr vec
			= FG_vector<T>::create(this->get_num_rows());
		for (size_t i = 0; i < this->get_num_rows(); i++)
			vec->set(i, this->get(i, col));
		return vec;
	}

	typename FG_vector<T>::ptr get_row(size_t row) const {
		assert(row < this->nrow);
		typename FG_vector<T>::ptr vec
			= FG_vector<T>::create(this->get_num_cols());
		for (size_t i = 0; i < this->get_num_cols(); i++)
			vec->set(i, this->get(row, i));
		return vec;
	}

	FG_eigen_matrix<T>::ptr transpose() const {
		return FG_eigen_matrix<T>::ptr(new FG_eigen_matrix<T>(
					this->matrix_store.get_matrix().transpose(),
					this->get_num_cols(), this->get_num_rows()));
	}

	graph::ptr conv2graph_wideM() const {
		// TODO we can generate an undirected graph for a symmetric matrix.
		typename directed_graph<T>::ptr g = directed_graph<T>::create(true);
		for (size_t i = 0; i < this->get_num_rows(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_cols(); j++) {
				edge<T> e(i, j, this->get(i, j));
				v.add_out_edge(e);
			}
			for (size_t j = 0; j < this->get_num_rows(); j++) {
				edge<T> e(j, i, this->get(j, i));
				v.add_in_edge(e);
			}
			g->add_vertex(v);
		}
		// There are more columns than rows in the matrix.
		for (size_t i = this->get_num_rows(); i < this->get_num_cols(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_rows(); j++) {
				edge<T> e(j, i, this->get(j, i));
				v.add_in_edge(e);
			}
			g->add_vertex(v);
		}

		return g;
	}

	graph::ptr conv2graph_tallM() const {
		// TODO we can generate an undirected graph for a symmetric matrix.
		typename directed_graph<T>::ptr g = directed_graph<T>::create(true);
		for (size_t i = 0; i < this->get_num_cols(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_cols(); j++) {
				edge<T> e(i, j, this->get(i, j));
				v.add_out_edge(e);
			}
			for (size_t j = 0; j < this->get_num_rows(); j++) {
				edge<T> e(j, i, this->get(j, i));
				v.add_in_edge(e);
			}
			g->add_vertex(v);
		}
		// There are more rows than cols in the matrix.
		for (size_t i = this->get_num_cols(); i < this->get_num_rows(); i++) {
			in_mem_directed_vertex<T> v(i, true);
			for (size_t j = 0; j < this->get_num_cols(); j++) {
				edge<T> e(i, j, this->get(i, j));
				v.add_out_edge(e);
			}
			g->add_vertex(v);
		}

		return g;
	}

	graph::ptr conv2graph() const {
		if (this->get_num_rows() > this->get_num_cols())
			return conv2graph_tallM();
		else
			return conv2graph_wideM();
	}
};

#endif
