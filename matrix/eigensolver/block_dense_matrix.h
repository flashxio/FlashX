#ifndef __BLOCK_DENSE_MATRIX_H__
#define __BLOCK_DENSE_MATRIX_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include "mem_dense_matrix.h"
#include "generic_type.h"
#include "mem_vector.h"
#include "sparse_matrix.h"

class sp_multiply
{
public:
	virtual void run(const fm::sparse_matrix &A, const fm::mem_vector &in,
			fm::mem_vector &out) const = 0;
	virtual void run(const fm::sparse_matrix &A, const fm::NUMA_vector &in,
			fm::NUMA_vector &out) const = 0;
	virtual void run(const fm::sparse_matrix &A, const fm::detail::matrix_store &in,
			fm::detail::matrix_store &out) const = 0;
};

template<class T>
class sp_multiply_impl: public sp_multiply
{
public:
	virtual void run(const fm::sparse_matrix &A, const fm::mem_vector &in,
			fm::mem_vector &out) const {
		A.multiply<T>(in, out);
	}
	virtual void run(const fm::sparse_matrix &A, const fm::NUMA_vector &in,
			fm::NUMA_vector &out) const {
		A.multiply<T>(in, out);
	}
	virtual void run(const fm::sparse_matrix &A, const fm::detail::matrix_store &in,
			fm::detail::matrix_store &out) const {
		A.multiply<T>(in, out);
	}
};

class block_multi_vector
{
	size_t block_size;
	size_t num_rows;
	size_t num_cols;
	const fm::scalar_type &type;

	static void sparse_matrix_multiply(const sp_multiply &multiply,
			const fm::sparse_matrix &A, const block_multi_vector &X,
			block_multi_vector &Y);
protected:
	std::vector<fm::dense_matrix::ptr> mats;

	block_multi_vector(size_t nrow, size_t ncol, size_t block_size,
			const fm::scalar_type &_type);
	block_multi_vector(const std::vector<fm::dense_matrix::ptr> &mats);
public:
	typedef std::shared_ptr<block_multi_vector> ptr;

	static ptr create(size_t nrow, size_t ncol, size_t block_size,
			const fm::scalar_type &type) {
		assert(ncol % block_size == 0);
		return ptr(new block_multi_vector(nrow, ncol, block_size, type));
	}

	bool resize_block(size_t new_block_size);

	size_t get_num_rows() const {
		return num_rows;
	}

	size_t get_num_cols() const {
		return num_cols;
	}

	size_t get_num_blocks() const {
		return num_cols / block_size;
	}

	size_t get_entry_size() const {
		return type.get_size();
	}

	size_t get_block_size() const {
		return block_size;
	}

	int get_num_nodes() const {
		return fm::matrix_conf.get_num_nodes();
	}

	const fm::scalar_type &get_type() const {
		return type;
	}

	fm::dense_matrix::const_ptr get_block(off_t block_idx) const {
		return mats[block_idx];
	}

	fm::dense_matrix::ptr get_block(off_t block_idx) {
		return mats[block_idx];
	}

	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat) {
		assert(mat->get_num_cols() == block_size);
		mats[block_idx] = mat->clone();
	}

	void set_block(const block_multi_vector &mv, const std::vector<int>& index);

	const char *get_col_raw(off_t col_idx) const;

	fm::dense_matrix::const_ptr get_col(off_t col_idx) const;

	block_multi_vector::ptr get_cols(const std::vector<int> &index);
	/*
	 * The difference between get_cols and get_cols_mirror is that any change
	 * to the multi_vector returned by get_cols() doesn't affect the original
	 * multi_vector, but any changes to the multi_vector returned by
	 * get_cols_mirror() affect the original multi_vector.
	 */
	block_multi_vector::ptr get_cols_mirror(const std::vector<int> &index);

	block_multi_vector::ptr clone() const;

	block_multi_vector::ptr gemm(const block_multi_vector &A,
			const fm::detail::mem_col_matrix_store &B,
			const fm::scalar_variable &alpha,
			const fm::scalar_variable &beta) const;

	void assign(const block_multi_vector &vecs);

	block_multi_vector::ptr add(const block_multi_vector &vecs) const;
	fm::mem_dense_matrix::ptr MvTransMv(const block_multi_vector &mv) const;

	template<class Type>
	void init_rand(Type min, Type max) {
		size_t num_blocks = get_num_blocks();
		for (size_t i = 0; i < num_blocks; i++)
			set_block(i, fm::mem_dense_matrix::create_rand<Type>(min, max,
						get_num_rows(), block_size, fm::matrix_layout_t::L_COL,
						get_num_nodes()));
	}

	template<class Type>
	block_multi_vector::ptr multiply_scalar(Type val) const {
		size_t num_blocks = get_num_blocks();
		block_multi_vector::ptr ret_vecs = block_multi_vector::create(
				get_num_rows(), get_num_cols(), block_size, type);
		for (size_t i = 0; i < num_blocks; i++)
			ret_vecs->set_block(i, get_block(i)->multiply_scalar(val));
		return ret_vecs;
	}

	template<class Type>
	block_multi_vector::ptr scale_cols(const std::vector<Type> &vec) const {
		size_t num_blocks = get_num_blocks();
		block_multi_vector::ptr ret_vecs = block_multi_vector::create(
				get_num_rows(), get_num_cols(), block_size, type);
		for (size_t i = 0; i < num_blocks; i++) {
			fm::mem_vector::ptr sub_vec = fm::mem_vector::create(block_size,
					fm::get_scalar_type<Type>());
			for (size_t k = 0; k < block_size; k++)
				sub_vec->set<Type>(k, vec[i * block_size + k]);
			ret_vecs->set_block(i, get_block(i)->scale_cols(*sub_vec));
		}
		return ret_vecs;
	}

	template<class Type>
	static void sparse_matrix_multiply(const fm::sparse_matrix &A,
			const block_multi_vector &X, block_multi_vector &Y) {
		sparse_matrix_multiply(sp_multiply_impl<Type>(), A, X, Y);
	}
};

#endif
