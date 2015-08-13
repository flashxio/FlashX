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

#include "dense_matrix.h"
#include "generic_type.h"
#include "vector.h"
#include "sparse_matrix.h"

#include "eigensolver.h"

namespace fm
{

namespace eigen
{

extern size_t num_col_writes;
extern size_t num_col_writes_concept;
extern size_t num_col_reads_concept;
extern size_t num_multiply_concept;
extern dense_matrix::ptr cached_mat;

class block_multi_vector
{
	size_t MAX_MUL_BLOCKS;

	bool in_mem;
	size_t block_size;
	size_t num_rows;
	size_t num_cols;
	const fm::scalar_type &type;
protected:
	std::vector<fm::dense_matrix::ptr> mats;

	block_multi_vector(size_t nrow, size_t ncol, size_t block_size,
			const fm::scalar_type &_type, bool in_mem);
	block_multi_vector(const std::vector<fm::dense_matrix::ptr> &mats,
			bool in_mem);
public:
	typedef std::shared_ptr<block_multi_vector> ptr;

	static void sparse_matrix_multiply(const spm_function &multiply,
			const block_multi_vector &X, block_multi_vector &Y);

	static ptr create(size_t nrow, size_t ncol, size_t block_size,
			const fm::scalar_type &type, bool in_mem) {
		assert(ncol % block_size == 0);
		return ptr(new block_multi_vector(nrow, ncol, block_size, type, in_mem));
	}

	void set_multiply_blocks(size_t num) {
		MAX_MUL_BLOCKS = num;
	}

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

	fm::dense_matrix::ptr get_block(off_t block_idx) const {
		return mats[block_idx];
	}

	virtual void set_block(off_t block_idx, fm::dense_matrix::const_ptr mat) {
		assert(mat->get_num_cols() == block_size);
		mats[block_idx] = mat->clone();
	}

	void set_block(const block_multi_vector &mv, const std::vector<int>& index);

	fm::dense_matrix::const_ptr get_col(off_t col_idx) const;

	fm::dense_matrix::ptr get_col_mat(const std::vector<off_t> &index) const;

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
			fm::detail::mem_col_matrix_store::const_ptr B,
			const fm::scalar_variable &alpha,
			const fm::scalar_variable &beta) const;

	void assign(const block_multi_vector &vecs);

	block_multi_vector::ptr add(const block_multi_vector &vecs) const;
	fm::dense_matrix::ptr MvTransMv(const block_multi_vector &mv) const;
	std::vector<double> MvDot(const block_multi_vector &mv) const;

	fm::dense_matrix::ptr conv2matrix() const;

	template<class Type>
	void init_rand(Type min, Type max) {
		size_t num_blocks = get_num_blocks();
		num_col_writes += this->get_num_cols();
		for (size_t i = 0; i < num_blocks; i++)
			set_block(i, fm::dense_matrix::create_randu<Type>(min, max,
						get_num_rows(), block_size, fm::matrix_layout_t::L_COL,
						get_num_nodes(), in_mem));
	}

	template<class Type>
	block_multi_vector::ptr multiply_scalar(Type val) const {
		size_t num_blocks = get_num_blocks();
		block_multi_vector::ptr ret_vecs = block_multi_vector::create(
				get_num_rows(), get_num_cols(), block_size, type, in_mem);
		for (size_t i = 0; i < num_blocks; i++)
			ret_vecs->set_block(i, get_block(i)->multiply_scalar(val));
		return ret_vecs;
	}

	template<class Type>
	block_multi_vector::ptr scale_cols(const std::vector<Type> &vec) const {
		size_t num_blocks = get_num_blocks();
		block_multi_vector::ptr ret_vecs = block_multi_vector::create(
				get_num_rows(), get_num_cols(), block_size, type, in_mem);
		for (size_t i = 0; i < num_blocks; i++) {
			fm::detail::smp_vec_store::ptr sub_vec
				= fm::detail::smp_vec_store::create(block_size,
						fm::get_scalar_type<Type>());
			for (size_t k = 0; k < block_size; k++)
				sub_vec->set<Type>(k, vec[i * block_size + k]);
			fm::vector::ptr sub_vec1 = fm::vector::create(sub_vec);
			ret_vecs->set_block(i, get_block(i)->scale_cols(sub_vec1));
		}
		return ret_vecs;
	}
};

}

}

#endif
