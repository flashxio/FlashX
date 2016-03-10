#ifndef __FM_COMBINED_MATRIX_STORE_H__
#define __FM_COMBINED_MATRIX_STORE_H__

/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "mapply_matrix_store.h"

namespace fm
{

namespace detail
{

/*
 * This class combines multiple matrices into a single one.
 * It always combines matrices in the shorter dimension.
 * e.g., given three n x m matrices, it generates a 3n x m matrix if n > m.
 * and a n x 3m matrix if m > n.
 * This matrix store supports writes.
 */
class combined_matrix_store: public mapply_matrix_store
{
	std::vector<matrix_store::const_ptr> mats;
protected:
	combined_matrix_store(const std::vector<matrix_store::const_ptr> &mats,
			matrix_layout_t layout);
public:
	typedef std::shared_ptr<combined_matrix_store> ptr;
	typedef std::shared_ptr<const combined_matrix_store> const_ptr;

	static ptr cast(matrix_store::ptr mat) {
		return std::dynamic_pointer_cast<combined_matrix_store>(mat);
	}

	static const_ptr cast(matrix_store::const_ptr mat) {
		return std::dynamic_pointer_cast<const combined_matrix_store>(mat);
	}

	static ptr create(const std::vector<matrix_store::const_ptr> &mats,
			matrix_layout_t layout);

	virtual std::string get_name() const;

	virtual matrix_store::const_ptr transpose() const;
	virtual matrix_store::const_ptr get_rows(const std::vector<off_t> &idxs) const;
	virtual matrix_store::const_ptr get_cols(
			const std::vector<off_t> &idxs) const {
		return matrix_store::get_cols(idxs);
	}

	size_t get_num_mats() const {
		return mats.size();
	}

	matrix_store::const_ptr get_mat(size_t off) const {
		return mats[off];
	}

	const matrix_store &get_mat_ref(size_t off) const {
		return *mats[off];
	}
};
}

}
#endif
