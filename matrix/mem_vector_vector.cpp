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

#include <string.h>
#include <assert.h>

#include "common.h"

#include "mem_vector_vector.h"
#include "mem_vector.h"
#include "data_frame.h"
#include "factor.h"
#include "local_vec_store.h"

namespace fm
{

vector::ptr mem_vector_vector::cat() const
{
	return mem_vector::create(detail::mem_vec_store::cast(store->cat()));
}

class sorted_gr_label_operate: public gr_apply_operate<local_vec_store>
{
	const factor_vector &labels;
	const gr_apply_operate<sub_vector_vector> &op;
	const mem_vector_vector &vv;
public:
	sorted_gr_label_operate(const factor_vector &_labels,
			const gr_apply_operate<sub_vector_vector> &_op,
			const mem_vector_vector &_vv): labels(_labels), op(_op), vv(_vv) {
	}

	virtual void run(const void *key, const local_vec_store &val,
			local_vec_store &vec) const;

	virtual const scalar_type &get_key_type() const {
		return labels.get_type();
	}

	virtual const scalar_type &get_output_type() const {
		return op.get_output_type();
	}

	virtual size_t get_num_out_eles() const {
		return 0;
	}
};

void sorted_gr_label_operate::run(const void *key, const local_vec_store &val,
		local_vec_store &vec) const
{
	off_t off_bytes = ((const char *) key) - labels.get_raw_arr();
	off_t off = off_bytes / labels.get_type().get_size();
	size_t len = val.get_length();
	std::vector<off_t> vec_idxs(len);
	for (size_t i = 0; i < len; i++)
		vec_idxs[i] = i + off;
	op.run(key, sub_vector_vector(vv, vec_idxs), vec);
}

vector_vector::ptr mem_vector_vector::groupby(const factor_vector &labels,
			const gr_apply_operate<sub_vector_vector> &op) const
{
	if (labels.get_length() != get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The label vector has a different length from the vector vector";
		return vector_vector::ptr();
	}
	assert(labels.is_sorted());

	// TODO this only works if the input labels have been sorted.
	data_frame::ptr df = labels.groupby(
			sorted_gr_label_operate(labels, op, *this), false);
	return mem_vector_vector::create(
			detail::mem_vv_store::cast(df->get_vec("agg")));
}

vector_vector::ptr mem_vector_vector::apply(const arr_apply_operate &op) const
{
	int num_parts = std::min((size_t) get_num_omp_threads(), get_num_vecs());
	std::vector<vector_vector::ptr> vvs(num_parts);
	size_t part_num_vecs = ceil(((double) get_num_vecs()) / num_parts);
#pragma omp parallel for
	for (int i = 0; i < num_parts; i++) {
		// TODO Here I divide the vector_vector evenly in terms of the number
		// of vectors in each partition. It potentially has load balancing
		// problems.
		off_t start = part_num_vecs * i;
		off_t end = std::min(part_num_vecs * (i + 1), get_num_vecs());
		// When there aren't many vectors, this may happen.
		if (end <= start)
			continue;
		mem_vector_vector::ptr sub_vv = mem_vector_vector::create(
				store->get_sub_vec_vec(start, end - start));
		vvs[i] = sub_vv->serial_apply(op);
	}

	if (num_parts == 1)
		return vvs[0];

	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(
			op.get_output_type());
	// It's possible that some threads didn't generate results.
	// It's guaranteed that the non-empty results are in the front.
	int num_non_empty = 0;
	std::vector<detail::vec_store::const_ptr> vv_stores(vvs.size());
	for (int i = 0; i < num_parts; i++) {
		if (vvs[i]) {
			assert(i == num_non_empty);
			num_non_empty++;
			vv_stores[i] = vvs[i]->get_raw_store();
		}
	}
	ret->append(vv_stores.begin(), vv_stores.begin() + num_non_empty);
	return mem_vector_vector::create(ret);
}

vector_vector::ptr mem_vector_vector::serial_apply(
		const arr_apply_operate &op) const
{
	const scalar_type &output_type = op.get_output_type();
	size_t out_size;
	// If the user can predict the number of output elements, we can create
	// a buffer of the expected size.
	if (op.get_num_out_eles() > 0)
		out_size = op.get_num_out_eles();
	else
		// If the user can't, we create a small buffer.
		out_size = 16;
	local_buf_vec_store buf(0, out_size, output_type, -1);

	detail::mem_vv_store::ptr ret = detail::mem_vv_store::create(output_type);
	for (size_t i = 0; i < get_num_vecs(); i++) {
		local_cref_vec_store lvec(this->get_raw_arr(i), 0,
				this->get_length(i), get_type(), -1);
		op.run(lvec, buf);
		ret->append(buf);
	}
	return mem_vector_vector::create(ret);
}

}
