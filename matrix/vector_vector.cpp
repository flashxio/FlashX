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

#include "vector.h"
#include "data_frame.h"
#include "factor.h"
#include "local_vec_store.h"
#include "local_vv_store.h"
#include "vector_vector.h"
#include "mem_vv_store.h"

namespace fm
{

vector_vector::ptr vector_vector::create(const detail::simple_raw_array &data,
		const std::vector<off_t> &offs, const scalar_type &type)
{
	detail::mem_vv_store::ptr vec = detail::mem_vv_store::create(data,
			offs, type);
	return ptr(new vector_vector(vec));
}

vector::ptr vector_vector::cat() const
{
	const detail::vv_store &store
		= static_cast<const detail::vv_store &>(get_data());
	return vector::create(store.cat());
}

class vv_gr_label_operate: public gr_apply_operate<sub_data_frame>
{
	const scalar_type &key_type;
	const gr_apply_operate<local_vv_store> &op;
public:
	vv_gr_label_operate(const scalar_type &type,
			const gr_apply_operate<local_vv_store> &_op): key_type(
				type), op(_op) {
	}

	virtual void run(const void *key, const sub_data_frame &val,
			local_vec_store &vec) const;

	virtual const scalar_type &get_key_type() const {
		return key_type;
	}

	virtual const scalar_type &get_output_type() const {
		return op.get_output_type();
	}

	virtual size_t get_num_out_eles() const {
		return 0;
	}
};

void vv_gr_label_operate::run(const void *key, const sub_data_frame &val,
		local_vec_store &vec) const
{
	assert(val.size() == 2);
	assert(val[0]->get_type() == key_type);
	assert(val[1]->get_entry_size() == 0);
	op.run(key, static_cast<const local_vv_store &>(*val[1]), vec);
}

vector_vector::ptr vector_vector::groupby(const factor_vector &labels,
			const gr_apply_operate<local_vv_store> &op) const
{
	struct vv_deleter {
		void operator()(detail::vec_store *vv) { }
	};
	if (labels.get_length() != get_num_vecs()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The label vector has a different length from the vector vector";
		return vector_vector::ptr();
	}
	assert(labels.is_sorted());

	data_frame::ptr df = data_frame::create();
	bool ret = df->add_vec("factor", const_cast<detail::vec_store *>(
				labels.get_raw_store().get())->shallow_copy());
	if (!ret)
		return vector_vector::ptr();
	ret = df->add_vec("vv", std::shared_ptr<detail::vec_store>(
				const_cast<detail::vec_store *>(&get_data()), vv_deleter()));
	if (!ret)
		return vector_vector::ptr();
	vv_gr_label_operate vv_op(labels.get_type(), op);
	return df->groupby("factor", vv_op);
}

vector_vector::ptr vector_vector::apply(const arr_apply_operate &op) const
{
	int num_parts = std::min((size_t) get_num_omp_threads(), get_num_vecs());
	std::vector<detail::vv_store::ptr> vv_stores(num_parts);
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
		local_vv_store::const_ptr lstore = local_vv_store::cast(
				get_data().get_portion(start, end - start));
		vv_stores[i] = detail::apply(*lstore, op);
	}

	detail::vv_store::ptr ret = detail::vv_store::create(
			op.get_output_type(), get_data().is_in_mem());
	// It's possible that some threads didn't generate results.
	// It's guaranteed that the non-empty results are in the front.
	int num_non_empty = 0;
	std::vector<detail::vec_store::const_ptr> tmp(vv_stores.size());
	for (int i = 0; i < num_parts; i++) {
		if (vv_stores[i]) {
			assert(i == num_non_empty);
			num_non_empty++;
			tmp[i] = vv_stores[i];
		}
	}
	ret->append(tmp.begin(), tmp.begin() + num_non_empty);
	return vector_vector::create(ret);
}

}
