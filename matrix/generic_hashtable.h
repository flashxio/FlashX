#ifndef __GENERIC_HASHTABLE_H__
#define __GENERIC_HASHTABLE_H__

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

#include<unordered_map>

#include "data_frame.h"
#include "mem_vec_store.h"
#include "bulk_operate_ext.h"

namespace fm
{

class agg_operate;

class generic_hashtable
{
public:
	typedef std::shared_ptr<generic_hashtable> ptr;

	/*
	 * This method inserts an array of keys and values.
	 * If a key exists, merge the old value and the new value with
	 * the provided operator.
	 */
	virtual void insert(size_t num, const void *keys, const void *vals,
			const agg_operate &op) = 0;
	/*
	 * Merge the input hashtable to this table.
	 * If a key exists in the original table, merge the corresponding values
	 * with the provided operator.
	 */
	virtual void merge(const generic_hashtable &table, const agg_operate &op) = 0;
	/*
	 * Convert the hashtable to a data frame.
	 */
	virtual data_frame::ptr conv2df() const = 0;
};

template<class KeyType, int ValSize>
class generic_hashtable_impl: public generic_hashtable
{
	const scalar_type &real_val_type;

	/*
	 * This isn't the real value type of the hashtable, but it has the same
	 * size as the real value, so we can use it to create a C++ hashtable.
	 */
	typedef struct {
		char data[ValSize];
	} ValType;
	std::unordered_map<KeyType, ValType> table;
public:
	generic_hashtable_impl(const scalar_type &type): real_val_type(type) {
	}

	void insert(size_t num, const void *pkeys, const void *pvals,
			const agg_operate &op) {
		assert(op.get_input_type().get_size() == ValSize);
		const KeyType *keys = (const KeyType *) pkeys;
		const ValType *vals = (const ValType *) pvals;
		assert(op.has_combine());
		for (size_t i = 0; i < num; i++) {
			auto ret = table.insert(std::pair<KeyType, ValType>(keys[i],
						vals[i]));
			// TODO this is going to be slow.
			if (!ret.second) {
				ValType tmp_vals[2];
				tmp_vals[0] = ret.first->second;
				tmp_vals[1] = vals[i];
				op.runCombine(2, tmp_vals, &ret.first->second);
			}
		}
	}

	virtual void merge(const generic_hashtable &gtable, const agg_operate &op) {
		assert(op.get_input_type().get_size() == ValSize);
		const generic_hashtable_impl<KeyType, ValSize> &gtable1
			= dynamic_cast<const generic_hashtable_impl<KeyType, ValSize> &>(
					gtable);
		assert(op.has_combine());
		for (auto it = gtable1.table.begin(); it != gtable1.table.end(); it++) {
			auto ret = table.insert(std::pair<KeyType, ValType>(it->first,
						it->second));
			if (!ret.second) {
				ValType tmp_vals[2];
				tmp_vals[0] = ret.first->second;
				tmp_vals[1] = it->second;
				op.runCombine(2, tmp_vals, &ret.first->second);
			}
		}
	}

	virtual data_frame::ptr conv2df() const {
		detail::smp_vec_store::ptr keys = detail::smp_vec_store::create(table.size(),
				get_scalar_type<KeyType>());
		detail::smp_vec_store::ptr vals = detail::smp_vec_store::create(table.size(),
				real_val_type);
		size_t vec_idx = 0;
		for (auto it = table.begin(); it != table.end(); it++) {
			keys->set<KeyType>(vec_idx, it->first);
			vals->set<ValType>(vec_idx, it->second);
			vec_idx++;
		}
		data_frame::ptr ret = data_frame::create();
		ret->add_vec("key", keys);
		ret->add_vec("val", vals);
		return ret;
	}
};

}

#endif
