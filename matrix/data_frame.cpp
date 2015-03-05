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

#include <boost/format.hpp>

#include "log.h"

#include "data_frame.h"
#include "bulk_operate.h"

namespace fm
{

bool data_frame::append(std::vector<data_frame::ptr>::const_iterator begin,
		std::vector<data_frame::ptr>::const_iterator end)
{
	std::unordered_map<std::string, std::vector<vector::ptr> > vecs;
	for (size_t i = 0; i < named_vecs.size(); i++) {
		std::vector<vector::ptr> _vecs;
		_vecs.push_back(get_vec(named_vecs[i].first));
		vecs.insert(std::pair<std::string, std::vector<vector::ptr> >(
					named_vecs[i].first, _vecs));
	}

	for (auto it = begin; it != end; it++) {
		data_frame::ptr df = *it;
		if (df->get_num_vecs() != get_num_vecs()) {
			BOOST_LOG_TRIVIAL(error)
				<< "The data frames have different numbers of vectors";
			return false;
		}
		for (auto vec_it = vecs.begin(); vec_it != vecs.end(); vec_it++) {
			vector::ptr vec = df->get_vec(vec_it->first);
			if (vec == NULL) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different names on the vectors";
				return false;
			}
			if (vec->get_type() != vec_it->second.front()->get_type()) {
				BOOST_LOG_TRIVIAL(error)
					<< "The data frames have different types for the vectors with the same name";
				return false;
			}
			vec_it->second.push_back(vec);
		}
	}

	for (auto it = vecs.begin(); it != vecs.end(); it++)
		it->second.front()->append(it->second.begin() + 1, it->second.end());
	return true;
}

bool data_frame::append(data_frame::ptr df)
{
	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++) {
		if (df->get_vec(it->first) == NULL) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("The new data frame doesn't have column %1%")
				% it->first;
			return false;
		}
	}

	for (auto it = named_vecs.begin(); it != named_vecs.end(); it++)
		it->second->append(df->get_vec(it->first));
	return true;
}

bool data_frame::expose_portion(off_t loc, size_t length)
{
	for (size_t i = 0; i < named_vecs.size(); i++) {
		bool ret = named_vecs[i].second->expose_sub_vec(loc, length);
		if (!ret)
			return false;
	}
	return true;
}

}
