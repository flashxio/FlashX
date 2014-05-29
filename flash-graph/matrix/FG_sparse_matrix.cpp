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

#include "FG_sparse_matrix.h"

void FG_sparse_matrix::group_by_mean(const FG_vector<int> &labels,
		bool row_wise, std::map<int, FG_vector<double>::ptr> &agg_results)
{
	group_by<double, double>(labels, row_wise, agg_results);

	count_map<int> cmap;
	labels.count_unique(cmap);

	typedef std::map<int, typename FG_vector<double>::ptr> sum_map_t;
	BOOST_FOREACH(typename sum_map_t::value_type v, agg_results) {
		int label = v.first;
		int count = cmap.get(label);
		FG_vector<double>::ptr sum = v.second;
		for (size_t i = 0; i < sum->get_size(); i++)
			sum->set(i, sum->get(i) / count);
	}
}
