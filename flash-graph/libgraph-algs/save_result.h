#ifndef __SAVE_RESULT_H__
#define __SAVE_RESULT_H__

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

#include "FG_vector.h"
#include "graph_engine.h"

namespace {

template<class T, class VertexType>
class save_query: public fg::vertex_query
{
	typename fg::FG_vector<T>::ptr vec;
public:
	save_query(typename fg::FG_vector<T>::ptr vec) {
		this->vec = vec;
	}

	virtual void run(fg::graph_engine &graph, fg::compute_vertex &v1) {
		VertexType &v = (VertexType &) v1;
		vec->set(graph.get_graph_index().get_vertex_id(v), v.get_result());
	}

	virtual void merge(fg::graph_engine &graph, fg::vertex_query::ptr q) {
	}

	virtual ptr clone() {
		return fg::vertex_query::ptr(new save_query(vec));
	}
};

}

#endif
