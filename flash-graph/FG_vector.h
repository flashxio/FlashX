#ifndef __FG_VECTOR_H__
#define __FG_VECTOR_H__

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

#include <memory>
#include <set>

#include "graph_engine.h"
#include "stat.h"

template<class T>
class FG_vector
{
	// TODO I might need to split the vector into partitions.
	std::vector<T> eles;

	FG_vector(graph_engine::ptr graph) {
		eles.resize(graph->get_num_vertices());
	}

	FG_vector(size_t size) {
		eles.resize(size);
	}

public:
	typedef typename std::shared_ptr<FG_vector<T> > ptr;

	static ptr create(graph_engine::ptr graph) {
		return ptr(new FG_vector<T>(graph));
	}

	static ptr create(size_t size) {
		return ptr(new FG_vector<T>(size));
	}

	void init(T v) {
#pragma omp parallel for
		for (size_t i = 0; i < eles.size(); i++)
			eles[i] = v;
	}

	void unique(std::set<T> &set) const {
		// TODO we need a parallel implementation.
		BOOST_FOREACH(T v, eles) {
			set.insert(v);
		}
	}

	void count_unique(count_map<T> &map) const {
		// TODO we need a parallel implementation.
		BOOST_FOREACH(T v, eles) {
			map.add(v);
		}
	}

	size_t get_size() const {
		return eles.size();
	}

	T dot_product(const FG_vector<T> &other) const {
		assert(this->get_size() == other.get_size());
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * other.get(i);
		return ret;
	}

	T norm2() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * get(i);
		return sqrt(ret);
	}

	T norm1() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += fabs(get(i));
		return ret;
	}

	T sum() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i);
		return ret;
	}

	T max() const {
		T ret = std::numeric_limits<T>::min();
		for (size_t i = 0; i < get_size(); i++)
			ret = ::max(get(i), ret);
		return ret;
	}

	void div_by_in_place(T v) {
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++)
			eles[i] /= v;
	}

	void normalize(int type) {
		T norm;
		switch(type) {
			case 2:
				norm = norm2();
				break;
			default:
				assert(0);
		}
		div_by_in_place(norm);
	}

	template<class ApplyFunc>
	void apply(ApplyFunc func, FG_vector<T> &output) {
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++)
			output.set(i, func(eles[i]));
	}

	// TODO these interfaces assume shared memory.

	void set(vertex_id_t id, const T &v) {
		eles[id] = v;
	}

	const T &get(vertex_id_t id) const {
		return eles[id];
	}

	T &get(vertex_id_t id) {
		return eles[id];
	}
};

template<class T, class ApplyFunc>
void multi_vec_apply(const std::vector<typename FG_vector<T>::ptr> &inputs,
		typename FG_vector<T>::ptr output, ApplyFunc apply)
{
	for (size_t i = 0; i < inputs.size(); i++)
		assert(output->get_size() == inputs[i]->get_size());
#pragma omp parallel for
	for (size_t i = 0; i < output->get_size(); i++)
		output->set(i, apply(i, inputs));
}

#endif
