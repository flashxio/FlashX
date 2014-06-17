#ifndef __FG_VECTOR_H__
#define __FG_VECTOR_H__

/*
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

/**
  * \brief FlashGraph vector that provides several parallelized methods
  *        when compared to an STL-vector.
  * **NOTE**: Not an STL-compatible data structure. This vector is also
  * ideally used with numeric data types.
  * Methods marked with the keyword **parallel** are parallelized impls
 */
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
	typedef typename std::shared_ptr<FG_vector<T> > ptr; /** Smart pointer for object access */

  /** 
    * \brief  Create graph using this static method. An object of this
    *         class should be created using this or the `create(size_t size)`
    *         method.
    * \param graph A shared pointer to a graph engine object. This is generally
    *       the graph for which you are creating the vector. Allows the avoidance
    *       of computing the number of vertices in the graph.
  */
	static ptr create(graph_engine::ptr graph) {
		return ptr(new FG_vector<T>(graph));
	}

  /**
    * \brief  Create graph using this static method. An object of this
    *         class should be created using this or the `create(graph_engine::ptr graph)`
    *         method.
    * \param size The length of the vector you desire.
  */
	static ptr create(size_t size) {
		return ptr(new FG_vector<T>(size));
	}

  /**
    * \brief Initialize the vector a single value as specified by parameter 1.
    *
    * \param v The initialization parameter for the vector data.
    * **parallel**
  */
	void init(T v) {
#pragma omp parallel for
		for (size_t i = 0; i < eles.size(); i++)
			eles[i] = v;
	}

  /**
    * \brief  Populate an [STL set](http://www.cplusplus.com/reference/set/set/)
    *         with the unique elements in the vector. All duplicates are ignored.
    *
    * \param set The *empty* STL set that will be populated with unique vector members.
    *
  */
	void unique(std::set<T> &set) const {
		// TODO we need a parallel implementation.
    
    assert(set.empty()); // FIXME: `new` a shared/unique ptr & remove param
		BOOST_FOREACH(T v, eles) {
			set.insert(v);
		}
	}

  /** 
    * FIXME: Verify
    * \brief  Count the number of unique items in the vector using a
    *         count map.
    * \param map An *empty* `count_map` object that is used to count
    *         the number of unique elements in the vector.
    *
  */
	void count_unique(count_map<T> &map) const {
		// TODO we need a parallel implementation.
    
    assert(map.get_size() == 0); // FIXME: `new` a shared/unique ptr & remove param
		BOOST_FOREACH(T v, eles) {
			map.add(v);
		}
	}

  /**
    * \brief Get the number of elements contained in the vector.
    *
    * \return The count of the number of elements in the vector
  */
	size_t get_size() const {
		return eles.size();
	}

  /**
    * \brief  Get adirect pointer to the memory array used internally by
    *         the vector to store its owned elements.
    * \return A pointer the underlying data memory array.
    *
  */
	T *get_data() {
		return eles.data();
	}

  /**
    * \brief  Const method to get adirect pointer to the memory array 
    *         used internally by the vector to store its owned elements.
    * \return A const pointer the underlying data memory array.
    *
    *
  */
	const T*get_data() const {
		return eles.data();
	}

  /**
    * \brief Compute the [dot product](http://en.wikipedia.org/wiki/Dot_product)
    *        of two FG vectors.
    * **parallel**
    *
    * \return A value of data type `T` value that is the dot product.
  */
	T dot_product(const FG_vector<T> &other) const {
		assert(this->get_size() == other.get_size());
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * other.get(i);
		return ret;
	}

  /**
    * \brief Compute the 
    *        [L2 Norm](http://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm) 
    *        (also know as Euclidean distance) of a vector.
    * **parallel**
    *
    * \return An object of type `T` with the value of the L2 norm. 
  */
	T norm2() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * get(i);
		return sqrt(ret);
	}

  /**
    * \brief Compute the 
    * [L1 Norm](http://en.wikipedia.org/wiki/Norm_(mathematics)#Taxicab_norm_or_Manhattan_norm) 
    * (also Taxicab norm) of an FG_vector.
    * **parallel**
    *
    * \return An object of type `T` with the L1 norm.
  */
	T norm1() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += fabs(get(i));
		return ret;
	}

  /**
    * \brief Compute the sum of all elements in the vector.
    * **parallel**
    * \return The sum of all items in the vector.
  */
	T sum() const {
		T ret = 0;
#pragma omp parallel for reduction(+:ret)
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i);
		return ret;
	}

  /**
    * \brief Find the index with the maximal value in the vector and
    *     return its value.
    * \return The maximal value in the vector.
  */
	T max() const {
		T ret = std::numeric_limits<T>::min();
		for (size_t i = 0; i < get_size(); i++)
			ret = ::max(get(i), ret);
		return ret;
	}

  /**
    * \brief In place division of vector by a single value.
    * \param v the value by which you want the array divided by
    * **parallel** 
  */
	void div_by_in_place(T v) {
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++)
			eles[i] /= v;
	}

  /**
    * \brief Normalize vector using an Lx form.
    * **parallel**
  */
	void normalize(int type) {
		T norm;
		switch(type) {
			case 2:
				norm = norm2();
				break;
      case 1:
        norm = norm1();
        break;
			default:
				assert(0);
		}
		div_by_in_place(norm);
	}

  /**
    * \brief Apply a function to every element in an FG_vector.
    *
    * \param func A user-defined function.
    * \param output The FG_vector that you want to apply the function to.
    *
    * **parallel**
  */
	template<class ApplyFunc>
	void apply(ApplyFunc func, FG_vector<T> &output) {
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++)
			output.set(i, func(eles[i]));
	}

	// TODO these interfaces assume shared memory.

  /**
    * Set a value of an index in the vector.
    *
    * **NOTE:** This function assumes a shared memory environment.
    * \param id The index where value is being set.
    * \param v The value that the index will be set to.
  */
	void set(vertex_id_t id, const T &v) {
		eles[id] = v;
	}

  /**
    * \brief Const get the value of a particular index.
    * \param id The index of the vector from where you want a value.
    * \return The value requested by param 1
    *
  */
	const T &get(vertex_id_t id) const {
		return eles[id];
	}

  /**
    * \brief Non-const get the value of a particular index.
    * \param id The index of the vector from where you want a value.
    * \return The value requested by param 1
    *
  */
	T &get(vertex_id_t id) {
		return eles[id];
	}
};

  /**
    * \brief Apply a user defined function to an STL vector of FG_vectors.
    * **parallel**
    * \param inputs A vector of FG_vectors that are the inputs.
    * \param output A vector of FG_vectors that are the outputs.
    * \param apply  The user-defined function that will be applied to all vecotros.
  */
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
