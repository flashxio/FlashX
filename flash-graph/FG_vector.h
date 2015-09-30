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
#include <fstream>

#include "graph_engine.h"
#include "stat.h"

namespace fg
{

/**
 * \brief FlashGraph vector that provides several parallelized methods
 *        when compared to an STL-vector. <br>
 * **NOTE**: Not an STL-compatible data structure. This vector is also
 * ideally used with numeric data types. <br>
 * Methods marked with the keyword **parallel** are parallelized implementations.
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
	 * \brief  Create a vector of the length the same as the number of vertices
	 *         in the graph. An object of this
	 *         class should be created using this or the `create(size_t size)`
	 *         method.
	 * \param graph A shared pointer to a graph engine object. This is generally
	 *       the graph for which you are creating the vector.
	 */
	static ptr create(graph_engine::ptr graph) {
		return ptr(new FG_vector<T>(graph));
	}

	/**
	 * \brief  Create a vector of the specified length. An object of this
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
	 * \brief Equivalent to += operator. Element by element
	 *     addition of one `FG_vector` to another.
	 * \param other An `FG_vector` smart pointer object.
	 *
	 */
	void plus_eq(FG_vector<T>::ptr other) {
		assert(get_size() == other->get_size());
		for (size_t i = 0; i < get_size(); i++) {
			eles[i] += other->get(i);
		}
	}

	/**
	 * \brief Assign a value `num` many times to the vector.
	 * \param num The number of elements to assign.
	 * \param val The value a user wnats to assign to vector positions.
	 */
	void assign(size_t num, T val) {
		eles.assign(num, val);
	}

	/** 
	 * \brief Make a shallow copy of the vector.
	 * \param other An `FG_vector` smart pointer.
	 * **paralel**
	 */
	void shallow_copy(FG_vector<T>::ptr other) {
		assert(this->get_size() == other->get_size());
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++) {
			this->eles[i] = other->eles[i];
		}
	}

	template<class T1>
	void copy_to(T1 *arr, size_t size) {
		size_t num = std::min(size, eles.size());
		for (size_t i = 0; i < num; i++)
			arr[i] = eles[i];
	}

	/**
	 * \brief Check for equality between two `FG_vector`s element by
	 *   element.
	 * \param other An `FG_vector` smart pointer.
	 */
	// TODO DM: Make parallel / smarter 
	bool eq_all(FG_vector<T>::ptr other) {
		return std::equal(this->eles.begin(), this->eles.end(), other->eles.begin());
	}

	void init_rand(long max = std::numeric_limits<T>::max(),
			unsigned int seed = 0) {
		if (seed > 0)
			srandom(seed);
		if (max >= std::numeric_limits<T>::max())
			max = std::numeric_limits<T>::max();
#pragma omp parallel for
		for (size_t i = 0; i < eles.size(); i++)
			eles[i] = random() % max;
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
	 * \return The number of elements in the vector
	 */
	size_t get_size() const {
		return eles.size();
	}

	/**
	 * \brief  Get a pointer to the memory array used internally by
	 *         the vector to store its owned elements.
	 * \return A pointer the underlying data memory array.
	 *
	 */
	T *get_data() {
		return eles.data();
	}

	/**
	 * \brief  Const method to get a pointer to the memory array 
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
	 *        of two FG vectors. <br>
	 * **parallel**
	 *
	 * \return A value of data type `T` value that is the dot product.
	 */
	T dot_product(const FG_vector<T> &other) const {
		assert(this->get_size() == other.get_size());
		T ret = 0;
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * other.get(i);
		return ret;
	}

	/**
	 * \brief Compute the 
	 *        [L2 Norm](http://en.wikipedia.org/wiki/Norm_(mathematics)#Euclidean_norm)
	 *        (also know as Euclidean distance) of a vector. <br>
	 * **parallel**
	 *
	 * \return An object of type `T` with the value of the L2 norm. 
	 */
	T norm2() const {
		T ret = 0;
		for (size_t i = 0; i < get_size(); i++)
			ret += get(i) * get(i);
		return sqrt(ret);
	}

	/**
	 * \brief Compute the 
	 * [L1 Norm](http://en.wikipedia.org/wiki/Norm_(mathematics)#Taxicab_norm_or_Manhattan_norm) 
	 * (also Taxicab norm) of an FG_vector. <br>
	 * **parallel**
	 *
	 * \return An object of type `T` with the L1 norm.
	 */
	T norm1() const {
		T ret = 0;
		for (size_t i = 0; i < get_size(); i++)
			ret += fabs(get(i));
		return ret;
	}

	/**
	 * \brief Compute the sum of all elements in the vector. <br>
	 * If the type is integer, the sum can overflow.
	 * **parallel**
	 * \return The sum of all items in the vector.
	 */
	T sum() const {
		return sum<T>();
	}

	/**
	 * \brief Compute the sum of all elements in the vector. <br>
	 * This sum() allows users to specify the type of the result, so users
	 * can avoid integer overflow.
	 * **parallel**
	 * \return The sum of all items in the vector.
	 */
	template<class ResType>
	ResType sum() const {
		struct identity_func {
			ResType operator()(T v) {
				return v;
			}
		};
		return aggregate<identity_func, ResType>(identity_func());
	}

	template<class Func, class ResType>
	ResType aggregate(Func func) const {
		ResType ret = 0;
		for (size_t i = 0; i < get_size(); i++)
			ret += func(eles[i]);
		return ret;
	}

	/**
	 * \brief Find the maximal value in the vector and return its value.
	 * \return The maximal value in the vector.
	 */
	T max() const {
		return max_val_loc().first;
	}

	/**
	 * \brief Find the maximal value in the vector and return its value
	 *        and its location.
	 * \return A pair that contains the maximal value and its location
	 *         in the vector.
	 */
	std::pair<T, off_t> max_val_loc() const {
		T ret = std::numeric_limits<T>::min();
		off_t idx = 0;
		for (size_t i = 0; i < get_size(); i++) {
			if (ret < get(i)) {
				ret = get(i);
				idx = i;
			}
		}
		return std::pair<T, off_t>(ret, idx);
	}

	void max_val_locs(size_t num, std::vector<std::pair<T, off_t> > &pairs) const {
		typedef std::pair<T, off_t> val_loc_t;
		struct comp_val {
			bool operator()(const val_loc_t &v1, const val_loc_t &v2) {
				return v1.first > v2.first;
			}
		};
		std::priority_queue<val_loc_t, std::vector<val_loc_t>, comp_val> queue;
		for (size_t i = 0; i < get_size(); i++) {
			T val = get(i);
			queue.push(val_loc_t(val, i));
			if (queue.size() > num)
				queue.pop();
		}
		while (!queue.empty()) {
			val_loc_t pair = queue.top();
			queue.pop();
			pairs.push_back(pair);
		}
	}

	/**
	 * \brief Find the index with the minmimal value in the vector and
	 *     return its value.
	 * \return The minimal value in the vector.
	 */
	T min() const {
		T ret = std::numeric_limits<T>::max();
		for (size_t i = 0; i < get_size(); i++)
			ret = std::min(get(i), ret);
		return ret;
	}

	/**
	 * \brief Find the index with the minimal value in the vector and
	 *     return *the index*.
	 * \return The minimal index value in the vector.
	 */
	size_t argmin() {
		typename std::vector<T>::iterator res = std::min_element(eles.begin(), eles.end());
		size_t ret = std::distance(eles.begin(), res);
		return ret;
	}

	/**
	 * \brief Serial element-wise print of the vector.
	 *	**Not intended for very large vectors**
	 */
	void print(vsize_t max_print_size=100) {
		vsize_t print_len = get_size() > max_print_size ? max_print_size : get_size();

		std::cout << "[";
		for (vsize_t i=0; i < print_len; i++) {
			std::cout << " " << get(i);
		}

		if (print_len == max_print_size && print_len != get_size()) {
			std::cout << "...";
		}
		std::cout <<  " ]\n\n";
	}

	/**
	 * \brief Write the space separated vector to file.
	 * \param fn The file name you wish written to file.
	 */
	void to_file(std::string fn) {
		std::ofstream f;
		f.open(fn);
		for (vsize_t i=0; i < get_size(); i++) {
			f << get(i) << " ";
		}
		f.close();
	}

	void neg_in_place() {
		for (size_t i = 0; i < get_size(); i++)
			eles[i] = -eles[i];
	}

	/**
	 * \brief In place division of vector by a single value.
	 * \param v The value by which you want the array divided.
	 * **parallel** 
	 */
	void div_by_in_place(T v) {
#pragma omp parallel for
		for (size_t i = 0; i < get_size(); i++)
			eles[i] /= v;
	}

	/**
	 * \brief element-wise merge with another vector and store the result
	 *        in this vector.
	 * \param vec The vector that you want to merge with.
	 * \param func The operator that you want to perform on each pair of
	 *             elements.
	 */
	template<class MergeFunc, class VecType>
		void merge_in_place(typename FG_vector<VecType>::ptr vec, MergeFunc func) {
			assert(this->get_size() == vec->get_size());
#pragma omp parallel for
			for (size_t i = 0; i < get_size(); i++)
				eles[i] = func(eles[i], vec->get(i));
		}

	/**
	 * \brief In place element-wise addition by another vector.
	 * \param vec The vector by which you want to add to this vector.
	 * **parallel**
	 */
	template<class T2>
	void add_in_place(typename FG_vector<T2>::ptr vec) {
		struct add_func {
			T operator()(const T &v1, const T2 &v2) {
				return v1 + v2;
			}
		};
		merge_in_place<add_func, T2>(vec, add_func());
	}

	/**
	 * \brief In place subtraction of the vector by another vector.
	 * \param vec The vector by which you want the array to be subtracted.
	 * **parallel** 
	 */
	template<class T2>
	void subtract_in_place(typename FG_vector<T2>::ptr &vec) {
		struct sub_func {
			T operator()(const T &v1, const T2 &v2) {
				return v1 - v2;
			}
		};
		merge_in_place<sub_func, T2>(vec, sub_func());
	}

	template<class T2>
	void multiply_in_place(T2 v) {
		for (size_t i = 0; i < get_size(); i++)
			eles[i] *= v;
	}

	template<class IN_TYPE, class OUT_TYPE>
	typename FG_vector<OUT_TYPE>::ptr multiply(IN_TYPE v) const {
		typename FG_vector<OUT_TYPE>::ptr ret = FG_vector<OUT_TYPE>::create(get_size());
		for (size_t i = 0; i < get_size(); i++)
			ret->set(i, this->eles[i] * v);
		return ret;
	}

	template<class IN_TYPE, class OUT_TYPE>
	typename FG_vector<OUT_TYPE>::ptr multiply(typename FG_vector<IN_TYPE>::ptr vec) const {
		if (vec->get_size() != this->get_size())
			return typename FG_vector<OUT_TYPE>::ptr();

		typename FG_vector<OUT_TYPE>::ptr ret = FG_vector<OUT_TYPE>::create(get_size());
		for (size_t i = 0; i < get_size(); i++)
			ret->eles[i] = this->eles[i] * vec->get(i);
		return ret;
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
				ABORT_MSG("normalize on wrong type");
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

	log_histogram log_hist(int power) const {
		T max_v = max();
		int num_buckets = ceil(log(max_v) / log(power));
		log_histogram hist(std::max(num_buckets, 1));
		for (size_t i = 0; i < get_size(); i++) {
			hist.add_value(eles[i]);
		}
		return hist;
	}
};

/**
 * \brief Apply a user defined function to multipl FG_vectors.
 * **parallel**
 * \param inputs A vector of FG_vectors that are the inputs.
 * \param output A FG_vector that are the outputs.
 * \param apply  The user-defined function that will be applied to all vecotors.
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

}

#endif
