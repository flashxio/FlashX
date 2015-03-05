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

#include "mem_vector.h"
#include "mem_data_frame.h"

namespace fm
{

static char *get_matrix_raw_data(mem_dense_matrix &data)
{
	char *arr = NULL;
	// TODO this may not work with submatrix.
	if (data.store_layout() == matrix_layout_t::L_ROW)
		arr = ((mem_row_dense_matrix &) data).get_row(0);
	else if (data.store_layout() == matrix_layout_t::L_COL)
		arr = ((mem_col_dense_matrix &) data).get_col(0);
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix layout";
	return arr;
}

mem_vector::mem_vector(mem_dense_matrix::ptr data): vector(
		// The length of the vector is the size of the dimension that isn't 1.
		data->get_num_rows() == 1 ? data->get_num_cols(): data->get_num_rows(),
		data->get_entry_size(), true)
{
	this->data = data;
	// The data buffer is the first row or column of the matrix, so
	// the shape of the matrix doesn't matter.
	arr = get_matrix_raw_data(*data);
}

mem_vector::mem_vector(std::shared_ptr<char> data, size_t length,
		size_t entry_size): vector(length, entry_size, true)
{
	// Maybe the column form may be more useful.
	mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(data, length,
			1, entry_size);
	this->arr = tmp->get_col(0);
	this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
}

mem_vector::mem_vector(size_t length, size_t entry_size): vector(length,
		entry_size, true)
{
	// Maybe the column form may be more useful.
	mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(length,
			1, entry_size);
	this->arr = tmp->get_col(0);
	this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
}

mem_vector::ptr mem_vector::cast(vector::ptr vec)
{
	if (!vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast a non-in-mem vector to in-mem vector";
		return mem_vector::ptr();
	}
	return std::static_pointer_cast<mem_vector>(vec);
}

mem_vector::const_ptr mem_vector::cast(vector::const_ptr vec)
{
	if (!vec->is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "can't cast a non-in-mem vector to in-mem vector";
		return mem_vector::const_ptr();
	}
	return std::static_pointer_cast<const mem_vector>(vec);
}

bool mem_vector::verify_groupby(const agg_operate &agg_op) const
{
	if (agg_op.input_entry_size() != get_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the operator accepts input types incompatible with the entry type in the vector";
		return false;
	}

	return true;
}

data_frame::ptr mem_vector::serial_groupby(const agg_operate &agg_op) const
{
	const scalar_type &output_type = agg_op.get_output_type();
	const agg_operate &find_next = get_type().get_agg_ops().get_find_next();
	if (!verify_groupby(agg_op))
		return data_frame::ptr();

	// If the vector hasn't been sorted, we need to sort it.
	vector::ptr vec;
	const mem_vector *sorted_vec;
	if (!is_sorted()) {
		// We don't want groupby changes the original vector.
		vec = this->deep_copy();
		vec->sort();
		sorted_vec = (const mem_vector *) vec.get();
	}
	else
		sorted_vec = this;

	size_t loc = 0;
	mem_vector::ptr agg = output_type.create_mem_vec(16);
	mem_vector::ptr val = this->create_int(16);
	size_t idx = 0;
	while (loc < sorted_vec->get_length()) {
		size_t curr_length = sorted_vec->get_length() - loc;
		const char *curr_ptr = sorted_vec->get_raw_arr() + get_entry_size() * loc;
		size_t rel_end;
		find_next.run(curr_length, curr_ptr, &rel_end);
		if (idx >= agg->get_length()) {
			agg->resize(agg->get_length() * 2);
			val->resize(val->get_length() * 2);
		}
		agg_op.run(rel_end, curr_ptr, agg->get(idx));
		memcpy(val->get(idx), curr_ptr, get_entry_size());
		idx++;
		loc += rel_end;
	}
	agg->resize(idx);
	val->resize(idx);
	mem_data_frame::ptr ret = mem_data_frame::create();
	ret->add_vec("val", std::static_pointer_cast<vector>(val));
	ret->add_vec("agg", std::static_pointer_cast<vector>(agg));
	return std::static_pointer_cast<data_frame>(ret);
}

data_frame::ptr mem_vector::groupby(const agg_operate &agg_op) const
{
	const agg_operate &find_next = get_type().get_agg_ops().get_find_next();
	if (!verify_groupby(agg_op))
		return data_frame::ptr();

	// If the vector hasn't been sorted, we need to sort it.
	vector::ptr vec;
	const mem_vector *sorted_vec;
	if (!is_sorted()) {
		// We don't want groupby changes the original vector.
		vec = this->deep_copy();
		vec->sort();
		sorted_vec = (const mem_vector *) vec.get();
	}
	else
		sorted_vec = this;

	// We need to find the start location for each thread.
	// The start location is where the value in the sorted array
	// first appears.
	int num_omp = get_num_omp_threads();
	off_t par_starts[num_omp + 1];
	for (int i = 0; i < num_omp; i++) {
		off_t start = sorted_vec->get_length() / num_omp * i;
		// This returns the relative start location of the next value.
		find_next.run(sorted_vec->get_length() - start,
				sorted_vec->get_raw_arr() + get_entry_size() * start,
				par_starts + i);
		// This is the absolute start location of this partition.
		par_starts[i] += start;
	}
	par_starts[num_omp] = sorted_vec->get_length();

	// It's possible that two partitions end up having the same start location
	// because the vector is small or a partition has only one value.
	assert(std::is_sorted(par_starts, par_starts + num_omp));
	off_t *end_par_starts = std::unique(par_starts, par_starts + num_omp);
	int num_parts = end_par_starts - par_starts;
	std::vector<data_frame::ptr> sub_results(num_parts);
	for (int i = 0; i < num_parts; i++) {
		off_t start = par_starts[i];
		off_t end = par_starts[i + 1];
		vector::const_ptr sub_vec = sorted_vec->get_sub_vec(start, end - start);
		sub_results[i] = mem_vector::cast(sub_vec)->serial_groupby(agg_op);
	}

	data_frame::ptr ret = sub_results[0];
	if (num_parts > 1)
		ret->append(sub_results.begin() + 1, sub_results.end());
	return ret;
}

bool mem_vector::append(std::vector<vector::ptr>::const_iterator vec_it,
		std::vector<vector::ptr>::const_iterator vec_end)
{
	// Get the total size of the result vector.
	size_t tot_res_size = this->get_length();
	for (auto it = vec_it; it != vec_end; it++) {
		tot_res_size += (*it)->get_length();
		if (!(*it)->is_in_mem()) {
			BOOST_LOG_TRIVIAL(error)
				<< "Not support appending an ext-mem vector to an in-mem vector";
			return false;
		}
	}

	// Merge all results to a single vector.
	off_t loc = this->get_length();
	this->resize(tot_res_size);
	for (auto it = vec_it; it != vec_end; it++) {
		assert(loc + (*it)->get_length() <= this->get_length());
		this->set_sub_vec(loc, **it);
		loc += (*it)->get_length();
	}
	return true;
}

bool mem_vector::append(vector::ptr vec)
{
	off_t loc = this->get_length();
	// TODO We might want to over expand a little, so we don't need to copy
	// the memory again and again.
	// TODO if this is a sub_vector, what should we do?
	this->resize(vec->get_length() + get_length());
	assert(loc + vec->get_length() <= this->get_length());
	this->set_sub_vec(loc, *vec);
	return true;
}

bool mem_vector::resize(size_t new_length)
{
	if (new_length == get_length())
		return true;

	size_t tot_len = data->get_num_rows() * data->get_num_cols();
	// We don't want to reallocate memory when shrinking the vector.
	if (new_length < tot_len) {
		return vector::resize(new_length);
	}

	// Keep the old information of the vector.
	mem_dense_matrix::ptr old_data = data;
	char *old_arr = arr;
	size_t old_length = get_length();

	// We realloate memory regardless of whether we increase or decrease
	// the length of the vector.
	mem_col_dense_matrix::ptr tmp = mem_col_dense_matrix::create(new_length,
			1, get_entry_size());
	if (tmp == NULL) {
		BOOST_LOG_TRIVIAL(error) << "can't allocate memory to resize the vector";
		return false;
	}

	this->arr = tmp->get_col(0);
	this->data = std::static_pointer_cast<mem_dense_matrix>(tmp);
	memcpy(arr, old_arr, std::min(old_length, new_length) * get_entry_size());
	return vector::resize(new_length);
}

bool mem_vector::set_sub_vec(off_t start, const vector &vec)
{
	if (!vec.is_in_mem()) {
		BOOST_LOG_TRIVIAL(error)
			<< "Not support setting a subvector from ext-mem vector";
		return false;
	}
	if (get_entry_size() != vec.get_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The two vectors don't have the same entry size";
		return false;
	}
	if (start + vec.get_length() > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "set_sub_vec: out of range";
		return false;
	}

	const mem_vector &mem_vec = (const mem_vector &) vec;
	memcpy(get_raw_arr() + start * get_entry_size(), mem_vec.get_raw_arr(),
			mem_vec.get_length() * mem_vec.get_entry_size());
	return true;
}

vector::const_ptr mem_vector::get_sub_vec(off_t start, size_t length) const
{
	if (start + length > get_length()) {
		BOOST_LOG_TRIVIAL(error) << "get_sub_vec: out of range";
		return vector::ptr();
	}

	// We need to discard the const from the "this" pointer.
	mem_vector *mutable_this = (mem_vector *) this;

	mem_vector::ptr mem_vec = mem_vector::cast(mutable_this->shallow_copy());
	mem_vec->resize(length);
	mem_vec->arr = mutable_this->get_raw_arr() + start * get_entry_size();
	mem_vec->data = mutable_this->get_data();
	return std::static_pointer_cast<const vector>(mem_vec);
}

bool mem_vector::expose_sub_vec(off_t start, size_t length)
{
	size_t tot_len = data->get_num_rows() * data->get_num_cols();
	if (start + length > tot_len) {
		exit(1);
		BOOST_LOG_TRIVIAL(error) << "expose_sub_vec: out of range";
		return false;
	}

	resize(length);
	arr = get_matrix_raw_data(*data) + start * get_entry_size();
	return true;
}

vector::ptr mem_vector::deep_copy() const
{
	// We need to discard the const from the "this" pointer.
	mem_vector *mutable_this = (mem_vector *) this;
	mem_vector::ptr mem_vec = mem_vector::cast(mutable_this->shallow_copy());
	mem_vec->data = mem_dense_matrix::cast(data->deep_copy());
	if (mem_vec->data->store_layout() == matrix_layout_t::L_ROW)
		mem_vec->arr = mem_row_dense_matrix::cast(mem_vec->data)->get_row(0);
	else if (mem_vec->data->store_layout() == matrix_layout_t::L_COL)
		mem_vec->arr = mem_col_dense_matrix::cast(mem_vec->data)->get_col(0);
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix layout";
	return std::static_pointer_cast<vector>(mem_vec);
}

bool mem_vector::equals(const mem_vector &vec) const
{
	if (vec.get_length() != this->get_length())
		return false;
	else
		return memcmp(this->arr, vec.arr,
				get_length() * get_entry_size()) == 0;
}

mem_vector::ptr mem_vector::get(type_mem_vector<off_t> &idxs) const
{
	mem_vector::ptr ret = create_int(idxs.get_length());
#pragma omp parallel for
	for (size_t i = 0; i < idxs.get_length(); i++) {
		off_t idx = idxs.get(i);
		// Check if it's out of the range.
		if (idx < 0 && idx >= this->get_length()) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("%1% is out of range") % idx;
			continue;
		}

		char *dst = ret->get(i);
		const char *src = this->get(idx);
		memcpy(dst, src, get_entry_size());
	}
	return ret;
}

bool mem_vector::export2(const std::string &file) const
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't open %1%, error: %1%") % file % strerror(errno);
		return false;
	}

	size_t ret = fwrite(arr, get_length(), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error) << boost::format(
				"can't write to %1%, error: %1%") % file % strerror(errno);
		return false;
	}

	fclose(f);
	return true;
}

template<>
vector::ptr create_vector<double>(double start, double end,
		double stride)
{
	// The result of division may generate a real number slightly smaller than
	// what we want because of the representation precision in the machine.
	// When we convert the real number to an integer, we may find the number
	// is smaller than exepcted. We need to add a very small number to
	// the real number to correct the problem.
	// TODO is it the right way to correct the problem?
	long n = (end - start) / stride + 1e-9;
	if (n < 0) {
		BOOST_LOG_TRIVIAL(error) <<"wrong sign in 'by' argument";
		return vector::ptr();
	}
	// We need to count the start element.
	n++;

	typename type_mem_vector<double>::ptr v
		= type_mem_vector<double>::create(n);
	v->get_data()->set_data(seq_set_operate<double>(n, start, stride));
	return std::static_pointer_cast<vector>(v);
}

}
