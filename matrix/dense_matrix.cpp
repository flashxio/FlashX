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

#include "log.h"

#include "dense_matrix.h"
#include "bulk_operate.h"
#include "mem_dense_matrix.h"
#include "EM_dense_matrix.h"
#include "generic_type.h"
#include "rand_gen.h"

namespace fm
{

bool dense_matrix::verify_inner_prod(const dense_matrix &m,
		const bulk_operate &left_op, const bulk_operate &right_op) const
{
	if (this->get_entry_size() != left_op.left_entry_size()
			|| m.get_entry_size() != left_op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The left operator isn't compatible with input matrices";
		return false;
	}

	if (left_op.output_entry_size() != right_op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The type of the left operator doesn't match the right operator";
		return false;
	}

	if (right_op.left_entry_size() != right_op.right_entry_size()
			|| right_op.left_entry_size() != right_op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output of the right operator has different types";
		return false;
	}

	if (get_num_cols() != m.get_num_rows()) {
		BOOST_LOG_TRIVIAL(error) << "The matrix size doesn't match";
		return false;
	}
	return true;
}

bool dense_matrix::verify_aggregate(const bulk_operate &op) const
{
	if (op.left_entry_size() != op.right_entry_size()
			|| op.left_entry_size() != op.output_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The input and output type of the operator is different";
		return false;
	}

	if (this->get_entry_size() != op.left_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The matrix entry size is different from the operator";
		return false;
	}
	return true;
}

bool dense_matrix::verify_mapply2(const dense_matrix &m,
			const bulk_operate &op) const
{
	if (this->get_num_rows() != m.get_num_rows()
			|| this->get_num_cols() != m.get_num_cols()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same shape";
		return false;
	}

	if (this->store_layout() != m.store_layout()) {
		BOOST_LOG_TRIVIAL(error)
			<< "two matrices in mapply2 don't have the same data layout";
		return false;
	}

	if (get_entry_size() != op.left_entry_size()
			|| m.get_entry_size() != op.right_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the element type in the matrices isn't compatible with the operator";
		return false;
	}

	return true;
}

bool dense_matrix::verify_apply(apply_margin margin, const arr_apply_operate &op) const
{
	if (get_entry_size() != op.input_entry_size()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the element type in the matrices isn't compatible with the operator";
		return false;
	}

	return true;
}

dense_matrix::ptr dense_matrix::create(size_t nrow, size_t ncol,
		const scalar_type &type, matrix_layout_t layout, bool in_mem)
{
	if (in_mem) {
		if (layout == matrix_layout_t::L_ROW)
			return mem_row_dense_matrix::create(nrow, ncol, type);
		else
			return mem_col_dense_matrix::create(nrow, ncol, type);
	}
	else {
		if (layout == matrix_layout_t::L_ROW) {
			fprintf(stderr, "EM row dense matrix isn't supported\n");
			return dense_matrix::ptr();
		}
		else
			return EM_col_dense_matrix::create(nrow, ncol, type);
	}
}

bool dense_matrix::write_header(FILE *f) const
{
	matrix_header header(DENSE, get_entry_size(), get_num_rows(),
			get_num_cols(), store_layout(), get_type().get_type());

	size_t ret = fwrite(&header, sizeof(header), 1, f);
	if (ret == 0) {
		BOOST_LOG_TRIVIAL(error)
			<< boost::format("can't write header: %2%") % strerror(errno);
		return false;
	}
	return true;
}

dense_matrix::ptr dense_matrix::load(const std::string &file_name)
{
	matrix_header header;

	FILE *f = fopen(file_name.c_str(), "r");
	if (f == NULL) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << boost::format("can't open %1%: %2%")
			% file_name % strerror(errno);
		return dense_matrix::ptr();
	}

	size_t ret = fread(&header, sizeof(header), 1, f);
	if (ret == 0) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << boost::format("can't read header from %1%: %2%")
			% file_name % strerror(errno);
		return dense_matrix::ptr();
	}

	header.verify();
	if (header.is_sparse()) {
		fclose(f);
		BOOST_LOG_TRIVIAL(error) << "The matrix to be loaded is sparse";
		return dense_matrix::ptr();
	}

	size_t nrow = header.get_num_rows();
	size_t ncol = header.get_num_cols();
	dense_matrix::ptr m;
	if (header.get_layout() == matrix_layout_t::L_ROW)
		m = std::static_pointer_cast<dense_matrix>(
				mem_row_dense_matrix::create(nrow, ncol,
					get_scalar_type(header.get_data_type()), f));
	else if (header.get_layout() == matrix_layout_t::L_COL)
		m = std::static_pointer_cast<dense_matrix>(
				mem_col_dense_matrix::create(nrow, ncol,
					get_scalar_type(header.get_data_type()), f));
	else
		BOOST_LOG_TRIVIAL(error) << "wrong matrix data layout";

	fclose(f);
	return m;
}

double dense_matrix::norm2() const
{
	// TODO this is an inefficient implementation.
	dense_matrix::ptr sq_mat = this->mapply2(
			*this, get_type().get_basic_ops().get_multiply());
	scalar_variable::ptr res = sq_mat->aggregate(
			sq_mat->get_type().get_basic_ops().get_add());
	double ret = 0;
	res->get_type().get_basic_uops().get_op(
			basic_uops::op_idx::SQRT)->runA(1, res->get_raw(), &ret);
	return ret;
}

namespace
{

/*
 * This class set elements in a container randomly.
 * set_operate can't change its own state and has to be thread-safe when
 * running on multiple threads. However, random generators aren't
 * thread-safe, so we have to create a random generator for each thread.
 */
class rand_init: public set_operate
{
	class rand_gen_wrapper {
		rand_gen::ptr gen;
	public:
		rand_gen_wrapper(rand_gen::ptr gen) {
			this->gen = gen;
		}

		rand_gen &get_gen() {
			return *gen;
		}
	};

	pthread_key_t gen_key;
	const scalar_type &type;
	const scalar_variable &min;
	const scalar_variable &max;

	rand_gen &get_rand_gen() const {
		void *addr = pthread_getspecific(gen_key);
		if (addr == NULL)
			addr = new rand_gen_wrapper(type.create_rand_gen(min, max));
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) addr;
		return wrapper->get_gen();
	}

	static void destroy_rand_gen(void *gen) {
		rand_gen_wrapper *wrapper = (rand_gen_wrapper *) gen;
		delete wrapper;
		printf("destroy rand gen\n");
	}
public:
	rand_init(const scalar_type &_type, const scalar_variable &_min,
			const scalar_variable &_max): type(_type), min(_min), max(_max) {
		int ret = pthread_key_create(&gen_key, destroy_rand_gen);
		assert(ret == 0);
	}

	~rand_init() {
		pthread_key_delete(gen_key);
	}

	virtual void set(void *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		get_rand_gen().gen(arr, num_eles);
	}
	virtual const scalar_type &get_type() const {
		return get_rand_gen().get_type();
	}
};

}

void dense_matrix::_init_rand(const scalar_variable &min,
		const scalar_variable &max)
{
	set_data(rand_init(get_type(), min, max));
}

}
