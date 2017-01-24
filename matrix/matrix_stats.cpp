/*
 * Copyright 2016 Open Connectome Project (http://openconnecto.me)
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

#include "matrix_stats.h"

namespace fm
{

namespace detail
{

size_t matrix_stats_t::inc_read_bytes(size_t bytes, bool in_mem)
{
#ifdef MATRIX_DEBUG
	if (in_mem) {
		mem_read_bytes += bytes;
		return mem_read_bytes;
	}
	else {
		EM_read_bytes += bytes;
		return EM_read_bytes;
	}
#else
	return 0;
#endif
}

size_t matrix_stats_t::get_read_bytes(bool in_mem) const
{
#ifdef MATRIX_DEBUG
	if (in_mem)
		return mem_read_bytes;
	else
		return EM_read_bytes;
#else
	return 0;
#endif
}

size_t matrix_stats_t::inc_write_bytes(size_t bytes, bool in_mem)
{
#ifdef MATRIX_DEBUG
	if (in_mem) {
		mem_write_bytes += bytes;
		return mem_write_bytes;
	}
	else {
		EM_write_bytes += bytes;
		return EM_write_bytes;
	}
#else
	return 0;
#endif
}

size_t matrix_stats_t::get_write_bytes(bool in_mem) const
{
#ifdef MATRIX_DEBUG
	if (in_mem)
		return mem_write_bytes;
	else
		return EM_write_bytes;
#else
	return 0;
#endif
}

size_t matrix_stats_t::inc_multiplies(size_t multiplies)
{
#ifdef MATRIX_DEBUG
	this->double_multiplies += multiplies;
	return double_multiplies;
#else
	return 0;
#endif
}

size_t matrix_stats_t::get_multiplies() const
{
#ifdef MATRIX_DEBUG
	return double_multiplies;
#else
	return 0;
#endif
}

void matrix_stats_t::print_diff(const matrix_stats_t &orig) const
{
#ifdef MATRIX_DEBUG
	if (this->mem_read_bytes != orig.mem_read_bytes)
		BOOST_LOG_TRIVIAL(info) << "in-mem read "
			<< (this->mem_read_bytes - orig.mem_read_bytes) << " bytes";
	if (this->mem_write_bytes != orig.mem_write_bytes)
		BOOST_LOG_TRIVIAL(info) << "in-mem write "
			<< (this->mem_write_bytes - orig.mem_write_bytes) << " bytes";
	if (this->EM_read_bytes != orig.EM_read_bytes)
		BOOST_LOG_TRIVIAL(info) << "ext-mem read "
			<< (this->EM_read_bytes - orig.EM_read_bytes) << " bytes";
	if (this->EM_write_bytes != orig.EM_write_bytes)
		BOOST_LOG_TRIVIAL(info) << "ext-mem write "
			<< (this->EM_write_bytes - orig.EM_write_bytes) << " bytes";
	if (this->double_multiplies != orig.double_multiplies)
		BOOST_LOG_TRIVIAL(info) << "multiply "
			<< (this->double_multiplies - orig.double_multiplies)
			<< " double float points";
#endif
}

matrix_stats_t matrix_stats;

}

}
