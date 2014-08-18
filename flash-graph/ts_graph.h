#ifndef __TS_ALGS_H__
#define __TS_ALGS_H__

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

#include <boost/date_time/posix_time/posix_time.hpp>

#include "vertex.h"

const int HOUR_SECS = 3600;
const int DAY_SECS = HOUR_SECS * 24;
const int MONTH_SECS = DAY_SECS * 30;

page_byte_array::seq_const_iterator<vertex_id_t> get_ts_iterator(
		const page_directed_vertex &v, edge_type type, time_t time_start,
		time_t time_interval);

static inline bool is_time_str(const std::string &str)
{
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	char *ret = strptime(str.c_str(), "%Y-%m-%d", &tm);
	return ret != NULL;
}

static inline time_t conv_str_to_time(const std::string &str)
{
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	char *ret = strptime(str.c_str(), "%Y-%m-%d", &tm);
	assert(ret);
	tm.tm_isdst = 1;
	return mktime(&tm);
}

#endif
