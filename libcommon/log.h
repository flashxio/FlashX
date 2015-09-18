#ifndef __MY_LOG_H__
#define __MY_LOG_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of libcommon.
 * This library is used by SAFS and FlashGraph, and potentially by other
 * systems.
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

#include <iostream>

enum c_log_level
{
	debug,
	info,
	warning,
	error,
	fatal,
};

#ifdef USE_BOOST_LOG

#include <boost/log/trivial.hpp>

static inline void set_log_level(enum c_log_level level)
{
	boost::log::core::get()->set_filter(
			boost::log::trivial::severity > boost::log::trivial::info);
}

#else

class simple_log_stream
{
	static c_log_level curr_log_level;
	c_log_level level;
	bool log_data;
public:
	simple_log_stream(c_log_level level) {
		this->level = level;
		log_data = false;
	}

	~simple_log_stream() {
		if (!log_data && level >= curr_log_level)
			std::cout << "\n";
	}

	static void set_global_log_level(c_log_level level) {
		curr_log_level = level;
	}

	template<class T>
	simple_log_stream operator<<(const T &v) {
		log_data = true;
		if (level >= curr_log_level)
			std::cout << v;
		return simple_log_stream(level);
	}
};

extern simple_log_stream log_debug;
extern simple_log_stream log_info;
extern simple_log_stream log_warning;
extern simple_log_stream log_error;
extern simple_log_stream log_fatal;

#define BOOST_LOG_TRIVIAL(x)		\
	log_##x

static inline void set_log_level(enum c_log_level level)
{
	simple_log_stream::set_global_log_level(level);
}

#endif

#endif
