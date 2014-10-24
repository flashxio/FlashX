#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of SAFSlib.
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
#include <string>
#include <exception>

class no_space_exception
{
};

/* out of memory exception */
class oom_exception
{
};

class unsupported_exception
{
};

class init_error: public std::exception
{
	std::string msg;
public:
	init_error(const std::string &msg) {
		this->msg = msg;
	}

	~init_error() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

#endif
