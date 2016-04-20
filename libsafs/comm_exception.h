#ifndef __COMM_EXCEPTION_H__
#define __COMM_EXCEPTION_H__

/*
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

/* out of memory exception */
class oom_exception: public std::exception
{
	std::string msg;
public:
	oom_exception() {
	}

	oom_exception(const std::string &msg) {
		this->msg = msg;
	}

	~oom_exception() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

class unsupported_exception: public std::exception
{
	std::string msg;
public:
	unsupported_exception() {
	}

	unsupported_exception(const std::string &msg) {
		this->msg = msg;
	}

	~unsupported_exception() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

class invalid_arg_exception: public std::exception
{
	std::string msg;
public:
	invalid_arg_exception(const std::string &msg) {
		this->msg = msg;
	}

	~invalid_arg_exception() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

#endif
