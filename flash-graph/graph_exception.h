#ifndef __GRAPH_EXCEPTION_H__
#define __GRAPH_EXCEPTION_H__

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

#include <string>
#include <exception>

namespace fg
{

class conf_exception: public std::exception
{
	std::string msg;
public:
	conf_exception(const std::string &msg) {
		this->msg = msg;
	}

	~conf_exception() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

class wrong_format: public std::exception
{
	std::string msg;
public:
	wrong_format(const std::string &msg) {
		this->msg = msg;
	}

	~wrong_format() throw() {
	}

	const char* what() const throw() {
		return msg.c_str();
	}
};

}

#endif
