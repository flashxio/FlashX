#ifndef __CONFIG_MAP_H__
#define __CONFIG_MAP_H__

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

#include <assert.h>

#include <map>
#include <memory>
#include <string>

/*
 * This maintains key-value pairs for options.
 */
class config_map
{
	std::map<std::string, std::string> configs;

	config_map() {
	}

public:
	typedef std::shared_ptr<config_map> ptr;

	static ptr create() {
		return ptr(new config_map());
	}

	static ptr create(const std::string &conf_file);

	void add_options(const char *opts[], int num);
	// Multiple options may exist in the input string and
	// the options are separated by space.
	void add_options(const std::string &opts);

	const std::string &get_option(const std::string &name) const {
		std::map<std::string, std::string>::const_iterator it
			= configs.find(name);
		assert(it != configs.end());
		return it->second;
	}

	bool has_option(const std::string &name) const {
		return configs.find(name) != configs.end();
	}

	bool read_option(const std::string &name, std::string &value) const {
		if (has_option(name)) {
			value = get_option(name);
			return true;
		}
		else
			return false;
	}

	bool read_option_int(const std::string &name, int &value) const {
		if (has_option(name)) {
			value = atoi(get_option(name).c_str());
			return true;
		}
		else
			return false;
	}

	bool read_option_long(const std::string &name, long &value) const {
		if (has_option(name)) {
			value = atol(get_option(name).c_str());
			return true;
		}
		else
			return false;
	}

	bool read_option_bool(const std::string &name, bool &value) const {
		if (has_option(name)) {
			value = true;
			return true;
		}
		else
			return false;
	}

	const std::map<std::string, std::string> &get_options() const {
		return configs;
	}
};

#endif
