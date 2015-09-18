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

#include <algorithm>
#include <iostream>
#include <sstream>

#include "common.h"
#include "config_map.h"

static bool read_config_file(const std::string &conf_file,
		std::map<std::string, std::string> &configs)
{
	FILE *f = fopen(conf_file.c_str(), "r");
	if (f == NULL) {
		perror("fopen");
		return false;
	}

	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	while ((read = getline(&line, &len, f)) > 0) {
		std::string str = line;
		if (str.length() == 1)
			continue;
		if (line[0] == '#')
			continue;

		size_t found = str.find("=");
		/* if there isn't `=', I assume it's a file name*/
		if (found == std::string::npos) {
			BOOST_LOG_TRIVIAL(error) << std::string("wrong format: ") + line;
			return false;
		}

		std::string value = str.substr(found + 1);
		value.erase(std::remove_if(value.begin(), value.end(), isspace),
				value.end());
		std::string key = str.substr(0, found);
		key.erase(std::remove_if(key.begin(), key.end(), isspace),
				key.end());
		bool res = configs.insert(std::pair<std::string, std::string>(key,
					value)).second;
		if (!res)
			configs[key] = value;
	}
	fclose(f);
	return true;
}

config_map::ptr config_map::create(const std::string &conf_file)
{
	config_map::ptr map = config_map::create();
	// If we can't read the config file, return an empty pointer.
	if (!read_config_file(conf_file, map->configs))
		return config_map::ptr(NULL);
	else
		return map;
}

/*
 * All options should have the following format:
 *		key=value.
 * All options that don't have the format are ignored.
 */
void config_map::add_options(const char *argv[], int argc)
{
	for (int i = 0; i < argc; i++) {
		std::string str = argv[i];

		size_t found = str.find("=");
		if (found == std::string::npos) {
			continue;
		}

		std::string value = str.substr(found + 1);
		value.erase(std::remove_if(value.begin(), value.end(), isspace),
				value.end());
		std::string key = str.substr(0, found);
		key.erase(std::remove_if(key.begin(), key.end(), isspace),
				key.end());
		bool res = configs.insert(std::pair<std::string, std::string>(key,
					value)).second;
		if (!res)
			configs[key] = value;
	}
}

size_t split_string(const std::string &str, const std::string &delimiter,
		std::vector<std::string> &splits)
{
	std::stringstream ss(str);
	size_t num = 0;
	while (ss.good()) {
		std::string s;
		ss >> s;
		splits.push_back(s);
		num++;
	}
	return num;
}

void config_map::add_options(const std::string &opts)
{
	std::vector<std::string> parts;
	split_string(opts, " ", parts);
	std::vector<const char *> part_chars;
	for (size_t i = 0; i < parts.size(); i++) {
		part_chars.push_back(parts[i].c_str());
	}
	config_map::add_options(part_chars.data(), part_chars.size());
}
