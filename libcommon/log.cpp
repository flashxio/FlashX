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

#include "log.h"

c_log_level simple_log_stream::curr_log_level;

simple_log_stream log_debug(c_log_level::debug);
simple_log_stream log_info(c_log_level::info);
simple_log_stream log_warning(c_log_level::warning);
simple_log_stream log_error(c_log_level::error);
simple_log_stream log_fatal(c_log_level::fatal);
