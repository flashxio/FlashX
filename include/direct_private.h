#ifndef __DIRECT_PRIVATE_H__
#define __DIRECT_PRIVATE_H__

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

#include "read_private.h"

class direct_io: public buffered_io
{
public:
	direct_io(const logical_file_partition &partition,
			thread *t): buffered_io(partition, t, O_DIRECT | O_RDWR) {
	}

	io_status access(char *buf, off_t offset, ssize_t size, int access_method);
};

#endif
