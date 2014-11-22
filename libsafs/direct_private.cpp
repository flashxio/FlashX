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

#include "direct_private.h"

namespace safs
{

io_status direct_io::access(char *buf, off_t offset,
		ssize_t size, int access_method) {
	assert(size >= MIN_BLOCK_SIZE);
	assert(size % MIN_BLOCK_SIZE == 0);
	assert(offset % MIN_BLOCK_SIZE == 0);
	assert((long) buf % MIN_BLOCK_SIZE == 0);
	return buffered_io::access(buf, offset, size, access_method);
}

}
