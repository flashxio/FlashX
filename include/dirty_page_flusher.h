#ifndef __DIRTY_PAGE_FLUSHER_H__
#define __DIRTY_PAGE_FLUSHER_H__

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

#include "common.h"

namespace safs
{

class thread_safe_page;
class io_interface;

class dirty_page_flusher
{
public:
	virtual ~dirty_page_flusher() {
	}

	/**
	 * This method flushes dirty pages with the specified io instance.
	 * The method has to be thread-safe, as it will be invoked in different
	 * thread contexts. It also shouldn't block threads.
	 */
	virtual void flush_dirty_pages(thread_safe_page *pages[], int num,
			io_interface *io) = 0;

	virtual int flush_dirty_pages(page_filter *filter, int max_num) = 0;
};

}

#endif
