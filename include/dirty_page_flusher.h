#ifndef __DIRTY_PAGE_FLUSHER_H__
#define __DIRTY_PAGE_FLUSHER_H__

/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "common.h"

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

#endif
