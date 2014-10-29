#ifndef __MATRIX_WORKER_THREAD_H__
#define __MATRIX_WORKER_THREAD_H__

/*
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
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

#include "thread.h"
#include "matrix_io.h"

class matrix_worker_thread: public thread
{
public:
	typedef std::shared_ptr<matrix_worker_thread> ptr;

	/*
	 */
	static ptr create(int node_id, file_io_factory::shared_ptr factory,
			matrix_io_generator::ptr gen, task_creator::ptr creator);
};

#endif
