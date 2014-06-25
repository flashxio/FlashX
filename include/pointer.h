#ifndef __POINTER_H__
#define __POINTER_H__

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

class ptr_interface
{
	int num_refs;
public:
	ptr_interface() {
		num_refs = 0;
	}

	void inc_ref() {
		num_refs++;
	}

	void dec_ref() {
		num_refs--;
	}

	int get_ref() const {
		return num_refs;
	}

	void wait4unref() {
		while (num_refs > 0) {}
	}
};

/**
 * This is a thread-safe version of ptr_interface.
 */
class TS_ptr_interface
{
	atomic_integer num_refs;
public:
	void inc_ref() {
		num_refs.inc(1);
	}

	void dec_ref() {
		num_refs.dec(1);
	}

	int get_ref() const {
		return num_refs.get();
	}

	void wait4unref() {
		while (get_ref() > 0) {}
	}
};

#endif
