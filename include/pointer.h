#ifndef __POINTER_H__
#define __POINTER_H__

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
