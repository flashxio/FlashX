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

#include "file_mapper.h"

int gen_RAID_rand_start(int num_files)
{
	int r = random();
	while (r == 0) {
		r = random();
	}
	return r % num_files;
}

int RAID0_mapper::rand_start;
int RAID5_mapper::rand_start;

atomic_integer file_mapper::file_id_gen;
