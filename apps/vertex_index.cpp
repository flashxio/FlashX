/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "common.h"
#include "native_file.h"

#include "vertex_index.h"

vertex_index *vertex_index::load(const std::string &index_file)
{
	native_file local_f(index_file);
	ssize_t size = local_f.get_size();
	assert((unsigned) size >= sizeof(vertex_index));
	char *buf = (char *) malloc(size);
	FILE *fd = fopen(index_file.c_str(), "r");
	size_t ret = fread(buf, size, 1, fd);
	assert(ret == 1);

	vertex_index *idx = (vertex_index *) buf;
	assert((unsigned) size >= sizeof(vertex_index)
			+ idx->get_num_vertices() * sizeof(idx->vertex_offs[0]));
	idx->header.verify();

	return idx;
}

void vertex_index::dump(const std::string &file)
{
	FILE *f = fopen(file.c_str(), "w");
	if (f == NULL) {
		perror("fopen");
		assert(0);
	}

	ssize_t ret = fwrite(this, get_serialize_size(), 1, f);
	assert(ret);

	fclose(f);
}
