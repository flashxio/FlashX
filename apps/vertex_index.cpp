#include "common.h"

#include "vertex_index.h"

vertex_index *vertex_index::load(const std::string &index_file)
{
	ssize_t size = get_file_size(index_file.c_str());
	assert((unsigned) size >= sizeof(vertex_index));
	char *buf = (char *) malloc(size);
	FILE *fd = fopen(index_file.c_str(), "r");
	size_t ret = fread(buf, size, 1, fd);
	assert(ret == 1);

	vertex_index *idx = (vertex_index *) buf;
	assert((unsigned) size >= sizeof(vertex_index)
			+ idx->num_vertices * sizeof(idx->vertex_offs[0]));
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
