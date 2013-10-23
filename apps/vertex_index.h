#ifndef __VERTEX_INDEX_H__
#define __VERTEX_INDEX_H__

#include <string.h>

#include <string>

#include "vertex.h"

/**
 * This vertex index maps a vertex id to the location of the vertex in a file.
 */
class vertex_index
{
	size_t num_vertices;
	// The total size of the graph in the form of adjacency list
	// in the external memory.
	size_t tot_size;
	off_t vertex_offs[0];

	vertex_index(size_t num) {
		num_vertices = 0;
		tot_size = 0;
		memset(vertex_offs, 0, sizeof(vertex_offs[0]) * num);
	}

	size_t get_serialize_size() const {
		return sizeof(vertex_index) + num_vertices * sizeof(vertex_offs[0]);
	}
public:
	static vertex_index *load(const std::string &index_file);
	static void destroy(vertex_index *index) {
		free(index);
	}
	template<class vertex_type>
	static vertex_index *create(const std::vector<vertex_type> &vertices) {
		void *addr = malloc(sizeof(vertex_index)
				+ sizeof(off_t) * vertices.size());
		vertex_index *index = new (addr) vertex_index(vertices.size());
		index->num_vertices = vertices.size();
		size_t tot_size = 0;
		for (size_t i = 0; i < vertices.size(); i++) {
			index->vertex_offs[i] = tot_size;
			tot_size += vertices[i].get_serialize_size();
		}
		index->tot_size = tot_size;
		return index;
	}

	void dump(const std::string &file);

	off_t get_vertex_off(vertex_id_t id) const {
		assert(id < num_vertices);
		return vertex_offs[id];
	}

	int get_vertex_size(vertex_id_t id) const {
		assert(id < num_vertices);
		if (id < num_vertices - 1)
			return vertex_offs[id + 1] - vertex_offs[id];
		else
			return tot_size - vertex_offs[id];
	}

	size_t get_num_vertices() const {
		return num_vertices;
	}
};

#endif
