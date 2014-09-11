#include "vertex_index.h"

int main(int argc, char *argv[])
{
	if (argc < 3) {
		fprintf(stderr,
				"compress_graph_index index_file compressed_index_file\n");
		exit(1);
	}

	std::string index_file = argv[1];
	std::string compressed_index_file = argv[2];

	vertex_index::ptr index = vertex_index::load(index_file);
	assert(index->get_graph_header().is_directed_graph());
	cdirected_vertex_index::ptr cindex = cdirected_vertex_index::construct(
			(directed_vertex_index &) *index);
	cindex->dump(compressed_index_file);

	vertex_index::ptr index1 = vertex_index::load(compressed_index_file);
	in_mem_cdirected_vertex_index::ptr in_mem_cindex
		= in_mem_cdirected_vertex_index::create(*index1);
	in_mem_cindex->verify_against(*directed_vertex_index::cast(index));
}
