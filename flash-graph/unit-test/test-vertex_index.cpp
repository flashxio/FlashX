#include <boost/foreach.hpp>

#include "graph_file_header.h"
#include "vertex.h"
#include "vertex_index.h"
#include "vertex_index_constructor.h"

using namespace fg;

size_t construct_large_vertex(in_mem_directed_vertex<empty_data> &v)
{
	int num_edges = compressed_vertex_entry::LARGE_VERTEX_SIZE * 2;
	for (int i = 0; i < num_edges; i++) {
		edge<empty_data> e1(i, v.get_id());
		v.add_in_edge(e1);
		edge<empty_data> e2(v.get_id(), i);
		v.add_out_edge(e2);
	}
	return num_edges;
}

size_t construct_large_vertex(in_mem_undirected_vertex<empty_data> &v)
{
	int num_edges = compressed_vertex_entry::LARGE_VERTEX_SIZE * 2;
	for (int i = 0; i < num_edges; i++) {
		edge<empty_data> e2(v.get_id(), i);
		v.add_edge(e2);
	}
	return num_edges;
}

bool operator==(const compressed_directed_vertex_entry &e1,
		const compressed_directed_vertex_entry &e2)
{
	if (e1.get_start_in_off() != e2.get_start_in_off()
			|| e1.get_start_out_off() != e2.get_start_out_off())
		return false;
	for (size_t i = 0; i < compressed_vertex_entry::ENTRY_SIZE; i++)
		if (e1.get_num_in_edges(i) != e2.get_num_in_edges(i)
				|| e1.get_num_out_edges(i) != e2.get_num_out_edges(i))
			return false;
	return true;
}

bool operator==(const compressed_undirected_vertex_entry &e1,
		const compressed_undirected_vertex_entry &e2)
{
	if (e1.get_start_off() != e2.get_start_off())
		return false;
	for (size_t i = 0; i < compressed_vertex_entry::ENTRY_SIZE; i++)
		if (e1.get_num_edges(i) != e2.get_num_edges(i))
			return false;
	return true;
}

bool operator==(const large_vertex_t &v1, const large_vertex_t &v2)
{
	return v1.first == v2.first && v1.second == v2.second;
}

template<class index_type>
void verify_vertex_index(typename index_type::ptr raw_index1,
		typename index_type::ptr raw_index2)
{
	raw_index1->verify();
	raw_index2->verify();
	assert(raw_index1->get_index_size() == raw_index2->get_index_size());
	assert(raw_index1->get_num_vertices() == raw_index2->get_num_vertices());
	assert(raw_index1->get_num_entries() == raw_index2->get_num_entries());
	assert(raw_index1->get_out_part_loc() == raw_index2->get_out_part_loc());
	assert(raw_index1->is_compressed());
	assert(raw_index2->is_compressed());

	char *raw_data1 = (char *) raw_index1.get();
	char *raw_data2 = (char *) raw_index2.get();
	size_t idx_size = raw_index1->get_index_size();
	for (size_t i = 0; i < vertex_index::get_header_size(); i++) {
		if (raw_data1[i] != raw_data2[i]) {
			printf("get different data at %ld", i);
			break;
		}
	}

	const typename index_type::entry_type *entries1 = raw_index1->get_entries();
	const typename index_type::entry_type *entries2 = raw_index2->get_entries();
	for (size_t i = 0; i < raw_index1->get_num_entries(); i++)
		assert(entries1[i] == entries2[i]);
}

void test_directed_vertex_index()
{
	std::vector<in_mem_directed_vertex<empty_data> > vertices;
	size_t num_edges = 0;
	for (int i = 0; i < 10000; i++) {
		in_mem_directed_vertex<empty_data> v(i, 0);
		if (random() % 5 == 0)
			num_edges += construct_large_vertex(v);
		vertices.push_back(v);
	}
	graph_header header(graph_type::DIRECTED, vertices.size(), num_edges, 0);
	printf("There are %ld edges\n", num_edges);

	in_mem_vertex_index::ptr cindex
		= in_mem_vertex_index::create_compressed(true, 0);
	in_mem_vertex_index::ptr index
		= in_mem_vertex_index::create(true);
	for (size_t i = 0; i < vertices.size(); i++) {
		cindex->add_vertex(vertices[i]);
		index->add_vertex(vertices[i]);
	}
	cdirected_vertex_index::ptr raw_index1
		= std::static_pointer_cast<cdirected_vertex_index>(cindex->dump(header, true));
	cdirected_vertex_index::ptr raw_index2
		= std::static_pointer_cast<cdirected_vertex_index>(index->dump(header, true));
	verify_vertex_index<cdirected_vertex_index>(raw_index1, raw_index2);

	assert(raw_index1->get_num_large_in_vertices() == raw_index2->get_num_large_in_vertices());
	assert(raw_index1->get_num_large_out_vertices() == raw_index2->get_num_large_out_vertices());
	const large_vertex_t *large_vertices1
		= ((const cdirected_vertex_index &) *raw_index1).get_large_in_vertices();
	const large_vertex_t *large_vertices2
		= ((const cdirected_vertex_index &) *raw_index2).get_large_in_vertices();
	for (size_t i = 0; i < raw_index1->get_num_large_in_vertices(); i++)
		assert(large_vertices1[i] == large_vertices2[i]);

	large_vertices1
		= ((const cdirected_vertex_index &) *raw_index1).get_large_out_vertices();
	large_vertices2
		= ((const cdirected_vertex_index &) *raw_index2).get_large_out_vertices();
	for (size_t i = 0; i < raw_index1->get_num_large_out_vertices(); i++)
		assert(large_vertices1[i] == large_vertices2[i]);
}

void test_undirected_vertex_index()
{
	std::vector<in_mem_undirected_vertex<empty_data> > vertices;
	size_t num_edges = 0;
	for (int i = 0; i < 10000; i++) {
		in_mem_undirected_vertex<empty_data> v(i, 0);
		if (random() % 5 == 0)
			num_edges += construct_large_vertex(v);
		vertices.push_back(v);
	}
	graph_header header(graph_type::DIRECTED, vertices.size(), num_edges, 0);
	printf("There are %ld edges\n", num_edges);

	in_mem_vertex_index::ptr cindex
		= in_mem_vertex_index::create_compressed(false, 0);
	in_mem_vertex_index::ptr index
		= in_mem_vertex_index::create(false);
	for (size_t i = 0; i < vertices.size(); i++) {
		cindex->add_vertex(vertices[i]);
		index->add_vertex(vertices[i]);
	}
	cundirected_vertex_index::ptr raw_index1
		= std::static_pointer_cast<cundirected_vertex_index>(cindex->dump(header, true));
	cundirected_vertex_index::ptr raw_index2
		= std::static_pointer_cast<cundirected_vertex_index>(index->dump(header, true));
	verify_vertex_index<cundirected_vertex_index>(raw_index1, raw_index2);

	assert(raw_index1->get_num_large_vertices() == raw_index2->get_num_large_vertices());
	const large_vertex_t *large_vertices1
		= ((const cundirected_vertex_index &) *raw_index1).get_large_vertices();
	const large_vertex_t *large_vertices2
		= ((const cundirected_vertex_index &) *raw_index2).get_large_vertices();
	for (size_t i = 0; i < raw_index1->get_num_large_vertices(); i++)
		assert(large_vertices1[i] == large_vertices2[i]);
}

int main()
{
	test_directed_vertex_index();
	test_undirected_vertex_index();
}
