#include "vertex.h"

int ext_mem_directed_vertex::serialize(const in_mem_directed_vertex &in_v,
		char *buf, int size)
{
	int mem_size = in_v.get_serialize_size();
	assert(size >= mem_size);
	ext_mem_directed_vertex *ext_v = (ext_mem_directed_vertex *) buf;
	ext_v->id = in_v.get_id();
	ext_v->num_in_edges = in_v.get_num_in_edges();
	ext_v->num_out_edges = in_v.get_num_out_edges();

	vertex_id_t *neighbors = ext_v->neighbors;
	for (int i = 0; i < in_v.get_num_in_edges(); i++) {
		neighbors[i] = in_v.get_in_edge(i).get_from();
	}
	for (int i = 0; i < in_v.get_num_out_edges(); i++) {
		neighbors[i + in_v.get_num_in_edges()] = in_v.get_out_edge(i).get_to();
	}
	return mem_size;
}

int ext_mem_undirected_vertex::serialize(const in_mem_undirected_vertex &v,
		char *buf, int size)
{
	int mem_size = v.get_serialize_size();
	assert(size >= mem_size);
	ext_mem_undirected_vertex *ext_v = (ext_mem_undirected_vertex *) buf;
	ext_v->id = v.get_id();
	ext_v->num_edges = v.get_num_edges();

	vertex_id_t *neighbors = ext_v->neighbors;
	for (int i = 0; i < v.get_num_edges(); i++) {
		neighbors[i] = v.get_edge(i).get_to();
	}

	return mem_size;
}

bool ext_mem_directed_vertex::is_edge_list_sorted(edge_type type) const
{
	if (type == edge_type::IN_EDGE || type == edge_type::BOTH_EDGES) {
		int num_in_edges = get_num_edges(edge_type::IN_EDGE);
		if (num_in_edges > 0) {
			vertex_id_t id = this->get_neighbor(edge_type::IN_EDGE, 0);
			for (int i = 1; i < num_in_edges; i++) {
				vertex_id_t i_neighbor = this->get_neighbor(edge_type::IN_EDGE, i);
				if (i_neighbor < id)
					return false;
				id = i_neighbor;
			}
		}
	}
	if (type == edge_type::OUT_EDGE || type == edge_type::BOTH_EDGES) {
		int num_out_edges = get_num_edges(edge_type::OUT_EDGE);
		if (num_out_edges > 0) {
			vertex_id_t id = this->get_neighbor(edge_type::OUT_EDGE, 0);
			for (int i = 1; i < num_out_edges; i++) {
				vertex_id_t i_neighbor = this->get_neighbor(edge_type::OUT_EDGE, i);
				if (i_neighbor < id)
					return false;
				id = i_neighbor;
			}
		}
	}

	return true;
}

bool ext_mem_undirected_vertex::is_edge_list_sorted(edge_type type) const
{
	int num_in_edges = get_num_edges(type);
	if (num_in_edges > 0) {
		vertex_id_t id = this->get_neighbor(type, 0);
		for (int i = 1; i < num_in_edges; i++) {
			vertex_id_t i_neighbor = this->get_neighbor(type, i);
			if (i_neighbor < id)
				return false;
			id = i_neighbor;
		}
	}

	return true;
}
