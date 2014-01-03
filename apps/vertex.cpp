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

#include "vertex.h"

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
