#ifndef __LOAD_BALANCER_H__
#define __LOAD_BALANCER_H__

/**
 * Copyright 2014 Da Zheng
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
#include "container.h"
#include "vertex.h"

class worker_thread;
class graph_engine;
class compute_vertex;

/**
 * This class is to help balance the load.
 * If the owner thread has finished the work originally assigned to it,
 * it can steal work from other threads through this class.
 */
class load_balancer
{
	worker_thread &owner;
	graph_engine &graph;

	// This is a local buffer that contains the completed stolen vertices.
	// All vertices here need to be returned to their owner threads.
	fifo_queue<vertex_id_t> *completed_stolen_vertices;
	int num_completed_stolen_vertices;
	// The thread where we should steal activated vertices from.
	int steal_thread_id;
public:
	load_balancer(graph_engine &_graph, worker_thread &_owner);

	~load_balancer();

	int steal_activated_vertices(compute_vertex *vertices[], int num);
	/**
	 * After the thread finishes processing the stolen vertices, it needs to
	 * return all the vertices to their owner threads.
	 */
	void return_vertices(vertex_id_t ids[], int num);

	// This method is to return all completed stolen vertices to their owner
	// threads.
	void process_completed_stolen_vertices();

	void reset();
};

#endif
