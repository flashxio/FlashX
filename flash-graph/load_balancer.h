#ifndef __LOAD_BALANCER_H__
#define __LOAD_BALANCER_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "container.h"
#include "vertex.h"

class worker_thread;
class graph_engine;
class compute_vertex;
class compute_vertex_pointer;

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

	int steal_activated_vertices(compute_vertex_pointer vertices[], int num);
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
