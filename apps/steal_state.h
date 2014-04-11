#ifndef __STEAL_STATE_H__
#define __STEAL_STATE_H__

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

#include <atomic>

#include "bitmap.h"

class compute_vertex;

/**
 * This class maintains the states to handle vertices being stolen by
 * other threads.
 * There are two stages:
 *	1. normal stage: the owner thread processes messages normally.
 *	2. stealing stage: the owner thread needs to indicate it is ready
 *	for vertices to be stolen. Whenever it processes messages to a vertex,
 *	it needs to check whether the vertex has been stolen by another thread.
 *	If it is stolen, postpone processing its messages.
 */
class steal_state_t
{
	// The id of the owner thread.
	const int worker_id;
	// Indicates the vertices stolen by other threads.
	thread_safe_bitmap stolen_bitmap;

	// The number of vertices stolen from the owner thread.
	std::atomic_ulong num_stolen;
	std::atomic_ulong num_returned;

	// All threads that want to steal vertices need to increase this counter.
	std::atomic_ulong prepare_steal;
	// The owner thread indicates that it has entered the stealing stage and
	// is ready for vertices to be stolen.
	std::atomic_ulong steal_state;
	// The owner thread indicates whether it's processing messages.
	// If the value is odd, it is processing messages.
	std::atomic_ulong guard;

	graph_engine &graph;
public:
	steal_state_t(graph_engine &_graph, worker_thread &owner): worker_id(
			owner.get_worker_id()), stolen_bitmap(owner.get_num_local_vertices(),
			owner.get_node_id()), graph(_graph) {
		num_stolen = 0;
		num_returned = 0;
		prepare_steal = 0;
		steal_state = 0;
		guard = 0;
	}

	/*
	 * The two methods below are used by other worker threads.
	 */

	void steal_vertices(compute_vertex *vertices[], int num);

	void return_vertices(vertex_id_t ids[], int num);

	bool steal_mode_enabled() const {
		// It should be fine to use the relaxed memory order. steal_state
		// can only be changed in the same thread.
		return steal_state.load(std::memory_order_relaxed);
	}

	/*
	 * The methods below are used by the owner thread.
	 */

	bool is_stolen(vertex_id_t id) const {
		int part_id;
		off_t off;
		graph.get_partitioner()->map2loc(id, part_id, off);
		return stolen_bitmap.get(off);
	}

	bool is_stolen(local_vid_t id) const {
		return stolen_bitmap.get(id.id);
	}

	void guard_msg_processing() {
		guard.fetch_add(1);
		if (prepare_steal.load() > 0)
			steal_state.store(1);
	}

	void unguard_msg_processing() {
		guard.fetch_add(1);
	}

	void reset() {
		assert(num_stolen == num_returned);
		num_stolen = 0;
		num_returned = 0;
		prepare_steal = 0;
		steal_state = 0;
		guard = 0;
		stolen_bitmap.clear();
	}

	size_t get_num_stolen() const {
		return num_stolen.load(std::memory_order_relaxed);
	}

	size_t get_num_returned() const {
		return num_returned.load(std::memory_order_relaxed);
	}
};


#endif
