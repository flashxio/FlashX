#ifndef __MESSAGE_PROCESSOR_H__
#define __MESSAGE_PROCESSOR_H__

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

#include <memory>

#include "container.h"

#include "messaging.h"
#include "vertex.h"

namespace fg
{

class graph_engine;
class worker_thread;
class steal_state_t;
class compute_vertex;
class compute_vertex_pointer;

/**
 * This class is to process the messages sent to the owner thread.
 * The complexity here is that some messages can't be processed when
 * we try to process them. Therefore, we need an extra buffer to keep them
 * and process them later.
 */
class message_processor
{
	graph_engine &graph;
	worker_thread &owner;

	std::shared_ptr<slab_allocator> msg_alloc;

	// The queue of messages sent from other threads.
	msg_queue msg_q;

	// The thread state of handling stealing vertices.
	std::unique_ptr<steal_state_t> steal_state;

	// This is a message buffer to keep all messages whose destination vertices
	// have been stolen by other threads.
	fifo_queue<message> stolenv_msgs;

	void buf_msg(vertex_message &msg);
	void buf_mmsg(local_vid_t id, multicast_message &mmsg);

	void process_msg(message &msg, bool check_steal);
	void process_multicast_msg(multicast_message &mmsg, bool check_steal);

public:
	message_processor(graph_engine &_graph, worker_thread &_owner,
			std::shared_ptr<slab_allocator> msg_alloc);

	void process_msgs();

	void steal_vertices(compute_vertex_pointer vertices[], int num);
	void return_vertices(vertex_id_t ids[], int num);

	msg_queue &get_msg_queue() {
		return msg_q;
	}

	void reset();
};

}

#endif
