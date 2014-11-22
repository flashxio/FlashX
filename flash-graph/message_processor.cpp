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

#include "message_processor.h"
#include "messaging.h"
#include "graph_engine.h"
#include "worker_thread.h"
#include "steal_state.h"

namespace fg
{

message_processor::message_processor(graph_engine &_graph,
		worker_thread &_owner, std::shared_ptr<slab_allocator> msg_alloc): graph(_graph),
	owner(_owner), msg_q(_owner.get_node_id(), "graph_msg_queue", 16, INT_MAX),
	stolenv_msgs(_owner.get_node_id(), PAGE_SIZE, true)
{
	if (graph_conf.use_serial_run())
		steal_state = std::unique_ptr<steal_state_t>(new steal_state_t(graph, owner));
	this->msg_alloc = msg_alloc;
}

void message_processor::buf_msg(vertex_message &vmsg)
{
	if (stolenv_msgs.is_empty() || stolenv_msgs.back().add(vmsg) == NULL) {
		message msg(msg_alloc.get());
		if (stolenv_msgs.is_full())
			stolenv_msgs.expand_queue(stolenv_msgs.get_size() * 2);
		stolenv_msgs.push_back(msg);
		// We have to make sure this is successful.
		BOOST_VERIFY(stolenv_msgs.back().add(vmsg));
	}
}

/**
 * This class converts a multicast message to a p2p message.
 */
class multicast_p2p_converter
{
	local_vid_t dest;
	multicast_message &mmsg;
public:
	multicast_p2p_converter(local_vid_t _dest,
			multicast_message &_mmsg): dest(_dest), mmsg(_mmsg) {
	}

	int get_serialized_size() const {
		return mmsg.get_body_size();
	}

	int serialize(char *buf, int size) const {
		int serialized_size = this->get_serialized_size();
		assert(serialized_size <= size);
		memcpy(buf, &mmsg, serialized_size);
		vertex_message *vmsg = (vertex_message *) buf;
		*vmsg = vertex_message(serialized_size, mmsg.is_activate());
		vmsg->set_dest(dest);
		return serialized_size;
	}
};

void message_processor::buf_mmsg(local_vid_t id, multicast_message &mmsg)
{
	multicast_p2p_converter converter(id, mmsg);
	if (stolenv_msgs.is_empty() || stolenv_msgs.back().add(converter) == NULL) {
		message msg(msg_alloc.get());
		if (stolenv_msgs.is_full())
			stolenv_msgs.expand_queue(stolenv_msgs.get_size() * 2);
		stolenv_msgs.push_back(msg);
		// We have to make sure this is successful.
		BOOST_VERIFY(stolenv_msgs.back().add(converter));
	}
}

void message_processor::process_multicast_msg(multicast_message &mmsg,
		bool check_steal)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	// Messages are always processed in the main vertex.
	vertex_program &curr_vprog = t->get_vertex_program(false);

	int num_dests = mmsg.get_num_dests();
	multicast_dest_list dest_list = mmsg.get_dest_list();

	if (mmsg.is_activation_msg()) {
		owner.activate_vertices(dest_list.get_dests(), num_dests);
		return;
	}

	if (!check_steal) {
		curr_vprog.run_on_multicast_message(mmsg);
		if (mmsg.is_activate())
			owner.activate_vertices(dest_list.get_dests(), num_dests);
		return;
	}

	for (int i = 0; i < num_dests; i++) {
		local_vid_t id = dest_list.get_dest(i);
		if (check_steal && steal_state && steal_state->is_stolen(id)) {
			buf_mmsg(id, mmsg);
		}
		else {
			compute_vertex &info = graph.get_vertex(owner.get_worker_id(), id);
			curr_vprog.run_on_message(info, mmsg);
		}
		if (mmsg.is_activate())
			owner.activate_vertex(id);
	}
}

void message_processor::process_msg(message &msg, bool check_steal)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	// Messages are always processed in the main vertex.
	vertex_program &curr_vprog = t->get_vertex_program(false);

	const int VMSG_BUF_SIZE = 128;
	vertex_message *v_msgs[VMSG_BUF_SIZE];
	while (!msg.is_empty()) {
		int num = msg.get_next(v_msgs, VMSG_BUF_SIZE);
		assert(num > 0);
		// If we aren't in the mode of load balancing and these aren't
		// multicast messages, we can use the fast path.
		// We only need to check the first message. All messages are
		// of the same type.
		if (!check_steal && !v_msgs[0]->is_multicast()) {
			curr_vprog.run_on_messages((const vertex_message **) v_msgs, num);
			for (int i = 0; i < num; i++) {
				local_vid_t id = v_msgs[i]->get_dest();
				if (v_msgs[i]->is_activate())
					owner.activate_vertex(id);
			}
			continue;
		}

		if (v_msgs[0]->is_multicast()) {
			for (int i = 0; i < num; i++)
				process_multicast_msg(*multicast_message::cast2multicast(
							v_msgs[i]), check_steal);
			continue;
		}

		// If we are here, we are in the load balancing mode.
		assert(check_steal);
		for (int i = 0; i < num; i++) {
			local_vid_t id = v_msgs[i]->get_dest();
			if (steal_state && steal_state->is_stolen(id)) {
				buf_msg(*v_msgs[i]);
			}
			else {
				compute_vertex &info = graph.get_vertex(owner.get_worker_id(), id);
				curr_vprog.run_on_message(info, *v_msgs[i]);
			}
			if (v_msgs[i]->is_activate())
				owner.activate_vertex(id);
		}
	}
}

void message_processor::process_msgs()
{
	if (steal_state && steal_state->get_num_returned() > 0 && !stolenv_msgs.is_empty()) {
		// TODO we might have to make sure that a lot of messages can be
		// processed. Otherwise, we are wasting time.
		// The vertices have been processed so no other threads will steal
		// them and process them.
		// TODO will it be still true if we allow a vertex to be activated
		// multiple times in the same iteration?
		stack_array<message> msgs(stolenv_msgs.get_num_entries());
		int num_fetched = stolenv_msgs.fetch(msgs.data(),
				stolenv_msgs.get_num_entries());
		for (int i = 0; i < num_fetched; i++)
			process_msg(msgs[i], true);
	}

	const int MSG_BUF_SIZE = 16;
	message msgs[MSG_BUF_SIZE];
	bool check_steal = false;
	if (steal_state) {
		steal_state->guard_msg_processing();
		check_steal = steal_state->steal_mode_enabled();
	}
	while (!msg_q.is_empty()) {
		int num_fetched = msg_q.fetch(msgs, MSG_BUF_SIZE);
		for (int i = 0; i < num_fetched; i++)
			process_msg(msgs[i], check_steal);
	}
	if (steal_state)
		steal_state->unguard_msg_processing();
}

void message_processor::steal_vertices(compute_vertex_pointer vertices[], int num)
{
	if (steal_state)
		steal_state->steal_vertices(vertices, num);
}

void message_processor::reset()
{
	if (steal_state)
		steal_state->reset();
	assert(msg_q.is_empty());
	assert(stolenv_msgs.is_empty());
}

void message_processor::return_vertices(vertex_id_t ids[], int num)
{
	if (steal_state)
		steal_state->return_vertices(ids, num);
}

void steal_state_t::steal_vertices(compute_vertex_pointer vertices[], int num)
{
	prepare_steal.fetch_add(1);
	int num_locals = 0;
	for (int i = 0; i < num; i++) {
		// If the vertex is vertically partitioned, the vertex state doesn't
		// need to process messages anyway. The messages are processed by
		// the main vertex.
		if (vertices[i].is_part())
			continue;

		// TODO we can avoid redudant computation here.
		vertex_id_t id = graph.get_graph_index().get_vertex_id(worker_id,
				*vertices[i].get());
		// If the stolen vertex doesn't belong to the worker thread, we'll get
		// INVALID_VERTEX_ID. This shouldn't happen.
		assert(id != INVALID_VERTEX_ID);
		int part_id;
		off_t off;
		graph.get_partitioner()->map2loc(id, part_id, off);
		stolen_bitmap.set(off);
		num_locals++;
	}
	// TODO I need a memory barrier here.
	// If the guard is odd, it means the owner thread is processing
	// messages. Wait for it to finish processing messages.
	// The thread can't proceed regardless of the state of the owner
	// thread if the owner thread is processing messages.
	while (guard.load() % 2 > 0);
	num_stolen += num_locals;
}

void steal_state_t::return_vertices(vertex_id_t ids[], int num)
{
	num_returned += num;
	for (int i = 0; i < num; i++) {
		int part_id;
		off_t off;
		graph.get_partitioner()->map2loc(ids[i], part_id, off);
		// The vertices have to be returned to the owner thread.
		assert(worker_id == part_id);
		stolen_bitmap.clear(off);
	}
}

}
