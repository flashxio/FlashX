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

#include "thread.h"

#include "vertex_program.h"
#include "messaging.h"
#include "worker_thread.h"
#include "message_processor.h"

vertex_program::~vertex_program()
{
	for (unsigned i = 0; i < msg_senders.size(); i++)
		simple_msg_sender::destroy(msg_senders[i]);
	for (unsigned i = 0; i < multicast_senders.size(); i++)
		multicast_msg_sender::destroy(multicast_senders[i]);
	for (unsigned i = 0; i < activate_senders.size(); i++)
		multicast_msg_sender::destroy(activate_senders[i]);
}

void vertex_program::init_messaging(const std::vector<worker_thread *> &threads,
		std::shared_ptr<slab_allocator> msg_alloc)
{
	vid_bufs = std::unique_ptr<std::vector<local_vid_t>[]>(
			new std::vector<local_vid_t>[graph->get_num_threads()]);
	vloc_size = graph->get_num_threads() * 2;
	vertex_locs = std::unique_ptr<vertex_loc_t[]>(new vertex_loc_t[vloc_size]);

	for (unsigned i = 0; i < threads.size(); i++) {
		msg_senders.push_back(simple_msg_sender::create(t->get_node_id(),
					msg_alloc, &threads[i]->get_msg_processor().get_msg_queue()));
		multicast_senders.push_back(multicast_msg_sender::create(msg_alloc,
					&threads[i]->get_msg_processor().get_msg_queue()));
		multicast_msg_sender *activate_sender = multicast_msg_sender::create(
				msg_alloc, &threads[i]->get_msg_processor().get_msg_queue());
		// All activation messages are the same. We can initialize the sender
		// here.
		activation_message msg;
		activate_sender->init(msg);
		activate_senders.push_back(activate_sender);
	}
}

void vertex_program::multicast_msg(vertex_id_t ids[], int num,
		vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(t == curr);

	if (num == 0)
		return;

	if (num < graph->get_num_threads() * 2) {
		for (int i = 0; i < num; i++)
			this->send_msg(ids[i], msg);
		return;
	}

	graph->get_partitioner()->map2loc(ids, num, vid_bufs.get(),
			graph->get_num_threads());
	for (int i = 0; i < graph->get_num_threads(); i++) {
		if (vid_bufs[i].empty())
			continue;

		multicast_msg_sender &sender = get_multicast_sender(i);
		sender.init(msg);
		int ret = sender.add_dests(vid_bufs[i].data(), vid_bufs[i].size());
		assert((size_t) ret == vid_bufs[i].size());
		vid_bufs[i].clear();
		sender.end_multicast();
	}
}

void vertex_program::multicast_msg(edge_seq_iterator &it, vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);

	int num_dests = it.get_num_tot_entries();
	if (num_dests == 0)
		return;

	if (num_dests < graph->get_num_threads() * 2) {
		PAGE_FOREACH(vertex_id_t, id, it) {
			this->send_msg(id, msg);
		} PAGE_FOREACH_END
		return;
	}

	graph->get_partitioner()->map2loc(it, vid_bufs.get(),
			graph->get_num_threads());
	for (int i = 0; i < graph->get_num_threads(); i++) {
		if (vid_bufs[i].empty())
			continue;

		multicast_msg_sender &sender = get_multicast_sender(i);
		sender.init(msg);
		int ret = sender.add_dests(vid_bufs[i].data(), vid_bufs[i].size());
		assert((size_t) ret == vid_bufs[i].size());
		vid_bufs[i].clear();
		sender.end_multicast();
	}
}

void vertex_program::send_msg(vertex_id_t dest, vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(t == curr);

	int part_id;
	// We are going to use the offset of a vertex in a partition as
	// the ID of the vertex.
	off_t local_id;
	graph->get_partitioner()->map2loc(dest, part_id, local_id);
	simple_msg_sender &sender = get_msg_sender(part_id);
	msg.set_dest(local_vid_t(local_id));
	sender.send_cached(msg);
}

void vertex_program::activate_vertices(vertex_id_t ids[], int num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);

	if (num == 0)
		return;

	// When there are very few destinations, this way is cheaper.
	if ((size_t) num <= vloc_size) {
		for (int i = 0; i < num; i++) {
			int part_id;
			// We are going to use the offset of a vertex in a partition as
			// the ID of the vertex.
			off_t local_id;
			graph->get_partitioner()->map2loc(ids[i], part_id, local_id);
			multicast_msg_sender &sender = get_activate_sender(part_id);
			bool ret = sender.add_dest((local_vid_t) local_id);
			assert(ret);
		}
		return;
	}

	graph->get_partitioner()->map2loc(ids, num, vid_bufs.get(),
			graph->get_num_threads());
	for (int i = 0; i < graph->get_num_threads(); i++) {
		multicast_msg_sender &sender = get_activate_sender(i);
		if (vid_bufs[i].empty())
			continue;
		int ret = sender.add_dests(vid_bufs[i].data(), vid_bufs[i].size());
		assert((size_t) ret == vid_bufs[i].size());
		vid_bufs[i].clear();
	}
}

void vertex_program::activate_vertices(edge_seq_iterator &it)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);

	size_t num_dests = it.get_num_tot_entries();
	if (num_dests == 0)
		return;

	// When there are very few destinations, this way is cheaper.
	if (num_dests <= vloc_size) {
		size_t ret = graph->get_partitioner()->map2loc(it, vertex_locs.get(),
				vloc_size);
		assert(ret == num_dests);
		for (size_t i = 0; i < ret; i++) {
			int part_id = vertex_locs[i].first;
			// We are going to use the offset of a vertex in a partition as
			// the ID of the vertex.
			local_vid_t local_id = vertex_locs[i].second;
			multicast_msg_sender &sender = get_activate_sender(part_id);
			bool ret = sender.add_dest(local_id);
			assert(ret);
		}
		return;
	}

	graph->get_partitioner()->map2loc(it, vid_bufs.get(),
			graph->get_num_threads());
	for (int i = 0; i < graph->get_num_threads(); i++) {
		multicast_msg_sender &sender = get_activate_sender(i);
		if (vid_bufs[i].empty())
			continue;
		int ret = sender.add_dests(vid_bufs[i].data(), vid_bufs[i].size());
		assert((size_t) ret == vid_bufs[i].size());
		vid_bufs[i].clear();
	}
}

void vertex_program::flush_msgs()
{
	for (size_t i = 0; i < msg_senders.size(); i++)
		msg_senders[i]->flush();
	for (size_t i = 0; i < multicast_senders.size(); i++)
		multicast_senders[i]->flush();
	for (size_t i = 0; i < activate_senders.size(); i++) {
		activate_senders[i]->flush();
		activation_message msg;
		activate_senders[i]->init(msg);
	}
}
