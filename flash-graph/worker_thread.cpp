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

#include <atomic>

#include "worker_thread.h"
#include "graph_engine.h"
#include "bitmap.h"
#include "vertex_compute.h"
#include "message_processor.h"
#include "load_balancer.h"
#include "steal_state.h"
#include "vertex_index_reader.h"

static void delete_val(std::vector<vertex_id_t> &vec, vertex_id_t val)
{
	size_t curr = 0;
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec[i] != val) {
			vec[curr] = vec[i];
			curr++;
		}
	}
	vec.resize(curr);
}

/*
 * This method split a list of vertices into a list of vertically
 * partitioned vertices and a list of unpartitioned vertices.
 * The input vertex list is sorted on vertex ID and will have
 * the unpartitioned vertices.
 */
static void split_vertices(graph_index &index, int part_id,
		std::vector<vertex_id_t> &vertices,
		std::vector<vpart_vertex_pointer> &vpart_ps)
{
	assert(std::is_sorted(vertices.begin(), vertices.end()));
	// Get the vertically partitioned vertices that are activated.
	std::vector<vpart_vertex_pointer> vpart_ps_tmp(
			index.get_num_vpart_vertices(part_id));
	index.get_vpart_vertex_pointers(part_id, vpart_ps_tmp.data(),
			vpart_ps_tmp.size());
	for (size_t i = 0, j = 0; i < vpart_ps_tmp.size() && j < vertices.size();) {
		vpart_vertex_pointer p = vpart_ps_tmp[i];
		vertex_id_t id = vertices[j];
		if (p.get_vertex_id() == id) {
			i++;
			vertices[j] = INVALID_VERTEX_ID;
			j++;
			vpart_ps.push_back(p);
		}
		else if (p.get_vertex_id() > id) {
			j++;
		}
		else {
			i++;
		}
	}
	delete_val(vertices, INVALID_VERTEX_ID);
}

void default_vertex_queue::init(const vertex_id_t buf[], size_t size, bool sorted)
{
	pthread_spin_lock(&lock);
	vertex_buf.clear();
	vpart_ps.clear();
	active_bitmap->clear();

	// The unpartitioned vertices
	std::vector<vertex_id_t> vertices;
	vertices.insert(vertices.end(), buf, buf + size);
	if (!sorted)
		std::sort(vertices.begin(), vertices.end());
	split_vertices(*index, part_id, vertices, vpart_ps);

	// The buffer contains the vertex Ids and we only store the location of
	// vertices in the local partition.
	vertex_buf.resize(vertices.size());
	index->get_vertices(vertices.data(), vertices.size(),
			compute_vertex_pointer::conv(vertex_buf.data()));

	buf_fetch_idx = scan_pointer(vertex_buf.size(), true);
	num_active = vertex_buf.size() + vpart_ps.size() * graph_conf.get_num_vparts();

	this->bitmap_fetch_idx = scan_pointer(0, true);
	pthread_spin_unlock(&lock);
	curr_vpart = 0;
}

void default_vertex_queue::init(worker_thread &t)
{
	pthread_spin_lock(&lock);
	vertex_buf.clear();
	vpart_ps.clear();
	assert(active_bitmap->get_num_set_bits() == 0);
	// This process only happens in a single thread, so we can swap
	// the two bitmap safely.
	std::unique_ptr<bitmap> tmp;
	tmp = std::move(active_bitmap);
	active_bitmap = std::move(t.next_activated_vertices);
	t.next_activated_vertices = std::move(tmp);
	num_active = active_bitmap->get_num_set_bits();

	// Get the vertically partitioned vertices that are activated.
	std::vector<vpart_vertex_pointer> vpart_ps_tmp(
			index->get_num_vpart_vertices(part_id));
	if (!vpart_ps_tmp.empty()) {
		index->get_vpart_vertex_pointers(part_id, vpart_ps_tmp.data(),
				vpart_ps_tmp.size());
		BOOST_FOREACH(vpart_vertex_pointer p, vpart_ps_tmp) {
			int part_id;
			off_t off;
			graph.get_partitioner()->map2loc(p.get_vertex_id(), part_id, off);
			if (active_bitmap->get(off)) {
				vpart_ps.push_back(p);
				active_bitmap->reset(off);
				num_active--;
			}
		}
		num_active += vpart_ps.size() * graph_conf.get_num_vparts();
	}

	bool forward = true;
	if (graph_conf.get_elevator_enabled())
		forward = graph.get_curr_level() % 2;
	bitmap_fetch_idx = scan_pointer(active_bitmap->get_num_longs(), forward);
	buf_fetch_idx = scan_pointer(0, true);
	pthread_spin_unlock(&lock);
	curr_vpart = 0;
}

void default_vertex_queue::fetch_from_map()
{
	assert(buf_fetch_idx.get_num_remaining() == 0);
	vertex_buf.clear();
	std::vector<vertex_id_t> ids;
	while (ids.size() < VERTEX_BUF_SIZE
			&& bitmap_fetch_idx.get_num_remaining() > 0) {
		size_t curr_loc = bitmap_fetch_idx.get_curr_loc();
		size_t new_loc = bitmap_fetch_idx.move(VERTEX_BUF_SIZE / NUM_BITS_LONG);
		// bitmap_fetch_idx points to the locations of longs.
		active_bitmap->get_reset_set_bits(min(curr_loc, new_loc) * NUM_BITS_LONG,
				max(curr_loc, new_loc) * NUM_BITS_LONG, ids);
	}
	std::vector<local_vid_t> local_ids(ids.size());
	for (size_t i = 0; i < ids.size(); i++)
		local_ids[i] = local_vid_t(ids[i]);
	vertex_buf.resize(local_ids.size());
	index->get_vertices(part_id, local_ids.data(), local_ids.size(),
			compute_vertex_pointer::conv(vertex_buf.data()));

	bool forward = true;
	if (graph_conf.get_elevator_enabled())
		forward = graph.get_curr_level() % 2;
	buf_fetch_idx = scan_pointer(vertex_buf.size(), forward);
}

void default_vertex_queue::fetch_vparts()
{
	if (curr_vpart >= graph_conf.get_num_vparts())
		return;

	assert(buf_fetch_idx.get_num_remaining() == 0);
	vertex_buf.clear();
	vertex_buf.resize(vpart_ps.size());
	index->get_vpart_vertices(part_id, curr_vpart, vpart_ps.data(),
			vpart_ps.size(), vertex_buf.data());
	curr_vpart++;

	// TODO Right now let's just scan the vertices in one direction.
	buf_fetch_idx = scan_pointer(vertex_buf.size(), true);
}

int default_vertex_queue::fetch(compute_vertex_pointer vertices[], int num)
{
	int num_fetched = 0;
	pthread_spin_lock(&lock);
	if (buf_fetch_idx.get_num_remaining() > 0) {
		int num_to_fetch = min(num, buf_fetch_idx.get_num_remaining());
		size_t curr_loc = buf_fetch_idx.get_curr_loc();
		size_t new_loc = buf_fetch_idx.move(num_to_fetch);
		memcpy(vertices, vertex_buf.data() + min(curr_loc, new_loc),
				num_to_fetch * sizeof(vertices[0]));
		num_fetched += num_to_fetch;
	}
	// We have fetched all we need.
	assert(num == num_fetched
			// the vertex buffer is empty.
			|| buf_fetch_idx.get_num_remaining() == 0);
	// If the vertex buffer is empty, let's get some from the bitmap.
	if (buf_fetch_idx.get_num_remaining() == 0) {
		// We start with unpartitioned vertices first and then
		// vertically partitioned vertices.
		fetch_from_map();
		if (vertex_buf.empty() && !vpart_ps.empty())
			fetch_vparts();
	}
	// If we still need some vertices.
	if (buf_fetch_idx.get_num_remaining() > 0 && num_fetched < num) {
		int fetch_again = min(num - num_fetched, buf_fetch_idx.get_num_remaining());
		size_t curr_loc = buf_fetch_idx.get_curr_loc();
		size_t new_loc = buf_fetch_idx.move(fetch_again);
		memcpy(vertices + num_fetched,
				vertex_buf.data() + min(curr_loc, new_loc),
				fetch_again * sizeof(vertices[0]));
		num_fetched += fetch_again;
	}
	num_active -= num_fetched;
	pthread_spin_unlock(&lock);
	return num_fetched;
}

void customized_vertex_queue::get_compute_vertex_pointers(
		const std::vector<vertex_id_t> &vertices,
		std::vector<vpart_vertex_pointer> &vpart_ps)
{
	sorted_vertices.resize(vertices.size()
			+ vpart_ps.size() * graph_conf.get_num_vparts());
	// Get unpartitioned vertices.
	index.get_vertices(vertices.data(), vertices.size(),
			compute_vertex_pointer::conv(sorted_vertices.data()));
	// Get vertically partitioned vertices.
	for (int i = 0; i < graph_conf.get_num_vparts(); i++) {
		off_t start = vertices.size() + i * vpart_ps.size();
		off_t end = start + vpart_ps.size();
		assert((size_t) end <= sorted_vertices.size());
		index.get_vpart_vertices(part_id, i, vpart_ps.data(), vpart_ps.size(),
				sorted_vertices.data() + start);
	}
}

void customized_vertex_queue::init(const vertex_id_t buf[],
		size_t size, bool sorted)
{
	pthread_spin_lock(&lock);
	sorted_vertices.clear();

	// The unpartitioned vertices
	std::vector<vertex_id_t> vertices;
	// The vertically partitioned vertices.
	std::vector<vpart_vertex_pointer> vpart_ps;
	vertices.insert(vertices.end(), buf, buf + size);
	if (!sorted)
		std::sort(vertices.begin(), vertices.end());
	split_vertices(index, part_id, vertices, vpart_ps);
	get_compute_vertex_pointers(vertices, vpart_ps);

	scheduler->schedule(sorted_vertices);
	bool forward = true;
	if (graph_conf.get_elevator_enabled())
		forward = graph.get_curr_level() % 2;
	this->fetch_idx = scan_pointer(sorted_vertices.size(), forward);
	pthread_spin_unlock(&lock);
}

void customized_vertex_queue::init(worker_thread &t)
{
	pthread_spin_lock(&lock);
	sorted_vertices.clear();
	std::vector<vertex_id_t> local_ids;
	t.next_activated_vertices->get_reset_set_bits(local_ids);

	// the bitmap only contains the locations of vertices in the bitmap.
	// We have to translate them back to vertex ids.
	std::vector<vertex_id_t> vertices(local_ids.size());
	for (size_t i = 0; i < local_ids.size(); i++) {
		vertex_id_t id;
		graph.get_partitioner()->loc2map(part_id, local_ids[i], id);
		vertices[i] = id;
	}
	std::vector<vertex_id_t>().swap(local_ids);
	std::vector<vpart_vertex_pointer> vpart_ps;
	split_vertices(index, part_id, vertices, vpart_ps);
	get_compute_vertex_pointers(vertices, vpart_ps);

	scheduler->schedule(sorted_vertices);
	bool forward = true;
	if (graph_conf.get_elevator_enabled())
		forward = graph.get_curr_level() % 2;
	fetch_idx = scan_pointer(sorted_vertices.size(), forward);
	pthread_spin_unlock(&lock);
}

worker_thread::worker_thread(graph_engine *graph,
		file_io_factory::shared_ptr graph_factory,
		file_io_factory::shared_ptr index_factory,
		vertex_program::ptr prog, vertex_program::ptr vpart_prog,
		int node_id, int worker_id, int num_threads,
		vertex_scheduler::ptr scheduler): thread("worker_thread", node_id)
{
	this->scheduler = scheduler;
	curr_compute = NULL;
	this->vprogram = std::move(prog);
	vprogram->init(graph, this);
	this->vpart_vprogram = std::move(vpart_prog);
	vpart_vprogram->init(graph, this);
	start_all = false;
	this->worker_id = worker_id;
	this->graph = graph;
	this->io = NULL;
	this->graph_factory = graph_factory;
	this->index_factory = index_factory;
	// We increase the allocator by 1M each time.
	// It shouldn't need to allocate much memory.
	msg_alloc = std::shared_ptr<slab_allocator>(new slab_allocator("graph-message-allocator",
			GRAPH_MSG_BUF_SIZE, 1024 * 1024, INT_MAX, get_node_id(),
			false /* init */, false /* pinned */, 5 /* local_buf_size*/));
	msg_processor = std::unique_ptr<message_processor>(new message_processor(
				*graph, *this, msg_alloc));
	balancer = std::unique_ptr<load_balancer>(new load_balancer(*graph, *this));
	switch(graph->get_graph_header().get_graph_type()) {
		case graph_type::DIRECTED:
			alloc = new vertex_compute_allocator<directed_vertex_compute>(graph, this);
			part_alloc = new vertex_compute_allocator<part_directed_vertex_compute>(
					graph, this);
			break;
		case graph_type::UNDIRECTED:
			alloc = new vertex_compute_allocator<vertex_compute>(graph, this);
			part_alloc = NULL;
			break;
#if 0
		case graph_type::TS_DIRECTED:
			alloc = new vertex_compute_allocator<ts_vertex_compute>(graph, this);
			part_alloc = new vertex_compute_allocator<part_ts_vertex_compute>(
					graph, this);
			break;
#endif
		default:
			assert(0);

	}
}

worker_thread::~worker_thread()
{
	delete alloc;
	delete part_alloc;
}

void worker_thread::init()
{
	// We should create these objects in the context of the worker thread,
	// so we can allocate memory for the objects on the same node as
	// the worker thread.
	next_activated_vertices = std::unique_ptr<bitmap>(
			new bitmap(graph->get_partitioner()->get_part_size(worker_id,
					graph->get_num_vertices()), get_node_id()));
	notify_vertices = std::unique_ptr<bitmap>(
			new bitmap(graph->get_partitioner()->get_part_size(worker_id,
					graph->get_num_vertices()), get_node_id()));
	if (scheduler)
		curr_activated_vertices = std::unique_ptr<active_vertex_queue>(
				new customized_vertex_queue(*graph, scheduler, worker_id));
	else
		curr_activated_vertices = std::unique_ptr<active_vertex_queue>(
				new default_vertex_queue(*graph, worker_id, get_node_id()));

	io = graph_factory->create_io(this);
	if (graph_conf.use_in_mem_index())
		index_reader = simple_index_reader::create(
				graph->get_in_mem_index(),
				graph->get_graph_header().get_graph_type() == graph_type::DIRECTED,
				this);
	else
		index_reader = simple_index_reader::create(
				index_factory->create_io(this),
				graph->get_graph_header().get_graph_type() == graph_type::DIRECTED,
				this);

	if (!started_vertices.empty()) {
		assert(curr_activated_vertices->is_empty());
		curr_activated_vertices->init(started_vertices, false);
		if (vinitializer) {
			BOOST_FOREACH(vertex_id_t id, started_vertices) {
				compute_vertex &v = graph->get_vertex(id);
				vinitializer->init(v);
			}
		}
		// Free the space used by the vector.
		started_vertices = std::vector<vertex_id_t>();
	}
	if (filter) {
		std::vector<vertex_id_t> local_ids;
		graph->get_partitioner()->get_all_vertices_in_part(worker_id,
				graph->get_num_vertices(), local_ids);

		std::vector<vertex_id_t> kept_ids;
		BOOST_FOREACH(vertex_id_t id, local_ids) {
			compute_vertex &v = graph->get_vertex(id);
			if (filter && filter->keep(v))
				kept_ids.push_back(id);
		}
		assert(curr_activated_vertices->is_empty());
		curr_activated_vertices->init(kept_ids, false);
		printf("worker %d has %ld vertices and activates %ld of them\n",
				worker_id, local_ids.size(), kept_ids.size());
	}
	// If a user wants to start all vertices.
	else if (start_all) {
		next_activated_vertices->set_all();
		assert(curr_activated_vertices->is_empty());
		curr_activated_vertices->init(*this);
		assert(next_activated_vertices->get_num_set_bits() == 0);
		if (vinitializer) {
			std::vector<vertex_id_t> local_ids;
			graph->get_partitioner()->get_all_vertices_in_part(worker_id,
					graph->get_num_vertices(), local_ids);
			BOOST_FOREACH(vertex_id_t id, local_ids) {
				compute_vertex &v = graph->get_vertex(id);
				vinitializer->init(v);
			}
		}
	}

	bool ret = graph->progress_first_level();
	assert(!ret);
}

void worker_thread::init_messaging(const std::vector<worker_thread *> &threads)
{
	vprogram->init_messaging(threads, msg_alloc);
	vpart_vprogram->init_messaging(threads, msg_alloc);
}

/**
 * This is to process the activated vertices in the current iteration.
 */
int worker_thread::process_activated_vertices(int max)
{
	if (max <= 0)
		return 0;

	process_vertex_buf.resize(max);
	int num = curr_activated_vertices->fetch(process_vertex_buf.data(), max);
	if (num == 0) {
		assert(curr_activated_vertices->is_empty());
		num = balancer->steal_activated_vertices(process_vertex_buf.data(),
				max);
	}
	if (num > 0) {
		num_activated_vertices_in_level.inc(num);
		graph->process_vertices(num);
	}

	for (int i = 0; i < num; i++) {
		compute_vertex_pointer info = process_vertex_buf[i];
		curr_vertex = info;
		// We execute the pre-run to determine if the vertex has completed
		// in the current iteration.
		vertex_program &curr_vprog = get_vertex_program(info.is_part());
		curr_compute = NULL;
		curr_vprog.run(*info);
		// If the user code doesn't generate a vertex_compute, we are done
		// with the vertex in this iteration.
		if (curr_compute == NULL)
			complete_vertex(info);
		curr_vertex = compute_vertex_pointer();
	}
	return num;
}

size_t worker_thread::enter_next_level()
{
	// We have to make sure all messages sent by other threads are processed.
	msg_processor->process_msgs();

	// If vertices have request the notification of the end of an iteration,
	// this is the place to notify them.
	if (notify_vertices->get_num_set_bits() > 0) {
		std::vector<vertex_id_t> vertex_buf;
		const size_t stride = 1024 * 64;
		for (size_t i = 0; i < notify_vertices->get_num_bits(); i += stride) {
			vertex_buf.clear();
			notify_vertices->get_reset_set_bits(i,
					min(i + stride, notify_vertices->get_num_bits()), vertex_buf);
			BOOST_FOREACH(vertex_id_t id, vertex_buf) {
				local_vid_t local_id(id);
				compute_vertex &v = graph->get_vertex(worker_id, local_id);
				vprogram->notify_iteration_end(v);
			}
		}
	}

	curr_activated_vertices->init(*this);
	assert(next_activated_vertices->get_num_set_bits() == 0);
	balancer->reset();
	msg_processor->reset();
	return curr_activated_vertices->get_num_vertices();
}

/**
 * This method is the main function of the graph engine.
 */
void worker_thread::run()
{
	while (true) {
		int num_visited = 0;
		int num;
		do {
			balancer->process_completed_stolen_vertices();
			num = process_activated_vertices(
					graph->get_max_processing_vertices()
					- max(get_num_vertices_processing(),
						io->num_pending_ios()));
			num_visited += num;
			msg_processor->process_msgs();
			if (index_reader->get_num_pending_tasks()
					< (size_t) graph->get_max_processing_vertices())
				index_reader->wait4complete(0);
			else
				index_reader->wait4complete(1);
			io->access(adj_reqs.data(), adj_reqs.size());
			adj_reqs.clear();
			if (io->num_pending_ios() == 0 && index_reader->get_num_pending_tasks() > 0)
				index_reader->wait4complete(1);
			io->wait4complete(min(io->num_pending_ios() / 10, 2));
			// If there are vertices being processed, we need to call
			// wait4complete to complete processing them.
		} while (get_num_vertices_processing() > 0
				// We still have vertices remaining for processing
				|| !curr_activated_vertices->is_empty()
				// Even if we have processed all activated vertices belonging
				// to this thread, we still need to process vertices from
				// other threads in order to balance the load.
				|| graph->get_num_remaining_vertices() > 0);
		assert(index_reader->get_num_pending_tasks() == 0);
		assert(io->num_pending_ios() == 0);
		assert(active_computes.size() == 0);
		assert(curr_activated_vertices->is_empty());
//		printf("worker %d visited %d vertices\n", worker_id, num_visited);
		assert(num_visited == num_activated_vertices_in_level.get());
		if (num_visited != num_completed_vertices_in_level.get()) {
			printf("worker %d: visits %d vertices and completes %ld\n",
					worker_id, num_visited, num_completed_vertices_in_level.get());
		}
		assert(num_visited == num_completed_vertices_in_level.get());

		// Now we have finished this level, we can progress to the next level.
		num_activated_vertices_in_level = atomic_number<long>(0);
		num_completed_vertices_in_level = atomic_number<long>(0);

		vprogram->flush_msgs();
		vpart_vprogram->flush_msgs();
		// We have to make sure all stolen vertices are returned to their owner
		// threads.
		balancer->process_completed_stolen_vertices();
		balancer->reset();
		bool completed = graph->progress_next_level();
//		printf("thread %d finish in a level, completed? %d\n", get_id(), completed);
		if (completed)
			break;
	}
	stop();
}

int worker_thread::steal_activated_vertices(compute_vertex_pointer vertices[], int num)
{
	// This method is called in the context of other worker threads,
	// curr_activated_vertices may not have been initialized. If so,
	// skip it.
	if (curr_activated_vertices == NULL)
		return 0;
	// We want to steal as much as possible, but we don't want
	// to overloaded by the stolen vertices.
	size_t num_steal = std::max(1UL,
			curr_activated_vertices->get_num_vertices() / graph->get_num_threads());
	num = curr_activated_vertices->fetch(vertices,
			std::min((size_t) num, num_steal));
	if (num > 0)
		// If the thread steals vertices from another thread successfully,
		// it needs to notify the thread of the stolen vertices.
		msg_processor->steal_vertices(vertices, num);
	return num;
}

void worker_thread::return_vertices(vertex_id_t ids[], int num)
{
	msg_processor->return_vertices(ids, num);
}

void worker_thread::complete_vertex(const compute_vertex_pointer v)
{
	std::unordered_map<compute_vertex *, vertex_compute *>::iterator it
		= active_computes.find(v.get());
	// It's possible that a vertex_compute isn't created for the active
	// compute_vertex.
	if (it != active_computes.end()) {
		vertex_compute *compute = it->second;
		// Since we have finished the computation on the vertex, we can
		// delete the vertex_compute now.
		active_computes.erase(it);
		compute->dec_ref();
		// It's possible that the vertex_compute is issued to SAFS.
		// In this case, SAFS will delete it.
		if (compute->get_ref() == 0) {
			assert(compute->get_num_pending() == 0);
			compute_allocator *alloc = compute->get_allocator();
			alloc->free(compute);
		}
	}

	num_completed_vertices_in_level.inc(1);
	// The vertex might be stolen from another thread. Now we have
	// finished processing it, we should return it to its owner thread.
	int part_id = graph->get_partitioner()->map(v->get_id());
	// We don't need to return verticalled partitioned vertices to their
	// owner because messages are processed in the main vertices and the
	// main vertices cannot be stolen by other threads.
	if (part_id != worker_id && !v.is_part()) {
		vertex_id_t id = v->get_id();
		balancer->return_vertices(&id, 1);
	}
}

void worker_thread::activate_vertex(local_vid_t id)
{
	next_activated_vertices->set(id.id);
}

void worker_thread::activate_vertex(vertex_id_t id)
{
	int part_id;
	off_t off;
	graph->get_partitioner()->map2loc(id, part_id, off);
	assert(part_id == worker_id);
	next_activated_vertices->set(off);
}

vertex_compute *worker_thread::get_vertex_compute(compute_vertex_pointer v)
{
	std::unordered_map<compute_vertex *, vertex_compute *>::const_iterator it
		= active_computes.find(v.get());
	if (it == active_computes.end()) {
		vertex_compute *compute = (vertex_compute *) alloc->alloc();
		compute->init(v);
		active_computes.insert(std::pair<compute_vertex *, vertex_compute *>(
					v.get(), compute));
		compute->inc_ref();
		curr_compute = compute;
	}
	else
		curr_compute = it->second;
	return curr_compute;
}
