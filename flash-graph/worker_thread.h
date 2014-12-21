#ifndef __WORKER_THREAD_H__
#define __WORKER_THREAD_H__

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

#include <pthread.h>

#include <vector>
#include <unordered_map>

#include "graph_engine.h"
#include "bitmap.h"
#include "scan_pointer.h"

namespace safs
{
	class file_io_factory;
	class io_interface;
}

namespace fg
{

static const size_t MAX_ACTIVE_V = 1024;

class worker_thread;

/*
 * This data structure contains two data structures to represent active
 * vertices in an iteration. The bitmap is used when there are many active
 * vertices in an iteration; the vector is used when there are only a few
 * vertices in an iteration.
 */
class active_vertex_set
{
	// The fetech index in the active bitmap. It indicates the index of longs.
	bitmap active_map;
	scan_pointer bitmap_fetch_idx;

	std::vector<local_vid_t> active_v;

	struct local_vid_less {
		bool operator()(local_vid_t id1, local_vid_t id2) {
			return id1.id < id2.id;
		}
	};

	struct local_vid_eq {
		bool operator()(local_vid_t id1, local_vid_t id2) {
			return id1.id == id2.id;
		}
	};

	void set_bitmap(const local_vid_t ids[], int num) {
		for (int i = 0; i < num; i++)
			active_map.set(ids[i].id);
	}
public:
	active_vertex_set(size_t num_vertices, int node_id): active_map(
			num_vertices, node_id), bitmap_fetch_idx(0, true) {
	}

	void activate_all() {
		active_map.set_all();
	}

	void activate_vertex(local_vid_t id) {
		if (active_map.get_num_set_bits() > 0)
			active_map.set(id.id);
		else if (active_v.size() < MAX_ACTIVE_V)
			active_v.push_back(id);
		else {
			active_map.set(id.id);
			set_bitmap(active_v.data(), active_v.size());
			active_v.clear();
		}
	}

	void activate_vertices(const local_vid_t ids[], int num) {
		if (active_map.get_num_set_bits() > 0) {
			set_bitmap(ids, num);
		}
		else if (active_v.size() + num < MAX_ACTIVE_V)
			active_v.insert(active_v.end(), ids, ids + num);
		else {
			set_bitmap(ids, num);
			set_bitmap(active_v.data(), active_v.size());
			active_v.clear();
		}
	}

	vsize_t get_num_active_vertices() const {
		if (active_v.empty())
			return active_map.get_num_set_bits();
		else
			return active_v.size();
	}

	void finalize() {
		if (!active_v.empty()) {
			std::sort(active_v.begin(), active_v.end(), local_vid_less());
			std::vector<local_vid_t>::iterator new_end
				= std::unique(active_v.begin(), active_v.end(), local_vid_eq());
			size_t num_eles = new_end - active_v.begin();
			assert(num_eles <= active_v.size());
			active_v.resize(num_eles);
		}
	}

	void force_bitmap() {
		set_bitmap(active_v.data(), active_v.size());
		active_v.clear();
	}

	void reset_active_vertex(local_vid_t id) {
		assert(active_v.empty());
		active_map.reset(id.id);
	}

	bool is_active(local_vid_t id) const {
		assert(active_v.empty());
		return active_map.get(id.id);
	}

	void clear() {
		active_v.clear();
		active_map.clear();
		bitmap_fetch_idx = scan_pointer(0, true);
	}

	void set_dir(bool forward) {
		bitmap_fetch_idx = scan_pointer(active_map.get_num_longs(), forward);
	}

	void fetch_reset_active_vertices(size_t max_num,
			std::vector<local_vid_t> &local_ids);
	void fetch_reset_active_vertices(std::vector<local_vid_t> &local_ids);
};

/**
 * The queue for active vertices.
 */
class active_vertex_queue
{
public:
	virtual ~active_vertex_queue() {
	}

	virtual void init(const vertex_id_t buf[], size_t size, bool sorted) = 0;
	// This is the common case for iterations.
	virtual void init(worker_thread &) = 0;
	virtual int fetch(compute_vertex_pointer vertices[], int num) = 0;
	virtual bool is_empty() = 0;
	virtual size_t get_num_vertices() = 0;

	void init(const std::vector<vertex_id_t> &vec, bool sorted) {
		init(vec.data(), vec.size(), sorted);
	}
};

/**
 * This vertex queue is sorted based on the vertex ID.
 */
class default_vertex_queue: public active_vertex_queue
{
	static const size_t VERTEX_BUF_SIZE = 64 * 1024;
	pthread_spinlock_t lock;
	// It contains the offset of the vertex in the local partition
	// instead of the real vertex Ids.
	std::vector<compute_vertex_pointer> vertex_buf;
	// Pointers to the vertically partitioned vertices that are activated
	// in this iteration.
	std::vector<vpart_vertex_pointer> vpart_ps;
	int curr_vpart;
	std::unique_ptr<active_vertex_set> active_vertices;
	// The fetch index in the vertex buffer.
	scan_pointer buf_fetch_idx;
	graph_engine &graph;
	const graph_index &index;
	size_t num_active;
	int part_id;

	void fetch_from_map();
	void fetch_vparts();
public:
	default_vertex_queue(graph_engine &_graph, int part_id,
			int node_id): buf_fetch_idx(0, true), graph(_graph), index(
				_graph.get_graph_index()) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		num_active = 0;
		this->part_id = part_id;
		size_t num_local_vertices = _graph.get_partitioner()->get_part_size(
				part_id, _graph.get_num_vertices());
		this->active_vertices = std::unique_ptr<active_vertex_set>(
				new active_vertex_set(num_local_vertices, node_id));
		curr_vpart = 0;
	}

	virtual void init(const vertex_id_t buf[], size_t size, bool sorted);
	virtual void init(worker_thread &);
	virtual int fetch(compute_vertex_pointer vertices[], int num);

	virtual bool is_empty() {
		return num_active == 0;
	}

	virtual size_t get_num_vertices() {
		return num_active;
	}
};

class customized_vertex_queue: public active_vertex_queue
{
	pthread_spinlock_t lock;
	std::vector<compute_vertex_pointer> sorted_vertices;
	scan_pointer fetch_idx;
	vertex_scheduler::ptr scheduler;
	vertex_program::ptr vprog;
	graph_engine &graph;
	const graph_index &index;
	int part_id;

	void get_compute_vertex_pointers(const std::vector<vertex_id_t> &vertices,
		std::vector<vpart_vertex_pointer> &vpart_ps);
public:
	customized_vertex_queue(vertex_program::ptr vprog,
			vertex_scheduler::ptr scheduler, int part_id): fetch_idx(0,
				true), graph(vprog->get_graph()), index(
				graph.get_graph_index()) {
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		this->scheduler = scheduler;
		this->part_id = part_id;
		this->vprog = vprog;
	}

	void init(const vertex_id_t buf[], size_t size, bool sorted);
	void init(worker_thread &);

	int fetch(compute_vertex_pointer vertices[], int num) {
		pthread_spin_lock(&lock);
		int num_fetches = min(num, fetch_idx.get_num_remaining());
		if (num_fetches > 0) {
			size_t curr_loc = fetch_idx.get_curr_loc();
			size_t new_loc = fetch_idx.move(num_fetches);
			memcpy(vertices, sorted_vertices.data() + min(curr_loc, new_loc),
					num_fetches * sizeof(vertices[0]));
		}
		pthread_spin_unlock(&lock);
		return num_fetches;
	}

	bool is_empty() {
		pthread_spin_lock(&lock);
		bool ret = fetch_idx.get_num_remaining() == 0;
		pthread_spin_unlock(&lock);
		return ret;
	}

	size_t get_num_vertices() {
		pthread_spin_lock(&lock);
		size_t num = fetch_idx.get_num_remaining();
		pthread_spin_unlock(&lock);
		return num;
	}
};

class vertex_compute;
class steal_state_t;
class message_processor;
class load_balancer;
class simple_index_reader;

class worker_thread: public thread
{
	int worker_id;
	std::shared_ptr<safs::file_io_factory> graph_factory;
	std::shared_ptr<safs::file_io_factory> index_factory;
	std::shared_ptr<safs::io_interface> io;
	graph_engine *graph;
	const graph_index &index;

	std::unique_ptr<safs::compute_allocator> alloc;
	std::unique_ptr<safs::compute_allocator> merged_alloc;
	std::unique_ptr<safs::compute_allocator> sparse_alloc;
	vertex_program::ptr vprogram;
	// Vertex program on the vertically partitioned vertices.
	vertex_program::ptr vpart_vprogram;
	std::shared_ptr<simple_index_reader> index_reader;

	// This buffers the I/O requests for adjacency lists.
	std::vector<safs::io_request> adj_reqs;

	// When a thread process a vertex, the worker thread should keep
	// a vertex compute for the vertex. This is useful when a user-defined
	// compute vertex needs to reference its vertex compute.
	std::unordered_map<compute_vertex *, vertex_compute *> active_computes;
	// Determine whether the current vertex issues requests.
	bool req_on_vertex;
	// This points to the vertex that is currently being processed.
	compute_vertex_pointer curr_vertex;

	/**
	 * A vertex is allowed to send messages to other vertices.
	 * The message passing scheme between vertices are implemented as follows:
	 * all non-empty vertices (with edges) have a message holder;
	 * vertices are partitioned and each worker thread is responsible for
	 * a certin number of vertices;
	 * when a vertex issues messages, the thread that processes the vertex will
	 * redirect them to the right threads;
	 * a thread that receives messages process them and place them in
	 * the right message holder.
	 */

	std::unique_ptr<message_processor> msg_processor;
	std::unique_ptr<load_balancer> balancer;

	// This indicates the vertices that request the notification of the end
	// of an iteration.
	std::unique_ptr<bitmap> notify_vertices;
	// This is to collect vertices activated in the next level.
	std::unique_ptr<active_vertex_set> next_activated_vertices;
	// This contains the vertices activated in the current level.
	std::unique_ptr<active_vertex_queue> curr_activated_vertices;
	vertex_scheduler::ptr scheduler;

	// Indicate that we need to start all vertices.
	bool start_all;
	std::vector<vertex_id_t> started_vertices;
	std::shared_ptr<vertex_filter> filter;
	vertex_initializer::ptr vinitializer;

	// The buffer for processing activated vertex.
	embedded_array<compute_vertex_pointer> process_vertex_buf;

	// The number of activated vertices processed in the current level.
	atomic_number<long> num_activated_vertices_in_level;
	// The number of vertices completed in the current level.
	atomic_number<long> num_completed_vertices_in_level;

	/**
	 * Get the number of vertices being processed in the current level.
	 */
	int get_num_vertices_processing() const {
		return num_activated_vertices_in_level.get()
			- num_completed_vertices_in_level.get();
	}
	int process_activated_vertices(int max);
public:
	worker_thread(graph_engine *graph, std::shared_ptr<safs::file_io_factory> graph_factory,
			std::shared_ptr<safs::file_io_factory> index_factory, vertex_program::ptr prog,
			vertex_program::ptr part_prog, int node_id, int worker_id,
			int num_threads, vertex_scheduler::ptr scheduler,
			std::shared_ptr<slab_allocator> msg_alloc);

	~worker_thread();

	void init_messaging(const std::vector<worker_thread *> &threads,
			std::shared_ptr<slab_allocator> msg_alloc,
			std::shared_ptr<slab_allocator> flush_msg_alloc);

	void run();
	void init();

	/**
	 * When a vertex has been completed for the current iteration, this
	 * method is invoked to notify the worker thread, so that the worker
	 * thread can update its statistics on the number of completed vertices.
	 */
	void complete_vertex(const compute_vertex_pointer v);

	size_t enter_next_level();

	void start_vertices(const std::vector<vertex_id_t> &vertices,
			vertex_initializer::ptr initializer) {
		this->vinitializer = initializer;
		started_vertices = vertices;
	}

	void start_all_vertices(vertex_initializer::ptr init) {
		start_all = true;
		this->vinitializer = init;
	}

	void start_vertices(std::shared_ptr<vertex_filter> filter) {
		this->filter = filter;
	}

	compute_vertex_pointer get_curr_vertex() const {
		return curr_vertex;
	}
	void start_run_vertex(compute_vertex_pointer v) {
		assert(!curr_vertex.is_valid());
		curr_vertex = v;
		req_on_vertex = false;
	}
	bool finish_run_vertex(compute_vertex_pointer v) {
		assert(curr_vertex.is_valid());
		assert(curr_vertex.get() == v.get());
		curr_vertex = compute_vertex_pointer();
		return req_on_vertex;
	}

	void request_on_vertex(vertex_id_t id) {
		req_on_vertex = true;
	}
	vertex_compute *get_vertex_compute(compute_vertex_pointer v);

	/**
	 * Activate the vertex in its own partition for the next iteration.
	 */
	void activate_vertex(local_vid_t id) {
		next_activated_vertices->activate_vertex(id);
	}

	void activate_vertices(const local_vid_t ids[], int num) {
		next_activated_vertices->activate_vertices(ids, num);
	}

	void request_notify_iter_end(local_vid_t id) {
		notify_vertices->set(id.id);
	}

	int steal_activated_vertices(compute_vertex_pointer vertices[], int num);
	void return_vertices(vertex_id_t ids[], int num);

	size_t get_num_local_vertices() const {
		return graph->get_partitioner()->get_part_size(worker_id,
					graph->get_num_vertices());
	}

	int get_worker_id() const {
		return worker_id;
	}

	vertex_program &get_vertex_program(bool part) {
		return part ? *vpart_vprogram : *vprogram;
	}

	graph_engine &get_graph() {
		return *graph;
	}

	message_processor &get_msg_processor() {
		return *msg_processor;
	}

	simple_index_reader &get_index_reader() {
		return *index_reader;
	}

	void issue_io_request(safs::io_request &req) {
		adj_reqs.push_back(req);
	}

	size_t get_activates() const {
		return curr_activated_vertices->get_num_vertices();
	}

	safs::compute_allocator &get_merged_compute_allocator() {
		return *merged_alloc;
	}

	safs::compute_allocator &get_sparse_compute_allocator() {
		return *sparse_alloc;
	}

	int get_stolen_vertex_part(const compute_vertex &v) const;

	friend class load_balancer;
	friend class default_vertex_queue;
	friend class customized_vertex_queue;
};

}

#endif
