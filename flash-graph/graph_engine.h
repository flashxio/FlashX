#ifndef __GRAPH_ENGINE_H__
#define __GRAPH_ENGINE_H__

/*
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

#include "thread.h"
#include "container.h"
#include "concurrency.h"
#include "slab_allocator.h"
#include "io_interface.h"

#include "vertex.h"
#include "vertex_index.h"
#include "trace_logger.h"
#include "messaging.h"
#include "vertex_interpreter.h"
#include "partitioner.h"
#include "graph_index.h"
#include "graph_config.h"
#include "vertex_request.h"
#include "vertex_program.h"

/**
 * The size of a message buffer used to pass vertex messages to other threads.
 */
const int GRAPH_MSG_BUF_SIZE = PAGE_SIZE * 4;

class graph_engine;
class vertex_request;

/**
  * \brief Class from which users' vertex-centric programs should inherit.
  *         Serial code written when implementing <I>run*</I>.
  *         methods here will be run in parallel within the graph engine.
*/
class compute_vertex
{
	vertex_id_t id;
public:
	/**
	 * \brief The constructor called by graph_engine to create vertex
	 * state.
	 *
	 * <I>Note that users never need to explictly call this ctor.</I>
	 *         
	 */
	compute_vertex(vertex_id_t id) {
		this->id = id;
	}

	/**
	 * \brief This allows a vertex to request the adjacency lists of vertices
	 *        in the graph.
     *
	 * \param ids The IDs of vertices.
     * \param num The number of vertex IDs you are requesting.
	 */
	void request_vertices(vertex_id_t ids[], size_t num);

	/**
	 * \brief This requests the vertex headers. It mainly contains the number
	 *        of edges (the number of in-edges and out-edges for directed
	 *        vertices).
	 * \param ids The IDs of vertices.
     * \param num The number of vertex IDs you are requesting.
	 */
	void request_vertex_headers(vertex_id_t ids[], size_t num);

    /**
     * \brief Get its own vertex ID.
     *
     * \return The current vertex's ID in the graph.
    */
	vertex_id_t get_id() const {
		return id;
	}
    
    /**
     * \brief Allows a vertex to perform a task at the end of every iteration. 
     * \param prog The vertex program associated with running the graph algorithm.
     */
	void notify_iteration_end(vertex_program & prog) {
	}

	void run_on_vertex_header(vertex_program &prog,
			const vertex_header &header) {
		assert(0);
	}
};

/**
 * \brief A directed version of the <c>compute_vertex</c> class that
 *         users inherit from when using the FlashGraph engine.
 *
 */
class compute_directed_vertex: public compute_vertex
{
public:
    /**
     * \brief The constructor callled by the grah engne.
     */
	compute_directed_vertex(vertex_id_t id): compute_vertex(id) {
	}

	/**
	 * \brief This allows a vertex to request partial vertices in the graph. <br>
     *  **Defn: ""partial vertices""** -- Part of a vertex's data. Specifically
     *          either a vertex's *in-edges* OR *out-edges*. This eliminates the overhead
     *          of bringing both *in* and *out* edges into the page cache when an algorithm
     *          only requires one of the two.
	 * \param reqs This is an array corresponding to the vertices you are requesting and defines
     *              which part of the vertex you want (e.g `IN_EDGE`).
     * \param num the number of elements in `reqs`.
	 */
	void request_partial_vertices(directed_vertex_request reqs[], size_t num);
};

/**
 *  \brief Time series compute vertex used for time series graph analytics.
 */
class compute_ts_vertex: public compute_vertex
{
public:
    /**
     * \brief The constructor. Nothing initialized here.
     */
	compute_ts_vertex(vertex_id_t id): compute_vertex(id) {
	}

#if 0
	/**
     * \brief This allows a vertex to request partial vertices in the graph. <br>
     *  **Defn: ""partial vertices""** -- Part of a vertex's data. Specifically
     *          either a vertex's edges in the specified time range. This eliminates the overhead
     *          of bringing both all edges into the page cache when an algorithm
     *          only requires edges in a certain time range.
	 * \param reqs This is an array corresponding to the vertices you are requesting and defines
     *              which time range you want.
     * \param num the number of elements in `reqs`.
	 * This allows a vertex to request partial vertices in the graph.
	 */
	void request_partial_vertices(ts_vertex_request reqs[], size_t num);
#endif
};

/**
 * \breif Order the position of vertex processing via this scheduler.
 */
class vertex_scheduler
{
public:
	typedef std::shared_ptr<vertex_scheduler> ptr; /** Smart pointer for object access.*/

    /**
     * \brief Implement this method in order to customize the vertex schedule.
     *
     *  \param vertices An `std::vector` of vertex IDs defining the order in which
     *          vertices in the graph should be processed.
     */
	virtual void schedule(std::vector<vertex_id_t> &vertices) = 0;
};

/**
 * \brief When the graph engine starts, a user can use this filter to decide
 * what vertices are activated for the first time.
 */
class vertex_filter
{
public:
    
    /**
     * \brief The method defines which vertices are active in the first iteration.
     *           If this method returns `true` the vertex will be activated, whereas
     *           `false` means the vertex will be inactive in the first iteration.
     * \param v A single `compute_vertex`. **Note:** This parameter is provided by the graph
     *          engine. A user need not ever provide it.
     */
	virtual bool keep(compute_vertex &v) = 0;
};

/**
 * \brief A user may be decide to initialize individual vertex state in a custom way not
 *          expressible via the vertex constructor. This provides that capability.
 */
class vertex_initializer
{
public:
	typedef std::shared_ptr<vertex_initializer> ptr; /** Type provides access to the object */
    
    /**
     * \brief Initialization method to initialize the given vertex.
     */
	virtual void init(compute_vertex &) = 0;
};


class graph_engine;

/**
 * \brief Parallized query of the vertex state of all vertices in the graph.
 *        Each worker thread gets an instance of the query and the per-thread
 *        query results will be merged in the end.
 *          Inherit from this class to run queries in parallel.
 */
class vertex_query
{
public:
	typedef std::shared_ptr<vertex_query> ptr; /** Type provides access to the object */
    
    /**
     * \brief This method is executed on vertices in parallel and contains any user defined
     *         code.
     */
	virtual void run(graph_engine &, compute_vertex &v) = 0;
    
    /**
     *  \brief All vertex results may be merged (not specially combined but any custom operation).
     *          <br>. This for instance can be used to aggregate (add, subtract, max etc.)
     *              a user defined data member for the class.
     * \param graph The graph engine you are concerned with.
     * \param q A pointer to the vertex query.
     */
	virtual void merge(graph_engine &graph, vertex_query::ptr q) = 0;
    
    /**
     * \brief Implements a copy constructor. The graph engine uses this method
	 *        to create an instance of this query for each thread.
     * Used internally by graph engine as a generic ctor in lieu of using a templated class.
     */
	virtual ptr clone() = 0;
};

class worker_thread;

/**
 * \brief This is the class that coordinates how & where algorithms are run.
 *          It can be seen as the central organ of FlashGraph.
*/
class graph_engine
{
	static std::atomic<size_t> init_count;

	int vertex_header_size;
	graph_header header;
	graph_index::ptr vertices;
	vertex_index::ptr vindex;
	std::unique_ptr<ext_mem_vertex_interpreter> interpreter;
	vertex_scheduler::ptr scheduler;

	// The number of activated vertices that haven't been processed
	// in the current level.
	atomic_number<size_t> num_remaining_vertices_in_level;
	atomic_integer level;
	volatile bool is_complete;

	// These are used for switching queues.
	pthread_mutex_t lock;
	pthread_barrier_t barrier1;
	pthread_barrier_t barrier2;

	int num_nodes;
	std::vector<worker_thread *> worker_threads;
	std::vector<vertex_program::ptr> vprograms;

	trace_logger::ptr logger;
	file_io_factory::shared_ptr graph_factory;
	file_io_factory::shared_ptr index_factory;
	int max_processing_vertices;

	// The time when the current iteration starts.
	struct timeval start_time;

	void init_threads(vertex_program_creater::ptr creater);
protected:
    /**
     * \brief Constructor usable by inheriting classes.
     * \param graph_file The path to the graph file on disk.
     * \param index The path to the graph index file on disk.
     * \param configs The path to the configuration file on disk.
     */
	graph_engine(const std::string &graph_file, graph_index::ptr index,
			const config_map &configs);
public:
	typedef std::shared_ptr<graph_engine> ptr; /** Smart pointer for object access.*/

	static void init_flash_graph(const config_map &configs);
	static void destroy_flash_graph();
    
    /**
     * \brief Use this method to in lieu of a constructor to create a graph object.
     * \param graph_file The path to the graph file on disk.
     * \param index The path to the graph index file on disk.
     * \param configs The path to the configuration file on
     */
	static graph_engine::ptr create(const std::string &graph_file,
			graph_index::ptr index, const config_map &configs) {
		return graph_engine::ptr(new graph_engine(graph_file, index, configs));
	}

    /**
     * \brief Class destructor
     */
	~graph_engine();

	/*
	 * The following four variants of get_vertex return the same compute_vertex
	 * but with slightly lower overhead.
	 * These can only be used in a shared machine.
	 */
    
    /**
	 * \brief Permits any vertex to pull any other vertex's state.<br>
	 * **NOTE:** *This can only be used in a shared machine.*
     * \param id the unique vertex ID.
     * \return The `compute_vertex` requested by id.
	 */
	compute_vertex &get_vertex(vertex_id_t id) {
		return vertices->get_vertex(id);
	}
    
    /** \internal */
	compute_vertex &get_vertex(int part_id, local_vid_t id) {
		return vertices->get_vertex(part_id, id);
	}

    /**
	 * \brief Permits any vertex to pull any set of the vertex state of
	 *        other vertices.<br>
	 * **NOTE:** *This can only be used in a shared machine.*
	 */
	size_t get_vertices(const vertex_id_t ids[], int num, compute_vertex *v_buf[]) {
		return vertices->get_vertices(ids, num, v_buf);
	}
    
    /** \internal */
	size_t get_vertices(int part_id, const local_vid_t ids[], int num,
			compute_vertex *v_buf[]) {
		return vertices->get_vertices(part_id, ids, num, v_buf);
	}
    
    /**
     * \brief Get the maximum vertex ID in the graph.
     * \return The maximum vertex ID in the graph.
     */
	vertex_id_t get_max_vertex_id() const {
		return vertices->get_max_vertex_id();
	}
    
    /**
     * \brief Get the minimum vertex ID in the graph.
     * \return The the minimum vertex ID in the graph.
     */
	vertex_id_t get_min_vertex_id() const {
		return vertices->get_min_vertex_id();
	}
    

    /**
     * \brief Get the number of vertices in the graph.
     * \return The number of vertices in the graph.
     */
	size_t get_num_vertices() const {
		return vertices->get_num_vertices();
	}
    
    /**
     * \brief Tell a user if the graph is directed or not.
     * \return true if the graph is directed, else false.
     */
	bool is_directed() const {
		return header.is_directed_graph();
	}
    
    /**
     * \brief Get the graph header info.
     * \return The graph header with all its associated metadata.
     */
	const graph_header &get_graph_header() const {
		return header;
	}
    
    /**
     * \brief Set the graph computation to use a custom vertex scheduler.
     * \param scheduler The user-defined vertex scheduler.
     */
	void set_vertex_scheduler(vertex_scheduler::ptr scheduler);
    
    /**
     * \brief Start the graph engine and begin computation on a subset of vertices.
     * \param filter A user defined `vertex_filter` which specifies which vertices
     *      are activated the next iteratoin.
	 * \param creater A creator that creates user-defined vertex program.
	 *                By default, a graph engine creates its own default
	 *                vertex program.
     */
	void start(std::shared_ptr<vertex_filter> filter,
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());
    
    /**
     * \brief Start the graph engine and begin computation on a subset of vertices.
     * \param ids The vertices that should be activated for the first iteration.
     * \param init An initializer used to alter the state of vertices activated
	 *             in the first iteration.
	 * \param creater A creator that creates user-defined vertex program.
	 *                By default, a graph engine creates its own default
	 *                vertex program.
     */
	void start(vertex_id_t ids[], int num,
			vertex_initializer::ptr init = vertex_initializer::ptr(),
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());
    
    /**
     * \brief Start the graph engine and begin computation on **all** vertices.
     * \param init An initializer used to alter the state of vertices activated
	 *             in the first iteration.
	 * \param creater A creator that creates user-defined vertex program.
	 *                By default, a graph engine creates its own default
	 *                vertex program.
     */
	void start_all(vertex_initializer::ptr init = vertex_initializer::ptr(),
			vertex_program_creater::ptr creater = vertex_program_creater::ptr());
    
    /**
     * \brief Synchronization barrier that waits for the graph algorithm to
	 *        complete.
     */
	void wait4complete();

	/**
	 * \brief This method preloads the entire graph to the page cache.
	 *        If the page cache is smaller than the graph, only the first part
	 *        of the graph image (the same size as the page cache) is preloaded
	 *        to the page cache.
	 */
	void preload_graph();

	/**
	 * \brief Allows users to initialize vertices to certain state.
     * \param ids The vertex ID for which you want initialize.
     * \param num The number of vertices you intend to initialize.
     * \param init An initializer used to alter the state of a vertex.
	 */
	void init_vertices(vertex_id_t ids[], int num, vertex_initializer::ptr init);
    
    /**
     * \brief Allows users to initialize **all** vertices to certain state.
     * \param init An initializer used to alter the state of a vertex.
     */
	void init_all_vertices(vertex_initializer::ptr init);

	/**
	 * \brief Allows users to query the information on the state of all vertices.
     * \param query The `vertex_query` you wish to apply to the graph.
	 */
	void query_on_all(vertex_query::ptr query);
    
    /**
     * \brief Return The per-thread vertex programs used by the graph engine.
     * \param programs An empty vector that will eventully contain all the vertex progrmas.
     */
	void get_vertex_programs(std::vector<vertex_program::ptr> &programs) {
		programs = vprograms;
	}

	/**
	 * \brief This returns the current iteration number in the graph engine.
     * \return The current iteration number.
	 */
	int get_curr_level() const {
		return level.get();
	}

	/**
	 * The methods below should be used internally.
	 */

	/**
     * \internal
	 * The algorithm progresses to the next level.
	 * It returns true if no more work can progress.
	 */
	bool progress_next_level();
    
    /** \internal*/
	trace_logger::ptr get_logger() const {
		return logger;
	}

	/**
     * \internal
	 * Get the file id where the graph data is stored.
	 */
	int get_file_id() const {
		return graph_factory->get_file_id();
	}
    
    /**\internal */
	ext_mem_vertex_interpreter &get_vertex_interpreter() const {
		return *interpreter;
	}
    
    /**\internal */
	const graph_partitioner *get_partitioner() const {
		return &vertices->get_partitioner();
	}
    
    /**\internal */
	int get_num_threads() const {
		return worker_threads.size();
	}
    
    /**\internal */
	worker_thread *get_thread(int idx) const {
		return worker_threads[idx];
	}

	/*
	 * The following two methods keep track of the number of active vertices
	 * globally in the current iteration.
	 */
    
    /**
    * \internal
	* We have processed the specified number of vertices.
    */
	void process_vertices(int num) {
		num_remaining_vertices_in_level.dec(num);
	}

	/**
     * \internal Get the number of activated vertices that still haven't been
     * processed in the current level.
     */
	size_t get_num_remaining_vertices() const {
		return num_remaining_vertices_in_level.get();
	}

	int get_vertex_header_size() const {
		return vertex_header_size;
	}

	/*
	 * \internal Get the in-memory vertex index.
	 */
	vertex_index::ptr get_in_mem_index() const {
		assert(vindex);
		return vindex;
	}

	void set_max_processing_vertices(int max) {
		max_processing_vertices = max;
	}

	int get_max_processing_vertices() const {
		return max_processing_vertices;
	}
};

#endif
