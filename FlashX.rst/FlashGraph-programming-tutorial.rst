FlashGraph provides a flexible vertex-centric programming interface. In
this programming model, each vertex performs user-defined tasks
independently and interacts with other vertices as defined by program
logic. A vertex affects the state of others by sending messages to them
as well as activating them. Notably, FlashGraph allows a vertex to send
messages to any vertex in the graph. A vertex can also read the vertex
information of any vertex from SSDs as well as the state of any vertex
in memory.

A graph algorithm usually progresses in iterations. In each iteration,
the graph engine executes a user-defined task on each activated vertex.
An iteration ends when there are no more active vertices in the
iteration and no vertices have pending requests in the graph engine. An
algorithm ends when there aren’t active vertices in the next iteration.

Vertex program
==============

The most commonly way of implementing a graph algorithm in FlashGraph is
to define computation vertices by inheriting the ``compute_vertex``
class . Users define vertex state and implement three ``run`` methods in
the computation vertices, as shown below. FlashGraph executes the
``run`` method exactly once for each active vertex in an iteration; the
order of execution of this method on vertices is subject to scheduling
by FlashGraph. The execution of the ``run_on_vertex`` and
``run_on_message`` methods is event-driven. FlashGraph executes
``run_on_vertex`` when the edge list of a vertex requested by the
current vertex is ready in the page cache. FlashGraph executes
``run_on_message`` if the vertex receives messages from other vertices.
The ``run_on_message`` method may be executed even if a vertex is
inactive in an iteration. All examples assume ``using namespace fg;`` is
declared.

.. code:: {.cpp}

    class compute_vertex
    {
      // run only on the vertex state.
      void run(vertex_program &prog);

      // run on the edge list of a vertex
      void run_on_vertex(vertex_program &prog, page_vertex &vertex);

      // process a message.
      void run_on_message(vertex_program &prog, vertex_message &msg);
    };

Given the programming interface, breadth-first search can be simply
expressed as the code below. If a vertex has not been visited, it issues
a request to read its neighbor list in ``run`` and activates its
neighbors in ``run_on_vertex``. In this example, vertices do not need to
send messages to one another so we do not need to implement
``run_on_message``.

.. code:: {.cpp}

    class bfs_vertex: public vertex
    {
      bool has_visited;
      bfs_vertex() {
        has_visited = false;
      }

      void run(vertex_program &prog) {
        if (!has_visited) {
          vertex_id_t id = prog.get_vertex_id(*this);

          // Request vertex neighbor list from SAFS
          request_vertices(&id, 1);
          set_visited = true;
        }
      }

      void run_on_vertex(vertex_program &prog, page_vertex &vertex) {
        vertex_id_t dest_buf[];
        vertex.read_edges(dest_buf);
        prog.activate_vertices(dest_buf, num_dests);
      }

      void run_on_message(vertex_program &prog, vertex_message &msg) {
      }
    };

Initialize vertex state
=======================

There are two ways of initializing vertex state. Programmers can
initialize vertex state in the constructor of the user-defined
computation vertex. In the example of BFS, programmers only need to
initialize ``has_visited`` in the constructor of bfs\_vertex. For simple
graph algorithms, this is usually enough.

In a more complex case, a graph algorithm may require to execute the
vertex program multiple times or execute multiple vertex programs.
Therefore, it needs to set some vertices to a certain state or reset all
vertices. FlashGraph provides another mechanism to initialize vertex
state. Programmers need to implement the ``vertex_initiator`` interface,
shown as below. Users can pass a customized vertex initializer to the
graph engine by invoking its ``init_all_vertices()`` or its start
function. An example of using a customized vertex initializer can be
found in `single source shortest
path <https://github.com/icoming/FlashGraph/blob/graph-release/flash-graph/sssp/sssp.cpp>`__.

.. code:: {.cpp}

    class vertex_initiator
    {
    public:
        typedef std::shared_ptr<vertex_initiator> ptr;
        virtual void init(compute_vertex &) = 0;
    };

Interaction with other vertices
===============================

There are four ways for a vertex to interact with other vertices: a
vertex can send messages to other vertices; a vertex can read in-memory
vertex state of other vertices directly; a vertex can read the adjacency
list of other vertices from SSDs.

message passing
---------------

FlashGraph provides two methods for message passing:
``vertex_program::send_msg()`` and ``vertex_program::multicast_msg()``.
The former method is point-to-point communication between two vertices
and the second method allows a vertex to send a message to multiple
vertices. In most of the cases, multicast is used because multicast has
much smaller overhead and most graph algorithms require a vertex to send
the same message to all of its neighbors. A vertex gets notified of the
messages sent from other vertices through ``run_on_message()``.

All messages need to be inherited from the ``vertex_message`` class. Its
constructor takes two arguments: the size of the user-defined message
and the ``activate`` flag. When the ``activate`` flag is set, the
recipient vertices will be activated.

To reduce memory consumption, FlashGraph delivers messages to vertices
whenever it receives messages. Therefore, there is no guarantee of the
execution order of the three run methods. It is programmers'
responsibility of maintaining the correctness of vertex state. By
delivering messages to vertices immediately, we enable asynchronous
execution of graph algorithms. That is, an update to vertex state can be
immediately exposed to other vertices. It has advantage for some graph
algorithms because asynchronous execution can accelerate some graph
algorithms. This is very different from Pregel, which only delivers
messages to vertices at the end of an iteration.

Vertex activation
-----------------

A vertex can activate other vertices to run in the next iteration. There
are two ways of activating other vertices: with the dedicated methods
``vertex_program::activate_vertex`` and
``vertex_program::activate_vertices``; with the activate flag in
messages sent to other vertices.

Directed memory read
--------------------

We can get a reference to a vertex of a specified ID with
``graph_engine::get_vertex()``. This interface only works in a
shared-memory machine and may cause significant random memory access.
Therefore, this interface is not favored and should be used with
caution.

Access adjacency list from SSDs
-------------------------------

It takes two steps to read adjacency lists from SSDs: a vertex issues
read requests; the user-defined computation vertex gets notified through
its ``run_on_vertex()``. A vertex can read entire adjacency lists with
``compute_vertex::request_vertices()``. A directed vertex can read
partial adjacency lists with
``compute_directed_vertex::request_partial_vertices()``. In a partial
request, a directed vertex can request an in-edge list or an out-edge
list or both.

Data iterators
--------------

| FlashGraph defines very useful iterators for neighbor lists and edge
  attributes (for graphs that contain them).
| FlashGraph implements both sequential (Java-style) iterators and
  traditional STL-style iterators. Java-style iterators will improve
  performance in sequential access tasks and can be parameterized with a
  ``start`` and ``end`` positions for partial edge list numeration. For
  both examples assume the vertex has requested it's edge list in the
  ``run(vertex_program &prog)`` method.

Java-style iterators
~~~~~~~~~~~~~~~~~~~~

The code below shows how the Java-style iterators can be used to iterate
an edge list and access a data item in an attributed graph.

.. code:: {.cpp}

    typedef safs::page_byte_array::const_iterator<edge_data_type> data_iterator;
    typedef safs::page_byte_array::seq_const_iterator<edge_count> data_seq_iterator;

    void nmf_vertex::run(vertex_program &prog, const page_vertex &vertex) {
        // Iterator for neighbor IDs
        edge_seq_iterator neigh_it = vertex.get_neigh_seq_it(IN_EDGE);
        // Iterator for egde count (weight) attribute
        data_seq_iterator count_it =
            ((const page_directed_vertex&)vertex).get_data_seq_it<edge_count>(IN_EDGE);

        while (neigh_it.has_next()) {
            vertex_id_t nid = neigh_it.next();
            edge_count e = count_it.next();

            // Make use of `nid` and `e`
            std::cout << "Neighbor = " << nid << " Edge count = "
                          << e.get_count() << std::endl;
        }
    }

STL-style iterators
~~~~~~~~~~~~~~~~~~~

| The code below shows how the STL-style iterators can be used to
  iterate an edge list and access a data item in an attributed graph.
| ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
| void nmf\_vertex::run(vertex\_program &prog, const page\_vertex
  &vertex) {
| // Iterator for neighbor IDs
| edge\_iterator neigh\_it =
  vertex.get\_neigh\_begin(edge\_type::OUT\_EDGE);
| edge\_iterator neigh\_end =
  vertex.get\_neigh\_end(edge\_type::OUT\_EDGE);

::

    // Iterator for edge count (weight) attribute
    data_iterator count_it = 
          ((const page_directed_vertex&)vertex).get_data_begin(OUT_EDGE);
    data_iterator count_end = 
             ((const page_directed_vertex&)vertex).get_data_end(OUT_EDGE);

    for (; neigh_it != neigh_end; ++neigh_it) {
        vertex_id_t nid = *it;
        ++count_it;
        edge_count e = *count_it;

        // Make use of `nid` and `e`
        std::cout << "Neighbor = " << nid << " Edge count = "
                              << e.get_count() << std::endl;
    }

| }
| ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execute vertex program
======================

The code below executes the BFS program shown above. We create a
``graph_index`` object that contains the user-defined vertex state for
all vertices and create a ``graph_engine`` object that executes the user
code for the graph algorithm. In the case of BFS, the algorithm starts
on a single vertex. When a graph engine starts, the user code runs in
the worker threads inside the graph engine. We can invoke
``wait4complete`` to wait the graph algorithm to complete.

.. code:: {.cpp}

    graph_index::ptr index = NUMA_graph_index<bfs_vertex>::create(index_file);
    graph_engine::ptr graph = graph_engine::create(graph_file, index, configs);

    graph->start(&start_vertex, 1);
    graph->wait4complete();

Synchronous execution
=====================

By default, FlashGraph executes user-defined vertex computation
asynchronously. That is, the update to the vertex state is immediately
exposed to all other vertices in the same iteration. The asynchronous
execution can accelerate the convergence of many graph algorithms.
However, it is not deterministic and some graph algorithms need to be
executed synchronously.

FlashGraph also allows synchronous execution. FIXME
