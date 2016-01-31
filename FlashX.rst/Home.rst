FlashGraph is a semi-external memory graph processing engine, optimized
for a high-speed SSD array. It enables the processing of billion-node
graphs in a single machine and has performance comparable to or
exceeding in-memory graph engines such as
`PowerGraph <http://graphlab.org/>`__. Furthermore, it notably provides
significantly shorter loading time than PowerGraph due to the nature of
semi-external memory graph engines. FlashGraph is developed and tested
in Ubuntu 12.04 and Ubuntu 14.04.

FlashGraph targets programmers who need to implement parallel graph
algorithms to process large graphs. FlashGraph provides a flexible
vertex-centric programming interface so that users may express many
graph algorithms as well as custom algorithm-specific optimizations with
minimum effort. Users' graph algorithms are implemented in vertex
programs and FlashGraph executes the vertex programs in parallel and out
of core. FlashGraph hides the complexity of parallel and out-of-core
execution of graph algorithms so that users only need to implement
serial code and access data in memory. FlashGraph is implemented with
C++ and only exposes C++ programming interface right now. Therefore,
users need to have experience of C++ programming and Linux environment
programming.

We further implement graph algorithms on top of FlashGraph and pack them
as a graph library to benefit end users, much like
`iGraph <http://igraph.org/>`__ has done for serial algorithm
implementations. Currently, the graph algorithms distributed with
FlashGraph are listed
`here <https://github.com/icoming/FlashGraph/wiki/Graph-algorithms-in-FlashGraph>`__.
Many more graph algorithms that can scale to a graph with billions of
vertices will be implemented and released in future versions!

FlashGraph maintains user-defined vertex state in memory and accesses
the edge lists of vertices from SSDs through SAFS (Set Associative File
System), as shown in the architecture below. It is tightly integrated
with SAFS to take advantage of the high I/O throughput of an SSD array.
FlashGraph loads the edge lists of vertices to memory only when they are
required by the graph algorithms at runtime. Part of vertex programs,
provided by users , is executed in FlashGraph and part of it is executed
in SAFS to overlap I/O with computation and also minimize the overhead
of memory copy.

|Architecture|

.. |Architecture| image:: http://www.cs.jhu.edu/~zhengda/FlashGraph/arch.png
