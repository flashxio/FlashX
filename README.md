

FlashGraph
===========

FlashGraph is an [SSD](http://en.wikipedia.org/wiki/Solid-state_drive)-based
semi-external memory graph processing engine, and is optimized for a high-speed
SSD array. It enables us to process a billion-node graph in a single machine
and has performance comparable to or exceed in-memory graph engines such as
[PowerGraph](http://graphlab.org/). It also has very short loading time.

As a semi-external memory graph engine, FlashGraph stores vertex state
in memory and the edge lists of vertices on SSDs. FlashGraph loads the edge
lists of vertices to memory whenever required by the graph algorithms.
Its current implementation is tightly
integrated with `SAFS` (Set Associative File System) to take advantage of
the high I/O throughput of an SSD array.

FlashGraph is deisnged to help users develop parallal and out-of-core
implementations of graph algorithms with minimum effort.
Flashgraph is implemented with C++ and provides a vertex-centric programming
interface for users to express graph algorithms. Users encapsulate their graph
algorithms in vertex programs and `FlashGraph` executes the vertex programs
in parallel and fetches data from SSDs for the vertex programs.

[FlashGraph Quickstart](https://github.com/icoming/FlashGraph/wiki/FlashGraph-Quick-Start-Guide)

[FlashGraph User manual](https://github.com/icoming/FlashGraph/wiki/User-manual-of-FlashGraph).

The detailed design of FlashGraph is documented in the paper
"FlashGraph: Processing Billion-Node Graphs on an Array of Commodity SSDs",
submitted to Supercomputing'14.

SAFSlib
========

`SAFSlib` is an open-source library that provides a filesystem-like interface
in the userspace to help users access a large SSD array in a [NUMA](http://en.wikipedia.org/wiki/Non-uniform_memory_access) machine.
It is designed to eliminate overhead in the block stack, without modifying the kernel,
and achieves the maximal performance of a large SSD array in a NUMA machine.
FlashGraph is an application to demonstrate the power of SAFS.

[SAFS user manual](https://docs.google.com/document/d/1OpsuLZw60MGCZAg4xO-j-1_AEWm3Yc2nqKKu8kXotkA/edit?usp=sharing).

The detailed design of SAFS is documented in the paper
"[Toward Millions of File System IOPS on Low-Cost, Commodity Hardware](http://dl.acm.org/citation.cfm?id=2503225&dl=ACM&coll=DL&CFID=350399128&CFTOKEN=49883861)".
