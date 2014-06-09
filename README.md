FlashGraph
===========

FlashGraph is an [SSD](http://en.wikipedia.org/wiki/Solid-state_drive)-based
semi-external memory graph processing engine. FlashGraph stores vertex state
in memory and the edge lists of vertices on SSDs. It is optimized for a high-speed
SSD array or other high-speed flash storage. Its current implementation is tightly
integrated with `SAFS` (Set Associative File System) to take advantage of
the high I/O throughput of an SSD array. It has very short load time and
small memory consumption when processing very large graphs.
Thanks to the high-speed SSDs, it also has performance comparable to
in-memory graph engines.

It exposes a vertex-centric programming interface for users to express
graph algorithms and `FlashGraph` executes the user graph applications
in parallel and fetches data from SSDs. 

SAFSlib
========

`SAFSlib` is an open-source library that implements the design described
in the paper "Toward Millions of File System IOPS on Low-Cost, Commodity Hardware".
The goal of the library is to provide a filesystem-like interface
in the userspace for accessing SSD arrays. It is designed to eliminate overhead
in the block stack, without modifying the kernel. It can achieve the optimal performance
of a large SSD array in a [NUMA](http://en.wikipedia.org/wiki/Non-uniform_memory_access) machine.

The SAFSlib user manual can be found [here](https://docs.google.com/document/d/1OpsuLZw60MGCZAg4xO-j-1_AEWm3Yc2nqKKu8kXotkA/edit?usp=sharing).

The FlashGraph user manual can be found [here](https://github.com/icoming/FlashGraph/wiki/User-manual-of-FlashGraph).

