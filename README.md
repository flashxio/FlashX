FlashGraph is a SSD-based semi-external memory graph processing engine.
It stores vertex state in memory and the edge lists of vertices on SSDs.
It exposes a vertex-centric programming interface for users to express
graph algorithms and FlashGraph executes the user graph applications
in parallel and fetch data from SSDs. Its goal is to perform graph algorithms
on SSDs without much performance loss compared with in-memory graph engines.
The current implementation is tightly integrated with SAFS to take advantage
of the high I/O throughput of an SSD array.

SAFSlib is an open-source library that implements the design described
in the paper "Toward Millions of File System IOPS on Low-Cost, Commodity Hardware".
The goal of the library is to provide a filesystem-like interface
in the userspace for accessing SSD arrays. It is designed to eliminate overhead
in the block stack, without modifying the kernel. It can achieve optimal performance
of a large SSD array in a NUMA machine.

The user manual of SAFSlib can be found [here](https://docs.google.com/document/d/1OpsuLZw60MGCZAg4xO-j-1_AEWm3Yc2nqKKu8kXotkA/edit?usp=sharing).

The user manual of FlashGraph can be found [here](https://github.com/icoming/FlashGraph/wiki/User-manual-of-FlashGraph).

