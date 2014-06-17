There are main components in the repository: FlashGraph and SAFS.

FlashGraph
===========

FlashGraph is a semi-external memory graph processing engine, optimized for a high-speed
SSD array. It helps us develop parallel graph algorithms with minimum efforts and
executes users' graph algorithms in parallel and out of core.
It enables us to process a billion-node graph in a single machine
and has performance comparable to or exceed in-memory graph engines such as
[PowerGraph](http://graphlab.org/).

[FlashGraph Quick start guide](https://github.com/icoming/FlashGraph/wiki/FlashGraph-Quick-Start-Guide)

[FlashGraph User manual](https://github.com/icoming/FlashGraph/wiki/FlashGraph-User-Manual).

SAFSlib
========

`SAFSlib` is an open-source library that provides a filesystem-like interface
in the userspace to help users access a large SSD array in a
[NUMA](http://en.wikipedia.org/wiki/Non-uniform_memory_access) machine.
It is designed to eliminate overhead in the block subsystem in Linux, without modifying the kernel,
and achieves the maximal performance of a large SSD array in a NUMA machine.
FlashGraph is an application to demonstrate the power of SAFS.

[SAFS user manual](https://docs.google.com/document/d/1OpsuLZw60MGCZAg4xO-j-1_AEWm3Yc2nqKKu8kXotkA/edit?usp=sharing).

The detailed design of SAFS is documented in the paper
"[Toward Millions of File System IOPS on Low-Cost, Commodity Hardware](http://dl.acm.org/citation.cfm?id=2503225&dl=ACM&coll=DL&CFID=350399128&CFTOKEN=49883861)".



Contact
========

Mailing list: flashgraph-dev@googlegroups.com
