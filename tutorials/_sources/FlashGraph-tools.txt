Convert graph format
====================

el2al takes a file or multiple files of edge lists in the text format
and converts them to the FlashGraph format, which is adjacency list in
the binary format. The edge list file(s) can be gzip'd. The gzip'd file
needs to have a filename extension ``.gz`` in order to be considered as
a gzip'd file. el2al generates two binary files: a graph file that
contains the adjacency lists of the graph and an index file that
contains the locations of vertices in the graph file. By default, el2al
keeps all intermediate data in memory.

el2al [options] adj\_list\_file index\_file edge\_list\_files (or
directories)

-  -u: undirected graph
-  -v: verify the created adjacency list
-  -t type: the type of edge data. Supported type: count, timestamp,
-  -m: merge multiple edge lists into a single graph.
-  -w: write the graph to a file
-  -T: the number of threads to process in parallel
-  -d: store intermediate data on disks

An example:

To convert a directed graph in the edge list format to the FlashGraph
format, we can run the following command:

``flash-graph/tools/el2al -w wiki-Vote.adj  wiki-Vote.index wiki-Vote.txt``

This converts an edge list in ``wiki-Vote.txt`` to the FlashGraph format
and generate two files: ``wiki-Vote.adj-v4`` and ``wiki-Vote.index-v4``.
``wiki-Vote.adj-v4`` contains the graph data and ``wiki-Vote.index-v4``
contains the index to the graph data.

To convert an undirected graph in the edge list format to the FlashGraph
format, we run:

``flash-graph/tools/el2al -w -u facebook.adj facebook.index facebook_combined.txt``

When converting a large graph, el2al uses a lot of memory if it keeps
all intermediate data in memory, and takes a long time to convert a
graph. el2al can keep intermediate data on disks and run in parallel.
When ``-d`` flag is used, el2al uses stxxl to store intermediate data on
disks. **Note: users may need to create .stxxl in the current directory
to configure the stxxl library where running el2al with -d flag.** The
`page <http://stxxl.sourceforge.net/tags/master/install_config.html>`__
describes how to configure stxxl. ``-T`` flag specifies the number of
threads used for graph format conversion.

For example, the following command converts a twitter graph with around
60 million vertices in the text edge list format to the adjacency list
format out of core and in parallel.

``el2al -w -T 32 -d twitter.adj twitter.index twitter_rv.net.gz``

FlashGraph configuration parameters
===================================

FIXME
