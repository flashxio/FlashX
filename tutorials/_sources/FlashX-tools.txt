FlashX tools
============

FlashX provides three tools:

-  ``utils/SAFS-util``: helps users to access files in SAFS.
-  ``matrix/utils/el2fg``: constructs a FlashGraph graph from the text
   edge list.
-  ``matrix/utils/fg2fm``: constructs a FlashMatrix sparse matrix from a
   FlashGraph graph.

Interact with SAFS
------------------

utils/SAFS-util is a tool provided by SAFS that helps users to interact
with SAFS. It provides a set of commands to access SAFS.

::

    SAFS-util conf_file command ...
    The supported commands are
        create file_name size: create a file with the specified size
        delete file_name: delete a file
        help: print the help info
        list: list existing files in SAFS
        load file_name [ext_file]: load data to the file
        load_part file_name ext_file part_id: load part of the file to SAFS
        verify file_name [ext_file]: verify data in the file
        export file_name ext_file: export an SAFS file to Linux filesystem
        info file_name: show the information of an SAFS file
        rename file_name new_name: rename an SAFS file

The most useful commands are:

-  create: create an SAFS file with the specified size
-  delete: delete an SAFS file
-  list: list all files in SAFS
-  load: load data to a SAFS file named ``file_name``. Users can load
   data from a Linux file to the SAFS file by providing a second
   argument. If a Linux file isn't provided, the SAFS file is filled
   with 0.
-  export: export an SAFS file to Linux filesystem
-  info file\_name: show the metadata information of an SAFS file
-  rename file\_name new\_name: rename an SAFS file

Construct FlashGraph graphs
---------------------------

| matrix/utils/el2fg takes a file of edge lists in the text format and
  converts them to the FlashGraph format, which is adjacency list in the
  binary format. The edge list file(s) can be gzip'd. The gzip'd file
  needs to have a filename extension ``.gz`` in order to be considered
  as a gzip'd file. el2fg takes a FlashX configuration file, a text edge
  list file and a graph name, and outputs two files: graph\_name.adj,
  which contains the adjacency lists of the graph, and
  graph\_name.index, which contains the locations of the adjacency lists
  in the graph file.
| By default, el2fg keeps all intermediate data in memory.

::

    matrix/utils/el2fg [options] conf_file edge_file graph_name
    -u: undirected graph
    -U: unqiue edges
    -e: use external memory
    -s size: sort buffer size
    -g size: groupby buffer size
    -t type: the edge attribute type

The most important flags are shown as follows:

-  -u: To construct an undirected graph, users have to explicitly enable
   this flag.
-  -U: enable this flag to remove redundant edges in a graph.
-  -e: enable this flag to use disks to construct a graph. It requires
   users to configure SAFS correctly.
-  -t type: users can specify the edge attribute type. When this option
   is specified, the input edge list has to have three columns and the
   third column provides the edge attribute. Currently, four attribute
   types are supported. "I": 32-bit integer attributes, "L": 64-bit
   integer attributes, "F": single-precision float-point, "D":
   double-precision float-point.

An example:

| To convert a directed graph
  `wiki-Vote <http://snap.stanford.edu/data/wiki-Vote.html>`__ in the
  edge list format to the FlashGraph format, we can run the following
  command:
| ``matrix/utils/el2fg matrix/conf/run_test.txt wiki-Vote.txt wiki-Vote``
| This converts an edge list in ``wiki-Vote.txt`` to the FlashGraph
  format and generate two files: ``wiki-Vote.adj`` and
  ``wiki-Vote.index``.

| To convert an undirected graph
  `facebook <http://snap.stanford.edu/data/egonets-Facebook.html>`__ in
  the edge list format to the FlashGraph format, we run:
| ``matrix/utils/el2fg matrix/conf/run_test.txt -u facebook_combined.txt facebook``

Construct FlashMatrix sparse matrices
-------------------------------------

matrix/utils/fg2fm converts a graph in the FlashGraph format to a sparse
matrix in the FlashMatrix format. The FlashMatrix is necessary to
support very efficient sparse matrix dense matrix multiplication (SpMM).
SpMM is much more efficient with the FlashMatrix format than the
FlashGraph format.

::

    fg2fm [options] conf_file graph_file index_file matrix_name
    -h height: the height of a 2D-partitioned matrix block
    -w width: the width of a 2D-partitioned matrix block
    -v: verify the generated matrix image
    -g size: the groupby buffer size on vector vectors
    -t type: the type of non-zero entries

Users don't need to provide the options to ``fg2fm`` most of time. The
only useful option for some users is ``-t``, which tells ``fg2fm`` the
type of the edge attribute in the FlashGraph format.

| ``fg2fm`` detects whether the provided graph is directed. For an
  undirected graph, it outputs two files: ``matrix_name.mat``, which
  contains the sparse matrix data, and ``matrix_name.mat_idx``, which
  contains the locations of some rows of the sparse matrix in
  ``matrix_name.mat``. For a directed graph, it outputs four files:
| ``matrix_name.mat``, ``matrix_name.mat_idx``, ``matrix_name_t.mat``
  and ``matrix_name_t.mat_idx``. The first two files contain the sparse
  matrix and the last two files contain the transpose of the sparse
  matrix.
