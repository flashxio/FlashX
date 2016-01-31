FlashGraph performance
======================

FlashGraph on a multi-million-node graph
----------------------------------

| FlashGraph was designed for specialized hardware (a large parallel
  machine with many SSDs).
  However, most people may not possess such powerful hardware, so we
  also evaluate its performance in the Amazon cloud. Table 1 shows the
  hardware configurations where FlashGraph is evaluated. We demonstrate
  that FlashGraph is able to outperform other well-known graph
  processing frameworks in both specialized hardware and Amazon cloud.

+------------+------------+------------+-----+-----+
|            | i2.xlarge  | i2.8xlarge | HW1 | HW2 |
+============+============+============+=====+=====+
| #CPU cores | 4          | 32         | 32  | 48  |
+------------+------------+------------+-----+-----+
| RAM (GB)   | 30         | 244        | 512 | 1024|
+------------+------------+------------+-----+-----+
| #SSDs      | 1          | 8          | 15  | 24  |
+------------+------------+------------+-----+-----+

| Table 1. The hardware configuration where FlashGraph is evaluated.

+----------+-----------+------------+-------+-------+
|          | i2.xlarge | i2.8xlarge | HW1   | HW2   |
+==========+===========+============+=======+=======+
| BFS      | 36.86     | 5.48       | 3.56  | 3.11  |
+----------+-----------+------------+-------+-------+
| PageRank | 438.74    | 67.34      | 54.07 | 33.58 |
+----------+-----------+------------+-------+-------+
| WCC      | 58.90     | 7.91       | 7.91  | 5.03  |
+----------+-----------+------------+-------+-------+
| Triangle | 4747.07   | 830.08     | 532.21| 442.65|
+----------+-----------+------------+-------+-------+
| SCC      | 141.62    | 22.39      | N/A   | 12.06 |
+----------+-----------+------------+-------+-------+

| Table 2. The runtime (seconds) of FlashGraph on different hardware
  with 1GB page cache.

| We evaluate the performance of FlashGraph on the Twitter graph with 42
  millions vertices and 1.5 billion edges (Table 2). The reason that we
  choose the Twitter graph is that this graph is used by many graph
  processing frameworks for performance evaluation. We run multiple
  graph algorithms: breadth-first search, triangle counting, weakly
  connected components and pagerank.
| The twitter graph is relatively small and all of the hardware has
  enough RAM to keep the entire graph in memory. We artificially limit
  the page cache size in FlashGraph to 1GB and 4GB to push most of data
  access to SSDs. As a result, FlashGraph only uses a small fraction of
  the RAM in the machines. Keep in mind that FlashGraph runs much faster
  if we do not limit its page cache size. To our surprise, the
  i2.8xlarge instance has performance very close to our specialized
  hardware.

For comparison, we pull performance results of distributed graph engines
from the `GraphX
paper <https://amplab.cs.berkeley.edu/wp-content/uploads/2014/09/graphx.pdf>`__.
The paper don't have performance results of all of the graph
applications we evaluated. Therefore, we only compare their performance
in PageRank and weakly connected components.
`Giraph <http://giraph.apache.org/>`__,
`GraphX <http://spark.apache.org/graphx/>`__ and
`PowerGraph <https://github.com/dato-code/PowerGraph>`__ were evaluated
in 16 m2.4xlarge instances, which each has 8 CPU cores and 68GB RAM.
Only PowerGraph is a little faster than FlashGraph in PageRank when
FlashGraph runs in a small Amazon instance.

.. image:: http://icoming.github.io/FlashX/tutorials/FlashGraph-runtime.png

Figure 1. The runtime of FlashGraph vs. distributed graph engines in
PageRank and weakly connected components.

Given such performance results, we can further demonstrate that
FlashGraph is much more economical than distributed graph engines in the
Amazon cloud. m2.4xlarge is the previous-generation instance type. A
similar instance (m4.2xlarge) of the current generation has cut the
price almost by half and we use the cost of m4.2xlarge for the
calculation. Figure 2 shows the runtime dollars (=runtime \* instance
cost) of the graph engines in the cloud. When comparing these graph
engines with this metric, FlashGraph in i2.xlarge is the most
economical. These distributed graph engines are one order of magnitude
more costly than FlashGraph in the cloud.

.. image:: http://icoming.github.io/FlashX/tutorials/FlashGraph-runtime-dollors.png

Figure 2. The runtime dollars of FlashGraph vs. distributed graph
engines.

FlashGraph on a billion-node graph
----------------------------------

FlashGraph enables us to process very large graphs in a single machine.
Here we demonstrate that FlashGraph can process a real-world `hyperlink
page graph <http://webdatacommons.org/hyperlinkgraph/>`__ (3.4B vertices
and 129B edges) in a single machine. The table below shows the
performance of the same graph applications as above. It also includes
the performance of Betweenness centrality (BC) from a single vertex.

+-------------+---------------+-------------+
|             | Runtime (sec) | Memory (GB) |
+=============+===============+=============+
| BFS         | 298           | 22          |
+-------------+---------------+-------------+
| Betweenness | 595           | 81          |
+-------------+---------------+-------------+
| Triangle    | 7818          | 55          |
+-------------+---------------+-------------+
| WCC         | 461           | 47          |
+-------------+---------------+-------------+
| PageRank    | 2041          | 46          |
+-------------+---------------+-------------+

| Table 3. The runtime of FlashGraph on the Page graph on HW1.

FlashGraph takes only six minutes to traverse the page graph (3.5B
vertices and 128B edges) with a cache size of 4GB.
`Pregel <http://dl.acm.org/citation.cfm?id=1807184>`__ used 300
multicore machines to run the shortest path algorithm on their largest
random graph (1B vertices and 127B edges) and took a little over ten
minutes. More recently,
`Trinity <http://research.microsoft.com/en-us/projects/trinity/>`__ took
over ten minutes to perform breadth-first search on a graph of one
billion vertices and 13 billion edges on 14 12-core machines.

FlashGraph takes 34 minutes to run 30 iterations of PageRank on the page
graph. Based on a recent talk by the main developer of
`Giraph <http://www.youtube.com/watch?v=b5Qmz4zPj-M>`__, Giraph running
with 20 workers can run five iterations of PageRank on a graph with 5
billion edges in approx. 75 minutes.

.. |runtime| image:: http://flashx.io/images/FlashGraph.vs.others.png
.. |runtime dollors| image:: http://flashx.io/images/FlashGraph.vs.others.dollor.png
