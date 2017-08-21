This repo contains the core of the FlashX project, which provides big data analytics tools
that perform data analytics in the form of graphs and matrices. As such, FlashX covers
a large range of data analysis tasks. All tools in FlashX utilize solid-state drives (SSDs) to
scale data analysis to large datasets in a single machine, while achieving
lightning speed (SSD-based solutions run almost as fast as in-memory solutions).
The main components in FlashX are FlashGraph and FlashMatrix.

FlashGraph
===========

FlashGraph is a general-purpose graph analysis framework that exposes
vertex-centric programming interface for users to express varieties of
graph algorithms. FlashGraph scales graph computation to large graphs by
keeping the edges of a graph on SSDs and computation state in memory.
With smart I/O scheduling, FlashGraph is able to achieve performance
comparable to state-of-art in-memory graph analysis frameworks and
significantly outperforms state-of-art distributed graph analysis frameworks
while being able to scale to graphs with billions of vertices and hundreds
of billions of edges. Please see
[the performance result](https://flashxio.github.io/FlashX-doc/FlashX-perf.html#flashgraph-vs-giraph-graphx-and-powergraph).

FlashMatrix
===========

FlashMatrix is a matrix computation engine that provides a small set of
generalized matrix operations on sparse matrices and dense matrices to express
varieties of data mining and machine learning algorithms. For certain graph
algorithms such as PageRank, which can be formulated as sparse matrix
multiplication, FlashMatrix is able to significantly outperform FlashGraph.

Programming interface
===========
FlashX exposes C++, R and Python programming interface. The R and Python programming interface
is highly compatible with the R base package and NumPy. As such, users can execute
R and Python machine learning code on FlashX with little or no modification. Our goal is to
eventually make the R and Python interface fully compatible with the ones in native R and NumPy.

* [FlashR](https://github.com/flashxio/FlashR) provides many matrix operations in the R base package.
* [FlashGraphR](https://github.com/flashxio/FlashGraphR) exposes many graph algorithms in FlashGraph to R.
* [FlashR-learn](https://github.com/flashxio/FlashR-learn) is a machine learning library implemented completely with FlashR.
* [FlashPy](https://github.com/flashxio/FlashPy) provides many array operations in NumPy.

[Documentation](https://flashxio.github.io/FlashX-doc/)
========

[FlashX Quick start guide](https://flashxio.github.io/FlashX-doc/FlashX-Quick-Start-Guide.html)

[FlashGraph programming tutorial](https://flashxio.github.io/FlashX-doc/FlashGraph-user-guide.html).

[FlashR programming tutorial](https://flashxio.github.io/FlashX-doc/FlashR-user-guide.html)

[FlashX performance and scalability](https://flashxio.github.io/FlashX-doc/FlashX-perf.html)

Publications
========
Da Zheng, Disa Mhembere, Joshua T. Vogelstein, Carey E. Priebe, and Randal Burns, “FlashMatrix: Parallel, scalable data analysis with generalized matrix operations using commodity ssds,” arXiv preprint arXiv:1604.06414, 2016 [[pdf](http://arxiv.org/pdf/1604.06414v3)]

Da Zheng, Disa Mhembere, Vince Lyzinski, Joshua Vogelstein, Carey E. Priebe, and Randal Burns, “Semi-external memory sparse matrix multiplication on billion-node graphs”, Transactions on Parallel and Distributed Systems, 2016. [[pdf](https://arxiv.org/pdf/1602.02864v3.pdf)]

Heng Wang, Da Zheng, Randal Burns, Carey Priebe, Active Community Detection in Massive Graphs, SDM-Networks 2015 [[pdf](http://arxiv.org/pdf/1412.8576v3.pdf)]

Da Zheng, Disa Mhembere, Randal Burns, Joshua Vogelstein, Carey E. Priebe, Alexander S. Szalay, FlashGraph: Processing Billion-Node Graphs on an Array of Commodity SSDs, FAST'15, [[pdf](https://www.usenix.org/system/files/conference/fast15/fast15-paper-zheng.pdf)][[bib](https://www.usenix.org/biblio/export/bibtex/188418)]

Da Zheng, Randal Burns, Alexander S. Szalay, Toward Millions of File System IOPS on Low-Cost, Commodity Hardware, in Proceeding SC '13 Proceedings of the International Conference on High Performance Computing, Networking, Storage and Analysis, [[pdf](http://www.cs.jhu.edu/~zhengda/sc13.pdf)][[bib](http://dl.acm.org/downformats.cfm?id=2503225&parent_id=2503210&expformat=bibtex&CFID=445591569&CFTOKEN=95321450)]

Contact
========

Mailing list: flashgraph-user@googlegroups.com

[![Join the chat at https://gitter.im/icoming/FlashGraph](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/icoming/FlashGraph?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
