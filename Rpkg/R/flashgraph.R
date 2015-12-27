# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashR.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Reconfigure FlashR
#'
#' This reconfigures FlashR with the settings in the configuration file.
#' The configuration file contains a list of key-value pairs. Each line in
#' the file is a key-value pair in the form of "key_name=value".
#' @param conf.file The configuration file.
#' @name fg.set.conf
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.set.conf <- function(conf.file)
{
	ret <- .Call("R_FG_set_conf", conf.file, PACKAGE="FlashR")
	stopifnot(ret);
}

fg.set.log.level <- function(level)
{
	.Call("R_FG_set_log_level", level, PACKAGE="FlashR")
}

#' List graphs loaded to FlashR
#'
#' This function lists all graphs that have been loaded to FlashR.
#' @return A list of graphs in a data frame. The first column of the data
#' frame is the graph name. The second column indicates whether a graph
#' is stored in memory or on disks.
#' @name fg.list.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.list.graphs <- function()
{
	.Call("R_FG_list_graphs", PACKAGE="FlashR")
}

#' Indicate whether a graph has been loaded to FlashR
#'
#' This function indicates whether a graph has been loaded to FlashR.
#' @param graph A graph name.
#' @return A boolean value: true if the graph has been loaded to FlashR.
#' @name fg.exist.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.exist.graph <- function(graph)
{
	.Call("R_FG_exist_graph", graph, PACKAGE="FlashR")
}

fg.get.params <- function(name)
{
	.Call("R_FG_get_params", name, PACKAGE="FlashR")
}

#' Load a graph to FlashR.
#'
#' Load a graph to FlashR from difference sources.
#' 
#' `fg.load.graph' loads a graph from the following sources: an edge list
#' file in the text format and an adjacency list file in the FlashGraph format.
#'
#' `fg.load.igraph' loads a graph from an iGraph object.
#'
#' `fg.get.graph' gets a FlashGraph object that references a graph
#' that has been loaded to FlashR.
#'
#' Once a graph is loaded to FlashR, FlashR will maintain it.
#' A user can use fg.list.graphs() to list all graphs that have been loaded
#' FlashR and use fg.get.graph() to get a reference to a graph. A user
#' should provide a name for the graph so later on he or she will be able
#' to identify the graph more easily. By default, the graph name is
#' the input graph file name.
#' 
#' A graph in the FlashGraph format is stored in two files: a FlashGraph
#' adjacency list file and a FlashGraph index file. When a user provides
#' an index file, the input graph file is considered as an adjacency list
#' file, otherwise, an edge list file.
#'
#' When loading a graph from an edge list file, FlashR will construct
#' it into the FlashGraph format. A user needs to indicate whether the edge
#' list represents a directed graph. A user can also use multiple threads
#' to accelerate constructing a graph.
#'
#' When loading a graph from iGraph, FlashR
#' will construct it into the FlashGraph format. A user can use multiple
#' threads to accelerate graph construction.
#'
#' @param graph The input graph file or the input iGraph object.
#' @param index.file The input index file for the graph. A user only needs
#'                   to provide an index file if the input graph uses
#'                   the FlashGraph format.
#' @param graph.name The graph name a user provides when a graph is
#'                   loaded to FlashR.
#' @param directed   Indicate whether the input graph is directed. This is
#'                   only used if the input graph use the edge list format.
#' @param nthreads   The number of threads used to construct a graph to
#'                   the FlashGraph format.
#' @return a FlashGraph object.
#' @name fg.load.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("edge_list.txt")
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' ig <- read.graph("edge_list.txt")
#' fg <- fg.load.igraph(ig)
fg.load.graph <- function(graph, index.file = NULL, graph.name=graph,
						  directed=TRUE)
{
	if (is.null(index.file)) {
		ret <- .Call("R_FG_load_graph_el", graph.name, graph,
			  as.logical(directed), PACKAGE="FlashR")
		if (is.null(ret))
			ret
		else
			structure(ret, class="fg")
	}
	else {
		ret <- .Call("R_FG_load_graph_adj", graph.name, graph, index.file,
			  PACKAGE="FlashR")
		if (is.null(ret))
			ret
		else
			structure(ret, class="fg")
	}
}

#' @rdname fg.load.graph
fg.load.igraph <- function(graph, graph.name=paste("igraph-v", vcount(graph),
												  "-e", ecount(graph), sep = ""))
{
	stopifnot(is.igraph(graph))
	df <- get.data.frame(graph)
	# iGraph is 1-based but FlashGraph is 0-based, so we need to subtract
	# vertex IDs by 1.
	df["from"] <- df["from"] - 1
	df["to"] <- df["to"] - 1
	ret <- .Call("R_FG_load_graph_el_df", graph.name, df,
				 as.logical(is.directed(graph)), PACKAGE="FlashR")
	if (is.null(ret))
		ret
	else
		structure(ret, class="fg")
}

#' @rdname fg.load.graph 
fg.get.graph <- function(graph.name)
{
	stopifnot(fg.exist.graph(graph.name))
	ret <- .Call("R_FG_get_graph_obj", graph.name, PACKAGE="FlashR")
	if (is.null(ret))
		ret
	else
		structure(ret, class="fg")
}

#' Export graph image.
#'
#' This function exports a graph image in FlashR to the local filesystem.
#'
#' @param graph The FlashGraph object
#' @param graph.file The graph file in the local filesystem to which
#' the adjacency lists of the graph is exported to.
#' @param index.file The index file in the local filesystem to which
#' the index of the graph is exported to.
#' @return true if the graph is exported to the local filesystem correctly;
#' false, otherwise.
fg.export.graph <- function(graph, graph.file, index.file)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	.Call("R_FG_export_graph", graph, graph.file, index.file,
		  PACKAGE="FlashR")
}

#' Graph information
#'
#' Functions for providing the basic information of a graph.
#'
#' `fg.vcount' gets the number of vertices in a graph.
#'
#' `fg.ecount' gets the number of edges in a graph.
#'
#' `fg.in.mem' indicates whether a graph is stored in memory.
#'
#' `fg.is.directed' indicates whether a graph is directed.
#'
#' @param The FlashGraph object
#' @return `fg.vcount' and `fg.ecount' returns integer constants.
#' `fg.in.mem' and `fg.is.directed' returns boolean constants.
#' @name fg.graph.info
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fg.graph.info
fg.vcount <- function(graph)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	graph$vcount
}

#' @rdname fg.graph.info
fg.ecount <- function(graph)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	graph$ecount
}

#' @rdname fg.graph.info
fg.in.mem <- function(graph)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	graph$in.mem
}

#' @rdname fg.graph.info
fg.is.directed <- function(graph)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	graph$directed
}

#' Connected components of a graph
#'
#' Compute all (weakly or strongly) connected components of a graph.
#'
#' For an undirected graph, this function only computes connected components
#' and ignores the argument `mode'.
#'
#' For strongly connected components, we use a customized version of the
#' algorithm described in the paper
#'
#' Sungpack Hong, Nicole C. Rodia, Kunle Olukotun, On Fast Parallel Detection
#' of Strongly Connected Components (SCC) in Small-World Graphs, Proceedings
#' of the International Conference on High Performance Computing, Networking,
#' Storage and Analysis, 2013
#'
#' @param graph The FlashGraph object
#' @param mode Character string, either "weak" or "strong". For directed
#'             graphs "weak" implies weakly, "strong" strongly c
#'             components to search. It is ignored for undirected graphs.
#' @return A numeric vector that indicates the cluster id to which each vertex
#'         blongs to.
#' @name fg.cc
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @references
#' Sungpack Hong, Nicole C. Rodia, Kunle Olukotun, On Fast Parallel Detection
#' of Strongly Connected Components (SCC) in Small-World Graphs, Proceedings
#' of the International Conference on High Performance Computing, Networking,
#' Storage and Analysis, 2013
fg.clusters <- function(graph, mode=c("weak", "strong"))
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	if (!graph$directed)
		.Call("R_FG_compute_cc", graph, PACKAGE="FlashR")
	else if (mode == "weak")
		.Call("R_FG_compute_wcc", graph, PACKAGE="FlashR")
	else if (mode == "strong")
		.Call("R_FG_compute_scc", graph, PACKAGE="FlashR")
	else
		stop("a wrong mode")
}

#fg.transitivity <- function(graph)
#{
#	stopifnot(!is.null(graph))
#	stopifnot(class(graph) == "fg")
#	stopifnot(graph$directed)
#	.Call("R_FG_compute_transitivity", graph, PACKAGE="FlashR")
#}

#' Degree of the vertices in a graph
#'
#' Get the degree of vertices in a graph.
#'
#' @param graph The FlashGraph object
#' @param mode Character string. "out" for out-degree, "in" for in-degree,
#'        "both" for the sum of the two. This argument is ignored for
#'        undirected graphs.
#' @return A numeric vector with the degree of each vertex in the graph.
#' @name fg.degree
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.degree <- function(graph, mode=c("both", "in", "out"))
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	.Call("R_FG_get_degree", graph, mode, PACKAGE="FlashR")
}

#' PageRank
#'
#' Compute the Google PageRank for a graph.
#'
#' This implementation computes PageRank values in the original PageRank
#' paper below and does not normalize PageRank values in each iteration.
#'
#' Sergey Brin and Larry Page: The Anatomy of a Large-Scale
#' Hypertextual Web Search Engine. Proceedings of the 7th World-Wide
#' Web Conference, Brisbane, Australia, April 1998.
#'
#' To improve performance, a vertex only sends the difference of its PageRank
#' value between the previous iteration and the current iteration to its
#' neighbors in each iteration. If the difference is smaller than a threshold,
#' a vertex does not send the difference to its neighbors. The algorithm
#' converges if all vertices stop sending messages.
#'
#' @param graph The FlashGraph object
#' @param no.iters The number of iterations
#' @param damping The damping factor ('d' in the original p)
#' @return A numeric vector that contains PageRank values of each vertex.
#' @name fg.pagerank
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @references
#' Sergey Brin and Larry Page: The Anatomy of a Large-Scale
#' Hypertextual Web Search Engine. Proceedings of the 7th World-Wide
#' Web Conference, Brisbane, Australia, April 1998.
fg.page.rank <- function(graph, no.iters=1000, damping=0.85)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(graph$directed)
	.Call("R_FG_compute_pagerank", graph, no.iters, damping,
		  PACKAGE="FlashR")
}

#' Triangle counting
#'
#' Count the number of triangles on each vertex in a graph.
#'
#' Triangle counting works for both directed and undirected graphs.
#' For directed graphs, this counts the number of cycle triangles, shown
#' below.
#' A -> B
#' ^   /
#' | v
#' C
#' @param graph The FlashGraph object
#' @param type The type of triangles. It is ignored for undirected graphs.
#' @return A numeric vector that contains the number of triangles associated
#'         with each vertex.
#' @name fg.triangle
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.triangles <- function(graph, type="cycle")
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	if (graph$directed) {
		.Call("R_FG_compute_directed_triangles", graph, type,
			  PACKAGE="FlashR")
	}
	else {
		.Call("R_FG_compute_undirected_triangles", graph, PACKAGE="FlashR")
	}
}

#' Locality statistic
#'
#' Compute locality statistic of vertices in a graph.
#'
#' `fg.topK.scan' finds the top K vertices with the largest locality
#' statistics in a graph and computes their locality statistic.
#'
#' `fg.local.scan' computes locality statistic of each vertex in a graph.
#'
#' Locality statistic is defined as the number of edges in the neighborhood
#' of a vertex.
#'
#' @param graph The FlashGraph object
#' @param order An integer scalar, the size of the local neighborhood for
#'              each vertex. Should be non-negative.
#' @return A numeric vector that contains locality statistic of each vertex.
#' @name fg.local.scan
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @references
#' C.E. Priebe, J.M. Conroy, D.J. Marchette, and Y. Park, Scan Statistics on
#' Enron Graphs," Computational and Mathematical Organization Theory, 2005.

#' @rdname fg.local.scan
fg.topK.scan <- function(graph, order=1, K=1)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(graph$directed)
	.Call("R_FG_compute_topK_scan", graph, order, K, PACKAGE="FlashR")
}

#' @rdname fg.local.scan
fg.local.scan <- function(graph, order=1)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	if (graph$directed) {
		.Call("R_FG_compute_local_scan", graph, as.integer(order),
			  PACKAGE="FlashR")
	}
	else if (order == 0) {
		fg.degree(graph)
	}
	else if (order == 1) {
		fg.triangles(graph) + fg.degree(graph)
	}
	else if (order == 2) {
		print("We don't support local scan in the order of 2 on an undirected graph");
		NULL
	}
}

#' Transitivity of a graph
#'
#' Transitivity measures the probability that the adjacent vertices
#' of a vertex are connected. This is sometimes also called the clustering
#' coefficient.
#'
#' There are essentially two types of transitivity measures, one is
#' a vertex-level, the other a graph level property. This function can compute
#' both types of transitivity measures and works for both directed and
#' undirected graphs.
#' @param graph The FlashGraph object
#' @param type The type of the transitivity measure
#' @return A numeric vector that contains transitivity of each vertex if
#' `type' is "local" and a single value if `type' is "global".
#' @name fg.transitivity
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.transitivity <- function(graph, type=c("global", "local"))
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	deg <- fg.degree(graph)
	if (type == "local") {
		if (graph$directed) {
			(fg.local.scan(graph) - deg) / (deg * (deg - 1))
		}
		else {
			2 * fg.triangles(graph) / (deg * (deg - 1))
		}
	}
	else {
		if (graph$directed) {
			sum(fg.local.scan(graph) - deg) / sum(deg * (deg - 1))
		}
		else {
			2 * sum(fg.triangles(graph)) / sum(deg * (deg - 1))
		}
	}
}

#' K-core decomposition of a graph.
#'
#'  The k-core of graph is a maximal subgraph in which each vertex has
#'  at least degree k. A vertex belongs to the k-th core if has degree >= k
#'	when all connected vertices with degree < k are recursively deleted.
#'
#' @param graph The FlashGraph object
#' @param k.start The lowest core that should be computed. Must be >= 2.
#' @param k.end The highest core that should be computed. Must be >= 2. 
#'		default is 10.
#'
#' @return A numeric vector that contains the core of each
#'			vertex up to `k.end`. Vertices in cores higher than
#'		   `k.end` will have entries with `-1` as their core.
#'
#' @name fg.kcore
#' @author Disa Mhembere <disa@@jhu.edu>
#' @rdname fg.kcore
fg.kcore <- function(graph, k.start=2, k.end=10)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	.Call("R_FG_compute_kcore", graph, k.start, k.end, PACKAGE="FlashR")
}

fg.overlap <- function(graph, vids)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(graph$directed)
	.Call("R_FG_compute_overlap", graph, vids, PACKAGE="FlashR")
}

#' Generate an induced subgraph
#'
#' Generate an induced subgraph that contains the specified vertices from
#' the input given graph.
#'
#' `fg.fetch.subgraph.igraph' generates the induced subgraph and converts
#' it into an iGraph object
#'
#' `fg.fetch.subgraph' generates the induced subgraph and returns
#' a FlashGraph object.
#'
#' @param graph The FlashGraph object
#' @param vertices A numeric vector that contains the ids of vertices in
#'                 the induced subgraph.
#' @param name The name of the FlashGraph object.
#' @return `fg.fetch.subgraph.igraph' returns an iGraph object,
#' `fg.fetch.subgraph' returns a FlashGraph object
#' @name fg.subgraph
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fg.subgraph
fg.fetch.subgraph.igraph <- function(graph, vertices)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	edge.list = .Call("R_FG_fetch_subgraph_el", graph, vertices,
					  PACKAGE="FlashR")
	dframe = data.frame(edge.list$src, edge.list$dst)
	graph.data.frame(dframe, graph$directed)
}

#' @rdname fg.subgraph
fg.fetch.subgraph <- function(graph, vertices,
							  name = paste(graph$name, "-sub", sep=""), compress=TRUE)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	ret <- .Call("R_FG_fetch_subgraph", graph, vertices, name, compress,
				 PACKAGE="FlashR")
	if (is.null(ret))
		ret
	else
		structure(ret, class="fg")
}

#' Diameter estimation
#'
#' Estimate the diameter of a graph, the longest distance between two vertices
#' in a graph.
#'
#' This implementation uses the double sweep method described in the paper
#' below to estimate the lower bound of the diameter.
#'
#' Cl茅mence Magnienand Matthieu Latapy and Michel Habib: Fast computation of
#' empirically tight bounds for the diameter of massive graphs, Journal of
#' Experimental Algorithmics (JEA), 2009
#'
#' @param graph The FlashGraph object
#' @param directed Indicates whether or not to respect the direction of edges
#'                 in a graph when traversing the graph. It is ignored for
#'                 undirected graphs.
#' @return A single number
#' @name fg.diameter
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @references
#' Cl茅mence Magnienand Matthieu Latapy and Michel Habib: Fast computation of
#' empirically tight bounds for the diameter of massive graphs, Journal of
#' Experimental Algorithmics (JEA), 2009
fg.diameter <- function(graph, directed=FALSE)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(graph$directed)
	.Call("R_FG_estimate_diameter", graph, directed, PACKAGE="FlashR")
}

#' Sparse matrix multiplication
#'
#' Multiply a sparse matrix with a dense vector or a dense matrix.
#'
#' `fg.multiply' multiplies a sparse matrix with a dense vector.
#'
#' `fg.multiply.matrix' multiplies a sparse matrix with a dense matrix.
#'
#' Note that the sparse matrix is represented by a graph. A symmetric matrix
#' is represented by an undirected graph and an asymmetric matrix is
#' represented by a directed graph. For an asymmetric matrix, we can multiply
#' the vector with the transpose of the matrix without really transposing
#' the matrix. Right now it only supports multiplication on the adjacency
#' matrix of the graph.
#'
#' @param graph The FlashGraph object
#' @param vec A numueric vector
#' @param m A numeric dense matrix.
#' @param transpose Indicate whether or not to multiply with the transpose of
#'                  the matrix. It is ignored by an undirected graph.
#' @return `fg.multiply' returns a numeric vector and `fg.multiply.matrix'
#' returns a numeric dense matrix.
#' @name fg.multiply
#' @author Da Zheng <dzheng5@@jhu.edu>

#' @rdname fg.multiply
fg.multiply <- function(graph, vec, transpose=FALSE)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(fg.vcount(graph) == length(vec))
#	stopifnot(graph$directed)
	.Call("R_FG_multiply_v", graph, vec, transpose, PACKAGE="FlashR")
}

#' @rdname fg.multiply
fg.multiply.matrix <- function(graph, m, transpose=FALSE)
{
	col.multiply <- function(x) {
		fg.multiply(graph, x, transpose)
	}
	apply(m, 2, col.multiply)
}

#' Perform spectral clustering
#'
#' Spectral clustering partitions an undirected graph into k groups based on
#' an embedding of the graph, which computes eigenvectors from a variant
#' of the adjacency matrix of the graph. Spectral clustering runs
#' KMeans clustering on eigenvectors.
#'
#' @param fg The FlashGraph object of the undirected graph.
#' @param k The number of clusters.
#' @param ase The function computes spectral embedding and generates
#'            eigenvectors.
#' @return A vector of integers (from '1:k') indicating the cluster to which
#'         each vertex belongs to.
#' @examples
#' fg <- fg.load.graph("edge_list.txt")
#' res <- fg.spectral.clusters(fg, 10)
#' @name fg.spectral
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.spectral.clusters <- function(fg, k, ase)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	stopifnot(!fg$directed)

	eigen <- ase(fg, k)
	kmeans(eigen$vectors, k, MacQueen)
}

print.fg <- function(fg)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	directed = "U"
	if (fg.is.directed(fg))
		directed = "D"
	cat("FlashGraph ", fg$name, " (", directed, "): ", fg.vcount(fg), " ", fg.ecount(fg),
		"\n", sep="")
}

#' Perform k-means clustering on a data matrix.
#'
#' Assign each row of a matrix to a cluster denoted by a numeric. The clusters
#' are based on the euclidean distance of each row to one another. The assigned
#' cluster will have the smallest distance from the cluster center mean.
#'
#' @param mat A numeric matrix of data.
#' @param k The number of clusters.
#' @param max.iters The maximum number of iterations allowed.
#' @param max.threads The maximum number of threads allowed (Default is 1 per core).
#' @param init The form of initialization to use when the algorithm begins.
#'              The default is "random". For a desciption of each see:
#'              http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods.
#' @param tol The tolerance for convergence. Between 0 and 1 and is the minimum fraction
#'				of cluster changes necessary to cause non-convergence. Default is -1 which
#'				represents no cluster changes. 
#'
#' @return A named list with the following members:
#'         iters: The number of (outer) iterations performed.
#'         centers: A matrix of cluster centers.
#'         cluster: A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
#'         sizes: The number of points in each cluster.
#'
#' @examples
#' num.clusts <- 3
#' mat <- replicate(5, rnorm(10))
#' kms <- fg.kmeans(mat, num.clusts)
#'
#' @name fg.kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
fg.kmeans <- function(mat, k, max.iters=10, max.threads=256,
					  init=c("random", "forgy","kmeanspp"), tol=-1)
{
    stopifnot(mat != NULL)
    stopifnot(class(mat) == "matrix")
	stopifnot(as.integer(max.threads) > 0)
    .Call("R_FG_kmeans", as.matrix(mat), as.integer(k), as.integer(max.iters),
		  as.integer(max.threads), init, as.double(tol), PACKAGE="FlashR")
}

#' Vertex betweenness centrality.
#'
#' The vertex betweenness centrality can be defined as the
#' number of geodesics (shortest paths) going through a vertex.
#'
#' @param fg The FlashGraph object.
#' @param vids A vector of vertex IDs. Default runs it on the entire
#'		graph.
#'
#' @return A vector with betweenness centrality values for all vertices
#'			with respect to `vids`.
#'
#' @examples
#' fg <- fg.load.graph("edge_list.txt")
#' res <- fg.betweenness(fg, c(1,10))
#'
#' @name fg.betweenness
#' @author Disa Mhembere <disa@@jhu.edu>
fg.betweenness <- function(fg, vids=0:(fg$vcount-1))
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
    .Call("R_FG_compute_betweenness", fg, vids, PACKAGE="FlashR")
}

.onLoad <- function(libname, pkgname)
{
	library(Rcpp)
	library.dynam("FlashR", pkgname, libname, local=FALSE);
	ret <- .Call("R_FG_init", paste(pkgname, ".conf", sep=""), PACKAGE="FlashR")
	stopifnot(ret)
	fm.init.basic.op()
}

#' Perform semi-extenal memory k-means clustering on a data matrix in FG format.
#'
#' Assign each row of a matrix to a cluster denoted by a numeric. The clusters
#' are based on the distance metric of each row to one another. The assigned
#' cluster will have the smallest distance from the cluster center mean.
#'
#' @param mat A numeric matrix of data in FG format.
#' @param k The number of clusters.
#' @param max.iters The maximum number of iterations allowed.
#' @param init The form of initialization to use when the algorithm begins.
#'              The default is "forgy". For a desciption of each see:
#'              http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods.
#' @param tol The tolerance for convergence. Between 0 and 1 and is the minimum fraction
#'				of cluster changes necessary to cause non-convergence. Default is -1 which
#'				represents no cluster changes.
#' @return A named list with the following members:
#'         iters: The number of (outer) iterations performed.
#'         centers: A matrix of cluster centers.
#'         cluster: A vector of integers (from 1:k) indicating the cluster to which each point is allocated.
#'         sizes: The number of points in each cluster.
#'
#' @examples
#' num.clusts <- 5
#' mat <- fg.load.graph("mat-adj", "mat-index") # TODO: how to load
#' kms <- fg.kmeans(mat, num.clusts)
#'
#' @name fg.sem.kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
fg.sem.kmeans <- function(mat, k, max.iters=10, init=c("random", "forgy","kmeanspp"),
                          tol=-1)
{
    stopifnot(mat != NULL)
    stopifnot(class(mat) == "fg")
    .Call("R_FG_sem_kmeans", mat, as.integer(k), init, as.integer(max.iters),
		  as.double(tol), PACKAGE="FlashR")
}
