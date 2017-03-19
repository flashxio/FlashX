# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashGraphR.
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

#' Reconfigure FlashGraphR
#'
#' This reconfigures FlashGraphR with the settings in the configuration file.
#' The configuration file contains a list of key-value pairs. Each line in
#' the file is a key-value pair in the form of "key_name=value".
#' @param conf.file The configuration file.
#' @name fg.set.conf
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.set.conf <- function(conf.file)
{
	.Call("R_FG_destroy", PACKAGE="FlashGraphR")
	fm.set.conf(conf.file)
	ret <- .Call("R_FG_init", conf.file, PACKAGE="FlashGraphR")
}

fg.set.log.level <- function(level)
{
	.Call("R_FG_set_log_level", level, PACKAGE="FlashGraphR")
}

#' List graphs loaded to FlashGraphR
#'
#' This function lists all graphs that have been loaded to FlashGraphR.
#' @return A list of graphs in a data frame. The first column of the data
#' frame is the graph name. The second column indicates whether a graph
#' is stored in memory or on disks.
#' @name fg.list.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.list.graphs <- function()
{
	.Call("R_FG_list_graphs", PACKAGE="FlashGraphR")
}

#' Indicate whether a graph has been loaded to FlashGraphR
#'
#' This function indicates whether a graph has been loaded to FlashGraphR.
#' @param graph A graph name.
#' @return A boolean value: true if the graph has been loaded to FlashGraphR.
#' @name fg.exist.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.exist.graph <- function(graph)
{
	.Call("R_FG_exist_graph", graph, PACKAGE="FlashGraphR")
}

fg.get.params <- function(name)
{
	.Call("R_FG_get_params", name, PACKAGE="FlashGraphR")
}

#' Load a graph to FlashGraphR.
#'
#' Load a graph to FlashGraphR from difference sources.
#' 
#' `fg.load.graph' loads a graph from the following sources: an edge list
#' file in the text format and an adjacency list file in the FlashGraph format.
#'
#' `fg.load.igraph' loads a graph from an iGraph object.
#'
#' `fg.get.graph' gets a FlashGraph object that references a graph
#' that has been loaded to FlashGraphR.
#'
#' Once a graph is loaded to FlashGraphR, FlashGraphR will maintain it.
#' A user can use fg.list.graphs() to list all graphs that have been loaded
#' FlashGraphR and use fg.get.graph() to get a reference to a graph. A user
#' should provide a name for the graph so later on he or she will be able
#' to identify the graph more easily. By default, the graph name is
#' the input graph file name.
#' 
#' A graph in the FlashGraph format is stored in two files: a FlashGraph
#' adjacency list file and a FlashGraph index file. When a user provides
#' an index file, the input graph file is considered as an adjacency list
#' file, otherwise, an edge list file.
#'
#' When loading a graph from an edge list file, FlashGraphR will construct
#' it into the FlashGraph format. A user needs to indicate whether the edge
#' list represents a directed graph. A user can also use multiple threads
#' to accelerate constructing a graph.
#'
#' When loading a graph from iGraph, FlashGraphR
#' will construct it into the FlashGraph format. A user can use multiple
#' threads to accelerate graph construction.
#'
#' @param graph The input graph file or the input iGraph object.
#' @param index.file The input index file for the graph. A user only needs
#'                   to provide an index file if the input graph uses
#'                   the FlashGraph format.
#' @param graph.name The graph name a user provides when a graph is
#'                   loaded to FlashGraphR.
#' @param directed   Indicate whether the input graph is directed. This is
#'                   only used if the input graph use the edge list format.
#' @param in.mem     Indicate whether to load a graph to SAFS.
#' @param delim		 The delimiter of separating elements in the text format.
#'					 When delim is "auto", FlashGraph will try to detect
#'					 the delimiter automatically.
#' @return a FlashGraph object.
#' @name fg.load.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("edge_list.txt")
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' ig <- read.graph("edge_list.txt")
#' fg <- fg.load.igraph(ig)
fg.load.graph <- function(graph, index.file = NULL, graph.name=graph,
						  directed=TRUE, in.mem=TRUE, delim="auto", attr.type="")
{
	# The graph name will becomes the file name in SAFS. It should contain
	# some special characters.
	graph.name <- gsub("/", "_", graph.name)
	graph.name <- gsub(" ", "_", graph.name)
	if (is.null(index.file)) {
		ret <- .Call("R_FG_load_graph_el", graph.name, graph,
			  as.logical(directed), as.logical(in.mem), as.character(delim),
			  as.character(attr.type), PACKAGE="FlashGraphR")
		if (is.null(ret))
			ret
		else
			structure(ret, class="fg")
	}
	else {
		ret <- .Call("R_FG_load_graph_adj", graph.name, graph, index.file,
			  PACKAGE="FlashGraphR")
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
	# The graph name will becomes the file name in SAFS. It should contain
	# some special characters.
	graph.name <- gsub("/", "_", graph.name)
	graph.name <- gsub(" ", "_", graph.name)

	stopifnot(is.igraph(graph))
	df <- get.data.frame(graph)
	# iGraph is 1-based but FlashGraph is 0-based, so we need to subtract
	# vertex IDs by 1.
	df["from"] <- df["from"] - 1
	df["to"] <- df["to"] - 1
	ret <- .Call("R_FG_load_graph_el_df", graph.name, df,
				 as.logical(is.directed(graph)), PACKAGE="FlashGraphR")
	if (is.null(ret))
		ret
	else
		structure(ret, class="fg")
}

#' @rdname fg.load.graph 
fg.get.graph <- function(graph.name)
{
	stopifnot(fg.exist.graph(graph.name))
	ret <- .Call("R_FG_get_graph_obj", graph.name, PACKAGE="FlashGraphR")
	if (is.null(ret))
		ret
	else
		structure(ret, class="fg")
}

#' Export graph image.
#'
#' This function exports a graph image in FlashGraphR to the local filesystem.
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
		  PACKAGE="FlashGraphR")
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
#' @param graph a FlashGraph object
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
		ret <- .Call("R_FG_compute_cc", graph, PACKAGE="FlashGraphR")
	else if (mode == "weak")
		ret <- .Call("R_FG_compute_wcc", graph, PACKAGE="FlashGraphR")
	else if (mode == "strong")
		ret <- .Call("R_FG_compute_scc", graph, PACKAGE="FlashGraphR")
	else
		stop("a wrong mode")
	new_fmV(ret)
}

#' Get the largest connected component in a graph
#'
#' Get the largest (weakly or strongly) connected component in a graph.
#'
#' For an undirected graph, this function returns the largest connected component
#' and ignores the argument `mode'.
#'
#' @param graph The FlashGraph object
#' @param mode Character string, either "weak" or "strong". For directed
#'             graphs "weak" implies weakly, "strong" strongly c
#'             components to search. It is ignored for undirected graphs.
#' @return a FlashGraph object that contains the largest connected component in
#'         the graph.
#' @name fg.lcc
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.get.lcc <- function(graph, mode=c("weak", "strong"))
{
	cc <- fg.clusters(graph, mode)
	res <- table(cc)
	df <- as.data.frame(res)
	lcc.id <- df$cc[which.max(df$Freq)]
	lccV <- which(cc == lcc.id)
	fg.fetch.subgraph(graph, vertices=lccV - 1, compress=TRUE)
}

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
	ret <- .Call("R_FG_get_degree", graph, mode, PACKAGE="FlashGraphR")
	new_fmV(ret)
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
	ret <- .Call("R_FG_compute_pagerank", graph, no.iters, damping,
				 PACKAGE="FlashGraphR")
	new_fmV(ret)
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
		ret <- .Call("R_FG_compute_directed_triangles", graph, type,
					 PACKAGE="FlashGraphR")
	}
	else {
		ret <- .Call("R_FG_compute_undirected_triangles", graph,
					 PACKAGE="FlashGraphR")
	}
	new_fmV(ret)
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
#' @param K The number of top scan statistics.
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
	.Call("R_FG_compute_topK_scan", graph, order, K, PACKAGE="FlashGraphR")
}

#' @rdname fg.local.scan
fg.local.scan <- function(graph, order=1)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	if (graph$directed) {
		ret <- .Call("R_FG_compute_local_scan", graph, as.integer(order),
					 PACKAGE="FlashGraphR")
		new_fmV(ret)
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
	ret <- .Call("R_FG_compute_kcore", graph, k.start, k.end,
				 PACKAGE="FlashGraphR")
	new_fmV(ret)
}

fg.overlap <- function(graph, vids)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(graph$directed)
	.Call("R_FG_compute_overlap", graph, vids, PACKAGE="FlashGraphR")
}

#' Generate an induced subgraph
#'
#' Generate an induced subgraph that contains the specified vertices from
#' the input given graph.
#'
#' `fg.fetch.subgraph' generates the induced subgraph and returns
#' a FlashGraph object.
#'
#' @param graph The FlashGraph object
#' @param vertices A numeric vector that contains the ids of vertices in
#'                 the induced subgraph.
#' @param compress This indicates whether to remove empty vertices in
#'                 the generated subgraph.
#' @param name The name of the FlashGraph object.
#' @return `fg.fetch.subgraph' returns a FlashGraph object
#' @name fg.subgraph
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.fetch.subgraph <- function(graph, vertices,
							  name = paste(graph$name, "-sub", sep=""), compress=TRUE)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	ret <- .Call("R_FG_fetch_subgraph", graph, vertices, name, compress,
				 PACKAGE="FlashGraphR")
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
	.Call("R_FG_estimate_diameter", graph, directed, PACKAGE="FlashGraphR")
}

print.fg <- function(x, ...)
{
	stopifnot(!is.null(x))
	stopifnot(class(x) == "fg")
	directed = "U"
	if (fg.is.directed(x))
		directed = "D"
	cat("FlashGraph ", x$name, " (", directed, "): ", fg.vcount(x), " ", fg.ecount(x),
		"\n", sep="")
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
    ret <- .Call("R_FG_compute_betweenness", fg, vids, PACKAGE="FlashGraphR")
	new_fmV(ret)
}

.onLoad <- function(libname, pkgname)
{
	library.dynam("FlashGraphR", pkgname, libname, local=FALSE);
	ret <- .Call("R_FG_init", NULL, PACKAGE="FlashGraphR")
	stopifnot(ret)
}

.new.fm <- function(fm)
{
	if (!is.null(fm))
		new("fm", pointer=fm$pointer, name=fm$name, nrow=fm$nrow, ncol=fm$ncol,
			type=fm$type, ele_type=fm$ele_type)
	else
		NULL
}

#' Load a sparse matrix to FlashGraphR.
#'
#' Load a sparse matrix to FlashGraphR from FlashGraph.
#' 
#' `fg.get.sparse.matrix' gets a FlashMatrix sparse matrix that references
#' a matrix represented by a FlashGraph object.
#'
#' @examples
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' fm <- fg.get.sparse.matrix(fg)
fg.get.sparse.matrix <- function(fg)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	stopifnot(fg.exist.graph(fg$name))
	m <- .Call("R_FG_get_matrix_fg", fg, PACKAGE="FlashGraphR")
	new_fm(m)
}

#' Print a graph into a file as an edge list.
#'
#' Print a graph in the FlashGraph format into a file as an edge list.
#'
#' @param fg the FlashGraph object
#' @param file a string of the output file name.
fg.print.graph <- function(fg, file, delim="\t", type="")
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	.Call("R_FG_print_graph", fg, as.character(file), as.character(delim),
		  as.character(type), PACKAGE="FlashGraphR")
}
