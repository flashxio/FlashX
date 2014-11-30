#package.skeleton('flashgraph', code_files = 'flashgraph.R')

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
	ret <- .Call("R_FG_set_conf", conf.file, PACKAGE="FlashGraphR")
	stopifnot(ret);
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
#' @param nthreads   The number of threads used to construct a graph to
#'                   the FlashGraph format.
#' @return a FlashGraphR object.
#' @name fg.load.graph
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @examples
#' fg <- fg.load.graph("edge_list.txt")
#' fg <- fg.load.graph("graph.adj", "graph.index")
#' ig <- read.graph("edge_list.txt")
#' fg <- fg.load.igraph(ig)
fg.load.graph <- function(graph, index.file = NULL, graph.name=graph,
						  directed=TRUE, nthreads=1)
{
	if (is.null(index.file)) {
		ret <- .Call("R_FG_load_graph_el", graph.name, graph,
			  as.logical(directed), as.integer(nthreads), PACKAGE="FlashGraphR")
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
												  "-e", ecount(graph), sep = ""),
						  nthreads=1)
{
	stopifnot(is.igraph(graph))
	df <- get.data.frame(graph)
	# iGraph is 1-based but FlashGraph is 0-based, so we need to subtract
	# vertex IDs by 1.
	df["from"] <- df["from"] - 1
	df["to"] <- df["to"] - 1
	ret <- .Call("R_FG_load_graph_el_df", graph.name, df,
				 as.logical(is.directed(graph)), as.integer(nthreads),
				 PACKAGE="FlashGraphR")
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
#' @param The FlashGraphR object
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
#' @param graph The FlashGraphR object
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
		.Call("R_FG_compute_cc", graph, PACKAGE="FlashGraphR")
	else if (mode == "weak")
		.Call("R_FG_compute_wcc", graph, PACKAGE="FlashGraphR")
	else if (mode == "strong")
		.Call("R_FG_compute_scc", graph, PACKAGE="FlashGraphR")
	else
		stop("a wrong mode")
}

#fg.transitivity <- function(graph)
#{
#	stopifnot(!is.null(graph))
#	stopifnot(class(graph) == "fg")
#	stopifnot(graph$directed)
#	.Call("R_FG_compute_transitivity", graph, PACKAGE="FlashGraphR")
#}

#' Degree of the vertices in a graph
#'
#' Get the degree of vertices in a graph.
#'
#' @param graph The FlashGraphR object
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
	.Call("R_FG_get_degree", graph, mode, PACKAGE="FlashGraphR")
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
#' @param graph The FlashGraphR object
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
		  PACKAGE="FlashGraphR")
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
#' @param graph The FlashGraphR object
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
			  PACKAGE="FlashGraphR")
	}
	else {
		.Call("R_FG_compute_undirected_triangles", graph, PACKAGE="FlashGraphR")
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
#' @param graph The FlashGraphR object
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
	.Call("R_FG_compute_topK_scan", graph, order, K, PACKAGE="FlashGraphR")
}

#' @rdname fg.local.scan
fg.local.scan <- function(graph, order=1)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	if (graph$directed) {
		.Call("R_FG_compute_local_scan", graph, order, PACKAGE="FlashGraphR")
	}
	else {
		fg.triangles(graph) + fg.degree(graph)
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
#' @param graph The FlashGraphR object
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
#' @param graph The FlashGraphR object
#' @param k.start The lowest core that should be computed. Must be >= 2.
#' @param k.end The highest core that should be computed. Must be >= 2. 
#'		default is 10.
#' @return A numeric vector that contains the core of each
#		vertex up to `k.end`. Vertices in cores higher than
#		`k.end` will have entries with `-1` as their core.
#' @name fg.kcore
#' @author Disa Mhembere <disa@@jhu.edu>
#' @rdname fg.kcore
fg.kcore <- function(graph, k.start=2, k.end=10)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	.Call("R_FG_compute_kcore", graph, k.start, k.end, PACKAGE="FlashGraphR")
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
#' `fg.fetch.subgraph.igraph' generates the induced subgraph and converts
#' it into an iGraph object
#'
#' `fg.fetch.subgraph' generates the induced subgraph and returns
#' a FlashGraphR object.
#'
#' @param graph The FlashGraphR object
#' @param vertices A numeric vector that contains the ids of vertices in
#'                 the induced subgraph.
#' @param name The name of the FlashGraphR object.
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
					  PACKAGE="FlashGraphR")
	dframe = data.frame(edge.list$src, edge.list$dst)
	graph.data.frame(dframe, graph$directed)
}

#' @rdname fg.subgraph
fg.fetch.subgraph <- function(graph, vertices, name = paste(graph$name, "-sub", sep=""))
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	ret <- .Call("R_FG_fetch_subgraph", graph, vertices, name, PACKAGE="FlashGraphR")
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
#' @param graph The FlashGraphR object
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
#' @param graph The FlashGraphR object
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
	.Call("R_FG_multiply_v", graph, vec, transpose, PACKAGE="FlashGraphR")
}

#' @rdname fg.multiply
fg.multiply.matrix <- function(graph, m, transpose=FALSE)
{
	col.multiply <- function(x) {
		fg.multiply(graph, x, transpose)
	}
	apply(m, 2, col.multiply)
}

#' Eigensolver
#'
#' Compute eigenvalues/vectors of the adjacency matrix of an undirected graph.
#'
#' This implements Implicitly Restart Lanczos method described in the paper
#'
#' D. Calvetti and L. Reichel and D. C. Sorensen: An Implicitly Restarted
#' Lanczos method for Large Symmetric Eigenvalue problems, Jounral of ETNA,
#' 1994.
#'
#' @param graph The FlashGraphR object
#' @param which Specify which eigenvalues/vectors to compute, character
#'              constant with exactly two characters.
#' Possible values for symmetric input matrices:
#' "LA' Compute 'nev' largest (algebraic) eigenvalues.
#' "SA" Compute "nev" smallest (algebraic) eigenvalues.
#' "LM" Compute `nev' largest (in magnitude) eigenvalues.
#' "SM" Compute `nev' smallest (in magnitude) eigenvalues.
#' @param nev Numeric scalar. The number of eigenvalues to be computed.
#' @param ncv Number of Lanczos vectors to be generated.
#' @return A named list with the following members:
#'         values: Numeric vector, the desired eigenvalues.
#'         vectors: Numeric matrix, the desired eigenvectors as columns.
#' @name fg.eigen
#' @author Da Zheng <dzheng5@@jhu.edu>
#' @references
#' D. Calvetti and L. Reichel and D. C. Sorensen: An Implicitly Restarted
#' Lanczos method for Large Symmetric Eigenvalue problems, Jounral of ETNA,
#' 1994.
fg.eigen <- function(graph, which="LM", nev=1, ncv=2)
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	stopifnot(!graph$directed)
	.Call("R_FG_eigen_uw", graph, which, as.integer(nev), as.integer(ncv),
		  PACKAGE="FlashGraphR")
}

fg.SVD <- function(graph, which="LM", nev=1, ncv=2, type="LS")
{
	stopifnot(!is.null(graph))
	stopifnot(class(graph) == "fg")
	.Call("R_FG_SVD_uw", graph, which, as.integer(nev), as.integer(ncv),
		  type, PACKAGE="FlashGraphR")
}

#' Spectral embedding
#'
#' This computes spectral embedding of a given graph.
#'
#' It supports five variants.
#'
#' Adjacency matrix (A),
#'
#' A - c * D (AcD), where c is a constant provided by a user.
#'
#' Laplacian matrix (L), L = D - A,
#'
#' Normalized Laplacian matrix (nL), nL = I - D^(-1/2) * A * D^(-1/2)
#'
#' Regularized Laplacian matrix (nL_tau), nL_tau = D_tau^(-1/2) A D_tau^(-1/2)
#'
#' @param fg The FlashGraphR object
#' @param num.eigen The number of eigenvalues/vectors required.
#' @param which The type of the embedding.
#' @param which.eigen Specify which eigenvalues/vectors to compute,
#'                    character constant with exactly two characters.
#' @param c The constant used in the AcD variant of embedding.
#' @param tau The value used in the regularized Laplacian embedding.
#' @param tol Numeric scalar. Stopping criterion: the relative accuracy
#' of the Ritz value is considered acceptable if its error is less than
#' tol times its estimated value. If this is set to zero then machine
#' precision is used.
#' @return A named list with the following members:
#'
#' values
#'
#'    Numeric vector, the desired eigenvalues.
#'
#' vectors
#'
#'    Numeric matrix, the desired eigenvectors as columns. If complex=TRUE
#'   (the default for non-symmetric problems), then the matrix is complex.
#'
#' options
#'
#'   A named list with the supplied options and some information about
#'   the performed calculation, including an ARPACK exit code.
#' @name fg.ase
#' @author Da Zheng <dzheng5@@jhu.edu>
fg.ASE <- function(fg, num.eigen, which=c("A, AcD, L, nL, nL_tau"),
				   which.eigen=c("LM, LA, SM, SA"), c=1, tau=1, tol=1.0e-12)
{
	stopifnot(!is.null(fg))
	stopifnot(class(fg) == "fg")
	# multiply function for eigen on the adjacency matrix
	# this is the default setting.
	multiply <- function(x, extra)
	{
		fg.multiply(fg, x)
	}
	multiply.right <- function(m)
	{
		fg.multiply.matrix(fg, m, TRUE)
	}
	multiply.left.diag <- function(v, m)
	{
		apply(m, 2, function(x) x * v)
	}

	directed = fg.is.directed(fg)
	if (which == "A") {
		if (directed) {
			multiply <- function(x, extra)
			{
				t <- fg.multiply(fg, x, TRUE)
				fg.multiply(fg, t)
			}
		}
		else {
			multiply <- function(x, extra)
			{
				fg.multiply(fg, x)
			}
		}
	}
	else if (which == "AcD") {
		cd <- fg.degree(fg) * c
		if (directed) {
			multiply <- function(x, extra)
			{
				t <- fg.multiply(fg, x, TRUE) + cd * x
				fg.multiply(fg, t) + cd * t
			}

			multiply.right <- function(m)
			{
				fg.multiply.matrix(fg, m, TRUE) + multiply.left.diag(cd, m)
			}
		}
		else {
			multiply <- function(x, extra)
			{
				fg.multiply(fg, x) + cd * x
			}
		}
	}
	else if (which == "L") {
		d <- fg.degree(fg)
		# multiply function for eigen on the laplacian matrix.
		if (directed) {
			multiply <- function(x, extra)
			{
				t <- d * x - fg.multiply(fg, x, TRUE)
				d * t - fg.multiply(fg, t)
			}

			multiply.right <- function(m)
			{
				multiply.left.diag(d, m) - fg.multiply.matrix(fg, m, TRUE)
			}
		}
		else {
			multiply <- function(x, extra)
			{
				d * x - fg.multiply(fg, x)
			}
		}
	}
	else if (which == "nL") {
		d <- 1/sqrt(fg.degree(fg))
		# multiply function for eigen on the normalized laplacian matrix.
		if (directed) {
			multiply <- function(x, extra)
			{
				# D^(-1/2) * L * D^(-1/2) * x
				t <- x - d * fg.multiply(fg, d * x, TRUE)
				t - d * fg.multiply(fg, d * t)
			}

			multiply.right <- function(m)
			{
				m - multiply.left.diag(d,
					fg.multiply.matrix(fg, multiply.left.diag(d, m), TRUE))
			}
		}
		else {
			multiply <- function(x, extra)
			{
				# D^(-1/2) * L * D^(-1/2) * x
				x - d * fg.multiply(fg, d * x)
			}
		}
	}
	else if (which == "nL_tau") {
		d <- fg.degree(fg)
		d_tau <- 1/sqrt(fg.degree(fg) + tau)
		# multiply function for eigen on the normalized laplacian matrix.
		if (directed) {
			multiply <- function(x, extra)
			{
				# D^(-1/2) * L * D^(-1/2) * x
				t <- d_tau * d * d_tau * x - d_tau * fg.multiply(fg, d_tau * x, TRUE)
				d_tau * d * d_tau * t - d_tau * fg.multiply(fg, d_tau * t)
			}

			multiply.right <- function(m)
			{
				multiply.left.diag(d_tau * d * d_tau, m) - multiply.left.diag(
						d_tau, fg.multiply.matrix(fg, multiply.left.diag(d_tau, m), TRUE))
			}
		}
		else {
			multiply <- function(x, extra)
			{
				# D^(-1/2) * L * D^(-1/2) * x
				d_tau * d * d_tau * x - d_tau * fg.multiply(fg, d_tau * x)
			}
		}
	}
	else {
		print("wrong option")
		stopifnot(FALSE)
	}

	norm.col <- function(x)
	{
		x / sqrt(sum(x * x))
	}
	arpack.opts <- list(n=fg.vcount(fg), which=which.eigen,
						nev=num.eigen, ncv=2 * num.eigen, tol=tol)
	ret <- arpack(multiply, sym=TRUE, options=arpack.opts)
	if (directed) {
		list(values=sqrt(ret$values), left=ret$vectors,
			 right=apply(multiply.right(ret$vectors), 2, norm.col), options=ret$options)
	}
	else {
		ret
	}
}

#' Perform spectral clustering
#'
#' Spectral clustering partitions an undirected graph into k groups based on
#' an embedding of the graph, which computes eigenvectors from a variant
#' of the adjacency matrix of the graph. Spectral clustering runs
#' KMeans clustering on eigenvectors.
#'
#' @param fg The FlashGraphR object of the undirected graph.
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
	cat("FlashGraphR ", fg$name, " (", directed, "): ", fg.vcount(fg), " ", fg.ecount(fg),
		"\n", sep="")
}

#' Perform k-means clustering on a data matrix.
#'
#' @param mat A numeric matrix of data.
#' @param k The number of clusters.
#' @param max.iters The maximum number of iterations allowed.
#' @param init The form of initialization to use when the algorithm begins.
#'              The default is "random". For a desciption of each see:
#'              http://en.wikipedia.org/wiki/K-means_clustering#Initialization_methods
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
#
#' @name fg.kmeans
#' @author Disa Mhembere <disa@@jhu.edu>
fg.kmeans <- function(mat, k, max.iters=.Machine$integer.max, init=c("random", "forgy","kmeanspp"))
{
    stopifnot(mat != NULL)
    stopifnot(class(mat) == "matrix")
    .Call("R_FG_kmeans", as.matrix(mat), as.integer(k),
          as.integer(max.iters), init, PACKAGE="FlashGraphR")
}

.onLoad <- function(libname, pkgname)
{
	library(Rcpp)
	library.dynam("FlashGraphR", pkgname, libname, local=FALSE);
	ret <- .Call("R_FG_init", paste(pkgname, ".conf", sep=""), PACKAGE="FlashGraphR")
	stopifnot(ret)
}
