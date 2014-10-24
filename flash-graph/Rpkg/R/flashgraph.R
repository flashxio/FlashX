#package.skeleton('flashgraph', code_files = 'flashgraph.R')

fg.set.conf <- function(conf.file)
{
	.Call("R_FG_destroy", PACKAGE="FlashGraph")
	ret <- .Call("R_FG_init", conf.file, PACKAGE="FlashGraph")
}

fg.set.log.level <- function(level)
{
	.Call("R_FG_set_log_level", level, PACKAGE="FlashGraph")
}

fg.list.graphs <- function()
{
	.Call("R_FG_list_graphs", PACKAGE="FlashGraph")
}

fg.exist.graph <- function(graph)
{
	.Call("R_FG_exist_graph", graph, PACKAGE="FlashGraph")
}

fg.get.params <- function(name)
{
	.Call("R_FG_get_params", name, PACKAGE="FlashGraph")
}

# The graph object has three members:
#	name: the graph name.
#	cindex: whether there exists a compressed index for the graph.
#	directed: whether the graph is directed.
# In the future, it may have more information kept in the object.
fg.get.graph <- function(graph)
{
	stopifnot(fg.exist.graph(graph))
	.Call("R_FG_get_graph_obj", graph, PACKAGE="FlashGraph")
}

fg.vcount <- function(graph)
{
	stopifnot(graph != NULL)
	graph$vcount
}

fg.ecount <- function(graph)
{
	stopifnot(graph != NULL)
	graph$ecount
}

fg.is.directed <- function(graph)
{
	stopifnot(graph != NULL)
	graph$directed
}

fg.clusters <- function(graph, mode=c("weak", "strong"))
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	if (mode == "weak")
		.Call("R_FG_compute_wcc", graph, PACKAGE="FlashGraph")
	else if (mode == "strong")
		.Call("R_FG_compute_scc", graph, PACKAGE="FlashGraph")
	else
		stop("a wrong mode")
}

fg.transitivity <- function(graph)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_transitivity", graph, PACKAGE="FlashGraph")
}

fg.degree <- function(graph, type="both")
{
	stopifnot(graph != NULL)
	.Call("R_FG_get_degree", graph, type, PACKAGE="FlashGraph")
}

fg.page.rank <- function(graph, no.iters=1000, damping.factor=0.85)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_pagerank", graph, no.iters, damping.factor,
		  PACKAGE="FlashGraph")
}

fg.directed.triangles <- function(graph, type="cycle")
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_directed_triangles", graph, type, PACKAGE="FlashGraph")
}

fg.undirected.triangles <- function(graph)
{
	stopifnot(graph != NULL)
	.Call("R_FG_compute_undirected_triangles", graph, PACKAGE="FlashGraph")
}

fg.topK.scan <- function(graph, order=1, K=1)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_topK_scan", graph, order, K, PACKAGE="FlashGraph")
}

fg.local.scan <- function(graph, order=1)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_local_scan", graph, order, PACKAGE="FlashGraph")
}

fg.coreness <- function(graph)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	# FIXME set the right parameter.
	k.start <- 1
	k.end <- 0
	.Call("R_FG_compute_kcore", graph, k.start, k.end, PACKAGE="FlashGraph")
}

fg.overlap <- function(graph, vids)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_compute_overlap", graph, vids, PACKAGE="FlashGraph")
}

fg.fetch.subgraph <- function(graph, vertices)
{
	stopifnot(graph != NULL)
	edge.list = .Call("R_FG_fetch_subgraph", graph, vertices,
					  PACKAGE="FlashGraph")
	dframe = data.frame(edge.list$src, edge.list$dst)
	graph.data.frame(dframe, graph$directed)
}

fg.diameter <- function(graph, directed=FALSE)
{
	stopifnot(graph != NULL)
	stopifnot(graph$directed)
	.Call("R_FG_estimate_diameter", graph, directed, PACKAGE="FlashGraph")
}

fg.multiply <- function(graph, vec, transpose=FALSE)
{
	stopifnot(graph != NULL)
#	stopifnot(graph$directed)
	.Call("R_FG_multiply_v", graph, vec, transpose, PACKAGE="FlashGraph")
}

.onLoad <- function(libname, pkgname)
{
	library(Rcpp)
	library.dynam("FlashGraph", pkgname, libname, local=FALSE);
	.Call("R_FG_init", paste(libname, "/", pkgname, "/", pkgname, ".conf",
							 sep=""), PACKAGE="FlashGraph")
}
