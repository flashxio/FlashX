#package.skeleton('flashgraph', code_files = 'flashgraph.R')

fg.set.conf <- function(conf.file)
{
	.Call("R_FG_destroy", PACKAGE="FlashGraph")
	.Call("R_FG_init", conf.file, PACKAGE="FlashGraph")
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

fg.clusters <- function(graph, mode=c("weak", "strong"))
{
	if (mode == "weak")
		.Call("R_FG_compute_wcc", graph, PACKAGE="FlashGraph")
	else if (mode == "strong")
		.Call("R_FG_compute_scc", graph, PACKAGE="FlashGraph")
	else
		stop("a wrong mode")
}

fg.transitivity <- function(graph)
{
	.Call("R_FG_compute_transitivity", graph, PACKAGE="FlashGraph")
}

fg.degree <- function(graph, type="both")
{
	.Call("R_FG_get_degree", graph, type, PACKAGE="FlashGraph")
}

fg.page.rank <- function(graph, no.iters=1000, damping.factor=0.85)
{
	.Call("R_FG_compute_pagerank", graph, no.iters, damping.factor,
		  PACKAGE="FlashGraph")
}

fg.directed.triangles <- function(graph, type="cycle")
{
	.Call("R_FG_compute_directed_triangles", graph, type, PACKAGE="FlashGraph")
}

fg.undirected.triangles <- function(graph)
{
	.Call("R_FG_compute_undirected_triangles", graph, PACKAGE="FlashGraph")
}

fg.topK.scan <- function(graph, order=1, K=1)
{
	.Call("R_FG_compute_topK_scan", graph, order, K, PACKAGE="FlashGraph")
}

fg.local.scan <- function(graph, order=1)
{
	.Call("R_FG_compute_local_scan", graph, order, PACKAGE="FlashGraph")
}

fg.coreness <- function(graph)
{
	# FIXME set the right parameter.
	k.start <- 1
	k.end <- 0
	.Call("R_FG_compute_kcore", graph, k.start, k.end, PACKAGE="FlashGraph")
}

fg.overlap <- function(graph, vids)
{
	.Call("R_FG_compute_overlap", graph, vids, PACKAGE="FlashGraph")
}

fg.fetch.subgraph <- function(graph, vertices)
{
	edge.list = .Call("R_FG_fetch_subgraph", graph, vertices,
					  PACKAGE="FlashGraph")
	dframe = data.frame(edge.list$src, edge.list$dst)
	graph.data.frame(dframe, graph$directed)
}

.onLoad <- function(libname, pkgname) {
	library(Rcpp)
	library.dynam("FlashGraph", pkgname, libname, local=FALSE);
	.Call("R_FG_init", paste(libname, "/", pkgname, "/", pkgname, ".conf",
							 sep=""), PACKAGE="FlashGraph")
}
