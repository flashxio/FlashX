#package.skeleton('flashgraph', code_files = 'flashgraph.R')

get_fg_graph <- function(graph.file, index.file, conf.file)
{
	list(graph.file, index.file, conf.file)
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

.onLoad <- function(libname, pkgname) {
	library.dynam("FlashGraph", pkgname, libname, local=FALSE);
	.Call("R_FG_init", PACKAGE="FlashGraph")
}
