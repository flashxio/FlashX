get_fg_graph <- function(graph.file, index.file, conf.file)
{
	list(graph.file, index.file, conf.file)
}

fg.wcc <- function(graph)
{
	.Call("R_FG_compute_wcc", graph)
}

fg.scc <- function(graph)
{
	.Call("R_FG_compute_scc", graph)
}

fg.transitivity <- function(graph)
{
	.Call("R_FG_compute_transitivity", graph)
}

fg.degree <- function(graph, type="both")
{
	.Call("R_FG_get_degree", graph, type)
}

fg.page.rank <- function(graph, no.iters, damping.factor=0.85)
{
	.Call("R_FG_compute_pagerank", graph, no.iters, damping.factor)
}

fg.directed.triangles <- function(graph, type="cycle")
{
	.Call("R_FG_compute_directed_triangles", graph, type)
}

fg.local.scan <- function(graph, k=1)
{
	.Call("R_FG_compute_local_scan", graph, k)
}
