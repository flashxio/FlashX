project (FlashGraph)

add_library(graph-algs STATIC
	diameter_graph.cpp
	local_scan_graph.cpp
	scan_graph.cpp
	topK_scan_graph.cpp
	wcc.cpp
	directed_triangle_graph.cpp
	page_rank.cpp
	scc.cpp
	undirected_triangle_graph.cpp
)
