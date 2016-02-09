project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

add_library(graph STATIC
	FGlib.cpp
	graph_engine.cpp
	graph.cpp
	in_mem_storage.cpp
	load_balancer.cpp
	message_processor.cpp
	messaging.cpp
	partitioner.cpp
	ts_graph.cpp
	vertex_compute.cpp
	vertex.cpp
	vertex_index.cpp
	vertex_index_reader.cpp
	worker_thread.cpp
	vertex_program.cpp
	utils.cpp
	vertex_index_constructor.cpp
	graph_config.cpp
)

subdirs(libgraph-algs
	matrix
	test-algs
	tools
)
