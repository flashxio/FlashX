project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(ZLIB)
if (ZLIB_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_GZIP")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GZIP")
endif()

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
)

subdirs(libgraph-algs
	matrix
	test-algs
	tools
)
