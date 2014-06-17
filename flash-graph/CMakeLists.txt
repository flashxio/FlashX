project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(graph STATIC
	FGlib.cpp
	graph_engine.cpp
	message_processor.cpp
	partitioner.cpp
	vertex.cpp
	worker_thread.cpp
	graph.cpp
	load_balancer.cpp
	messaging.cpp
	vertex_compute.cpp
	vertex_program.cpp
)

subdirs(graph-bfs
	graph-stat
	k-core
	libgraph-algs
	LOL
	matrix
	overlap
	sssp
	sstsg
	test-algs
	tools
)
