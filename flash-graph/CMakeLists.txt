project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(stxxl)
if (STXXL_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_STXXL")
endif()

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
	utils_io.cpp
)

subdirs(libgraph-algs
	matrix
	test-algs
	tools
)
