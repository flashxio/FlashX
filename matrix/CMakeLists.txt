project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(FMatrix STATIC
	dense_matrix.cpp
	EM_vector.cpp
	matrix_worker_thread.cpp
	sparse_matrix.cpp
	EM_dense_matrix.cpp
	matrix_io.cpp
	mem_dense_matrix.cpp
	mem_vector.cpp
	data_frame.cpp
)
