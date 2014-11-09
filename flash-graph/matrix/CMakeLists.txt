project (FlashGraph)

find_package(Eigen3)
if (EIGEN3_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_EIGEN")

	add_library(matrix STATIC
		matrix_eigensolver.cpp
	)
endif()
