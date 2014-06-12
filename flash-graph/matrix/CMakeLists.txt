project (FlashGraph)

include_directories("/usr/include/eigen3")

add_library(matrix STATIC
	matrix_eigensolver.cpp
)
