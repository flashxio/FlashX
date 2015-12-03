project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

add_library(eigen STATIC
	block_dense_matrix.cpp
	eigensolver.cpp
	collected_col_matrix_store.cpp
)
