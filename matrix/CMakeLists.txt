project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(BLAS REQUIRED)
message( STATUS BLAS found: ${BLAS_LIBRARIES} )

add_library(FMatrix STATIC
	dense_matrix.cpp
	EM_vector.cpp
	sparse_matrix.cpp
	EM_dense_matrix.cpp
	matrix_io.cpp
	data_frame.cpp
	generic_type.cpp
	matrix_header.cpp
	sparse_matrix_format.cpp
	hilbert_curve.cpp
	NUMA_dense_matrix.cpp
	NUMA_vector.cpp
	raw_data_array.cpp
	mem_worker_thread.cpp
	local_matrix_store.cpp
	data_io.cpp
	local_vec_store.cpp
	mem_matrix_store.cpp
	mapply_matrix_store.cpp
	mem_vec_store.cpp
	one_val_matrix_store.cpp
	rand_gen.cpp
	fm_utils.cpp
	matrix_store.cpp
	mem_vv_store.cpp
	local_vv_store.cpp
	vec_store.cpp
	vector.cpp
	vector_vector.cpp
	vv_store.cpp
	EM_object.cpp
	local_mem_buffer.cpp
	bulk_operate.cpp
	factor.cpp
	col_vec.cpp
	combined_matrix_store.cpp
	cached_matrix_store.cpp
)

if(ENABLE_TRILINOS)
	subdirs(eigensolver)
endif()
subdirs(utils)
