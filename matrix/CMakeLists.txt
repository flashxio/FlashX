project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

find_package(BLAS REQUIRED)
message( STATUS BLAS found: ${BLAS_LIBRARIES} )

find_package(ZLIB)
if (ZLIB_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_GZIP")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GZIP")
endif()


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
	matrix_config.cpp
	agg_matrix_store.cpp
	matrix_stats.cpp
	block_matrix.cpp
	IPW_matrix_store.cpp
	sink_matrix.cpp
	project_matrix_store.cpp
	EM_vv_store.cpp
	groupby_matrix_store.cpp
	materialize.cpp
	sub_matrix_store.cpp
	set_data_matrix_store.cpp
	fm_utils.cpp
	set_rc_matrix_store.cpp
	cum_matrix.cpp
)

if(ENABLE_TRILINOS)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DENABLE_TRILINOS")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_TRILINOS")
	subdirs(eigensolver)
endif()
