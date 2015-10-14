project (FlashGraph)

add_executable(el2al el2al.cpp)
add_executable(rmat-gen rmat-gen.cpp)
# ext_mem_vertex_iterator.cpp

target_link_libraries(el2al graph safs common pthread numa aio)

if (ZLIB_FOUND)
    target_link_libraries(el2al z)
endif()

if (hwloc_FOUND)
    target_link_libraries(el2al hwloc)
endif()
