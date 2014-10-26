project (FlashGraph)

add_executable(el2al el2al.cpp ext_mem_vertex_iterator.cpp)
add_executable(rmat-gen rmat-gen.cpp)
# ext_mem_vertex_iterator.cpp

target_link_libraries(el2al graph safs common pthread numa aio stxxl z)
