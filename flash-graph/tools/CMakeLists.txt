project (FlashGraph)

add_executable(el2al el2al.cpp ext_mem_vertex_iterator.cpp)
add_executable(graph_merger graph_merger.cpp ext_mem_vertex_iterator.cpp)
add_executable(print_ts_graph print_ts_graph.cpp)
add_executable(rmat-gen rmat-gen.cpp)
# ext_mem_vertex_iterator.cpp

target_link_libraries(el2al graph safs pthread numa aio)
target_link_libraries(graph_merger graph safs pthread numa aio)
target_link_libraries(print_ts_graph graph safs pthread numa aio)
