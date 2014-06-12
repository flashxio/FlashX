project (FlashGraph)

add_executable(graph-bfs bfs_graph.cpp)
add_executable(undirected-graph-bfs bfs_undirected_graph.cpp)

target_link_libraries(graph-bfs graph safs pthread numa aio profiler)
target_link_libraries(undirected-graph-bfs graph safs pthread numa aio profiler)
