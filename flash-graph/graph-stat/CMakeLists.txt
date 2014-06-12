project (FlashGraph)

add_executable(graph-stat stat_graph.cpp)

target_link_libraries(graph-stat graph safs pthread numa aio profiler)
