project (FlashGraph)

add_executable(scan-statistics scan_graph.cpp)

target_link_libraries(scan-statistics graph safs pthread numa aio profiler)
