project (FlashGraph)

add_executable(overlap overlap.cpp)

target_link_libraries(overlap graph safs pthread numa aio profiler)
