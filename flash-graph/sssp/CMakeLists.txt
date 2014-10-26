project (FlashGraph)

add_executable(sssp sssp.cpp)

target_link_libraries(sssp graph safs common pthread numa aio)
