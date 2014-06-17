project (FlashGraph)

add_executable(k-core k_core.cpp)

target_link_libraries(k-core graph safs pthread numa aio)
