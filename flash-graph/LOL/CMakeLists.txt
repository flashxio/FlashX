project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}/matrix")
include_directories("/usr/include/eigen3")

add_executable(LOL LOL.cpp mnist_io.cpp)
add_executable(m2g m2g.cpp mnist_io.cpp)
add_executable(print_mnist print_mnist.cpp mnist_io.cpp)

target_link_libraries(LOL graph safs pthread numa aio profiler matrix)
target_link_libraries(m2g graph safs pthread numa aio profiler matrix)
target_link_libraries(print_mnist graph safs pthread numa aio profiler matrix)
