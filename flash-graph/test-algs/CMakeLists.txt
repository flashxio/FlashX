project (FlashGraph)

add_executable(test_algs test_algs.cpp)
add_executable(test_cc test_cc.cpp)
add_executable(test_eigen test_eigen.cpp)
add_executable(test_lcc test_lcc.cpp)

target_link_libraries(test_algs graph-algs graph safs pthread numa aio profiler)
target_link_libraries(test_cc graph-algs graph safs pthread numa aio profiler)
target_link_libraries(test_eigen matrix graph safs pthread numa aio profiler)
target_link_libraries(test_lcc graph-algs graph safs pthread numa aio profiler)
