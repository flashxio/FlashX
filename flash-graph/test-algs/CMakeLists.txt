project (FlashGraph)

add_executable(test_algs test_algs.cpp)
target_link_libraries(test_algs graph-algs graph safs pthread numa aio)

if (hwloc_FOUND)
    target_link_libraries(test_algs hwloc)
endif()
