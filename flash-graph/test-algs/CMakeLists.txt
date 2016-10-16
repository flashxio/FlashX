project (FlashGraph)

add_executable(test_algs test_algs.cpp)
target_link_libraries(test_algs graph-algs graph safs pthread aio)

if (LIBNUMA_FOUND)
    target_link_libraries(test_algs numa)
endif()

if (hwloc_FOUND)
    target_link_libraries(test_algs hwloc)
endif()
