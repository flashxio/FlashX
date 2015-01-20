project (FlashGraph)

add_executable(test_algs test_algs.cpp)
target_link_libraries(test_algs graph-algs graph safs common pthread numa aio)

if (STXXL_FOUND)
    target_link_libraries(test_algs stxxl)
endif()

if (ZLIB_FOUND)
    target_link_libraries(test_algs z)
endif()
