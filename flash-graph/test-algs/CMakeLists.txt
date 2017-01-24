project (FlashGraph)

add_executable(test_algs test_algs.cpp)
target_link_libraries(test_algs graph-algs graph FMatrix safs cblas pthread)

find_package(ZLIB)
if (ZLIB_FOUND)
	target_link_libraries(test_algs z)
endif()

if (LIBNUMA_FOUND)
    target_link_libraries(test_algs numa)
endif()

if (LIBAIO_FOUND)
    target_link_libraries(test_algs aio)
endif()

if (hwloc_FOUND)
    target_link_libraries(test_algs hwloc)
endif()
