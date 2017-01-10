project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

add_executable(el2fg el2fg.cpp)
target_link_libraries(el2fg graph FMatrix safs pthread cblas)

add_executable(fg2fm fg2fm.cpp)
target_link_libraries(fg2fm graph FMatrix safs pthread cblas)

if (LIBNUMA_FOUND)
    target_link_libraries(el2fg numa)
    target_link_libraries(fg2fm numa)
endif()

if (LIBAIO_FOUND)
    target_link_libraries(el2fg aio)
    target_link_libraries(fg2fm aio)
endif()

find_package(hwloc)
if (hwloc_FOUND)
	target_link_libraries(el2fg hwloc)
	target_link_libraries(fg2fm hwloc)
endif()

if (ZLIB_FOUND)
	target_link_libraries(el2fg z)
	target_link_libraries(fg2fm z)
endif()
