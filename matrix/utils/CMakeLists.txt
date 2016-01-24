project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

add_executable(el2fg el2fg.cpp)
target_link_libraries(el2fg FMatrix graph safs common pthread numa aio cblas)

add_executable(fg2fm fg2fm.cpp fg_utils.cpp)
target_link_libraries(fg2fm FMatrix graph safs common pthread numa aio cblas)

find_package(hwloc)
if (hwloc_FOUND)
	target_link_libraries(el2fg hwloc)
	target_link_libraries(fg2fm hwloc)
endif()

if (ZLIB_FOUND)
	target_link_libraries(el2fg z)
	target_link_libraries(fg2fm z)
endif()
