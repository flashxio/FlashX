project (FlashGraph)

include_directories("${PROJECT_SOURCE_DIR}")

add_executable(el2al el2al.cpp)
target_link_libraries(el2al FMatrix graph safs common pthread numa aio cblas)

add_executable(al22d al22d.cpp fg_utils.cpp)
target_link_libraries(al22d FMatrix graph safs common pthread numa aio cblas)

find_package(hwloc)
if (hwloc_FOUND)
	target_link_libraries(el2al hwloc)
	target_link_libraries(al22d hwloc)
endif()
