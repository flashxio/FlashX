project (FlashGraph)

link_directories (${PROJECT_BINARY_DIR}/libsafs) 

add_executable(SAFS-util SAFS-util.cpp)

target_link_libraries (SAFS-util safs pthread numa aio)

if (hwloc_FOUND)
	target_link_libraries(SAFS-util hwloc)
endif()
