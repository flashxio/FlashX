project (FlashGraph)

link_directories (${PROJECT_BINARY_DIR}/libsafs) 

add_executable(SAFS-util SAFS-util.cpp)

target_link_libraries (SAFS-util safs pthread)

if (LIBNUMA_FOUND)
    target_link_libraries(SAFS-util numa)
endif()

if (LIBAIO_FOUND)
    target_link_libraries(SAFS-util aio)
endif()

if (hwloc_FOUND)
	target_link_libraries(SAFS-util hwloc)
endif()
