cmake_minimum_required (VERSION 2.6)

project (FlashGraph)

include(CheckCCompilerFlag)

# The version number.
set (FlashGraph_VERSION_MAJOR 0)
set (FlashGraph_VERSION_MINOR 3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6.3")
		message(FATAL_ERROR "Insufficient gcc version")
	endif()
else()
	message(FATAL_ERROR "We require GCC")
endif()

find_package(Boost REQUIRED)
find_package(Threads REQUIRED)
find_package(LibRt REQUIRED)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++0x -O3 -fPIC" )

execute_process(COMMAND bash ${PROJECT_SOURCE_DIR}/get_num_nodes.sh OUTPUT_VARIABLE num_nodes)
message("#node: ${num_nodes}")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (${num_nodes} GREATER 1)
	find_package(hwloc)
	if (hwloc_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HWLOC")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HWLOC")
	endif()
endif()

if (${num_nodes} GREATER 1)
	find_package(LibNuma)
	if (LIBNUMA_FOUND)
		message(STATUS "Find libnuma.")
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_NUMA")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_NUMA")
	endif()
endif()

find_package(LibAio)
if (LIBAIO_FOUND)
	message(STATUS "Find libaio.")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_LIBAIO")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_LIBAIO")
endif()

check_c_compiler_flag("-mavx" HAVE_FLAG_M_AVX)
if(HAVE_FLAG_M_AVX)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
endif()

#set(CMAKE_BUILD_TYPE Release)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}"  libsafs flash-graph matrix)

subdirs(libsafs utils flash-graph matrix)
