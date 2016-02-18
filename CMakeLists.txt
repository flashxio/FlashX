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

find_package(Boost COMPONENTS system filesystem	REQUIRED)
find_package(Threads REQUIRED)
find_package(LibNuma REQUIRED)
find_package(LibAio REQUIRED)
find_package(LibRt REQUIRED)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++0x -O3 -g -fPIC" )

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(hwloc)
if (hwloc_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_HWLOC")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_HWLOC")
endif()

check_c_compiler_flag("-mavx" HAVE_FLAG_M_AVX)
if(HAVE_FLAG_M_AVX)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
endif()

#set(CMAKE_BUILD_TYPE Release)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}" include libcommon flash-graph)

subdirs(libcommon libsafs utils flash-graph matrix)
