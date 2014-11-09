cmake_minimum_required (VERSION 2.6)

project (FlashGraph)

# The version number.
set (FlashGraph_VERSION_MAJOR 0)
set (FlashGraph_VERSION_MINOR 1)

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

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++0x -O3 -fPIC" )
#set(CMAKE_BUILD_TYPE Release)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}" include libcommon)

subdirs(libcommon libsafs utils flash-graph)
