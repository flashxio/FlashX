# Copyright (c) 2012 Marwan Abdellah <marwan.abdellah@epfl.ch>
#                    Daniel Nachbaur <daniel.nachbaur@epfl.ch>
#               2013 Stefan.Eilemann@epfl.ch

# Use pkg-config to fetch the contents of the .pc file
# After that, use the directories refer to the libraries and
# also the headers

if(NOT PKGCONFIG_FOUND)
  find_package(PkgConfig QUIET)
endif()

if(HWLOC_ROOT)
  set(ENV{PKG_CONFIG_PATH} "${HWLOC_ROOT}/lib/pkgconfig")
else()
  foreach(PREFIX ${CMAKE_PREFIX_PATH})
    set(PKG_CONFIG_PATH "${PKG_CONFIG_PATH}:${PREFIX}/lib/pkgconfig")
  endforeach()
  set(ENV{PKG_CONFIG_PATH} "${PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
endif()

if(hwloc_FIND_REQUIRED)
  set(_hwloc_OPTS "REQUIRED")
endif()
if(hwloc_FIND_QUIETLY)
  set(_hwloc_OPTS "QUIET")
endif()
if(hwloc_FIND_REQUIRED AND hwloc_FIND_QUIETLY)
  set(_hwloc_OPTS "REQUIRED QUIET")
endif()
if(NOT hwloc_FIND_QUIETLY)
  set(_hwloc_output 1)
endif()

if(hwloc_FIND_VERSION)
  if(hwloc_FIND_VERSION_EXACT)
    pkg_check_modules(hwloc ${_hwloc_OPTS} hwloc=${hwloc_FIND_VERSION})
  else()
    pkg_check_modules(hwloc ${_hwloc_OPTS} hwloc>=${hwloc_FIND_VERSION})
  endif()
else()
  pkg_check_modules(hwloc ${_hwloc_OPTS} hwloc)
endif()

if(hwloc_FOUND)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(hwloc DEFAULT_MSG hwloc_LIBRARIES
    hwloc_INCLUDE_DIR)

  if(NOT ${hwloc_VERSION} VERSION_LESS 1.7.0)
    set(hwloc_GL_FOUND 1)
  endif()

  if(_hwloc_output)
    message(STATUS
      "Found hwloc ${hwloc_VERSION} in ${hwloc_INCLUDE_DIR}:${hwloc_LIBRARIES}")
  endif()
endif()
