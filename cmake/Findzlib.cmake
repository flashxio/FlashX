find_path(LIBZ_INCLUDE_DIRS NAMES zlib.h)
find_library(LIBZ_LIBRARIES NAMES z)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibZ DEFAULT_MSG
	LIBZ_LIBRARIES
	LIBZ_INCLUDE_DIRS)

mark_as_advanced(LIBZ_INCLUDE_DIRS LIBZ_LIBRARIES)
