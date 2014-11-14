find_path(LIBSTXXL_INCLUDE_DIRS NAMES stxxl.h)
find_library(LIBSTXXL_LIBRARIES NAMES stxxl)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libstxxl DEFAULT_MSG
	LIBSTXXL_LIBRARIES
	LIBSTXXL_INCLUDE_DIRS)

if (LIBSTXXL_INCLUDE_DIRS)
	set(STXXL_FOUND TRUE)
endif()

mark_as_advanced(LIBSTXXL_INCLUDE_DIRS LIBSTXXL_LIBRARIES)
