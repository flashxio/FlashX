find_path(LIBAIO_INCLUDE_DIRS NAMES libaio.h)
find_library(LIBAIO_LIBRARIES NAMES aio)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibAio DEFAULT_MSG
  LIBAIO_LIBRARIES
  LIBAIO_INCLUDE_DIRS)

mark_as_advanced(LIBAIO_INCLUDE_DIRS LIBAIO_LIBRARIES)
