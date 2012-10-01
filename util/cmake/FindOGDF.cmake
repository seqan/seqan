# - Try to find OGDF (Open Graph Drawing Framework)
# Once done, this will define
#
#  OGDF_FOUND - system has OGDF
#  OGDF_INCLUDE_DIRS - the OGDF include directories
#  OGDF_LIBRARIES - link these to use OGDF

include(LibFindMacros)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(OGDF_PKGCONF OGDF)

# Include dir
find_path(OGDF_INCLUDE_DIR
  NAMES ogdf/basic/Graph.h
  PATHS ${OGDF_PKGCONF_INCLUDE_DIRS}
        $ENV{OGDF_INCLUDE_PATH}
)

# Finally the library itself
find_library(OGDF_LIBRARY
  NAMES OGDF
  PATHS ${OGDF_PKGCONF_LIBRARY_DIRS}
        $ENV{OGDF_LIBRARY_PATH}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(OGDF_PROCESS_INCLUDES OGDF_INCLUDE_DIR OGDF_INCLUDE_DIRS)
set(OGDF_PROCESS_LIBS OGDF_LIBRARY OGDF_LIBRARIES)
message(STATUS ${OGDF_LIBRARIES})
libfind_process(OGDF)

