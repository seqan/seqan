# include (CMakePackageConfigHelpers)
# write_basic_package_version_file( ${CMAKE_CURRENT_LIST_DIR}/seqan-config-version-test.cmake VERSION 2.4.1 COMPATIBILITY ExactVersion )

# This is a basic version file for the Config-mode of find_package().
# It is used by write_basic_package_version_file() as input file for configure_file()
# to create a version-file which can be installed along a config.cmake file.
#
# The created file sets PACKAGE_VERSION_EXACT if the current version string and
# the requested version string are exactly the same and it sets
# PACKAGE_VERSION_COMPATIBLE if the current version is >= requested version,
# but only if the requested major version is the same as the current one.
# The variable CVF_VERSION must be set before calling configure_file().

find_file(SEQAN_VERSION_FILE seqan/version.h HINTS ${CMAKE_CURRENT_LIST_DIR}/../../include)

file(STRINGS "${SEQAN_VERSION_FILE}" SEQAN_VERSION_FILE_H REGEX "#define SEQAN_VERSION_(MAJOR|MINOR|PATCH)")
string(REGEX REPLACE "#define SEQAN_VERSION_(MAJOR|MINOR|PATCH) " "" SEQAN_VERSION_STRING "${SEQAN_VERSION_FILE_H}")
string(REGEX REPLACE ";" "." SEQAN_VERSION_STRING "${SEQAN_VERSION_STRING}")

set(PACKAGE_VERSION "${SEQAN_VERSION_STRING}")

if(PACKAGE_VERSION VERSION_LESS PACKAGE_FIND_VERSION)
  set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()

  if("${PACKAGE_VERSION}" MATCHES "^([0-9]+)\\.")
    set(CVF_VERSION_MAJOR "${CMAKE_MATCH_1}")
  else()
    set(CVF_VERSION_MAJOR "${PACKAGE_VERSION}")
  endif()

  if(PACKAGE_FIND_VERSION_MAJOR STREQUAL CVF_VERSION_MAJOR)
    set(PACKAGE_VERSION_COMPATIBLE TRUE)
  else()
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
  endif()

  if(PACKAGE_FIND_VERSION STREQUAL PACKAGE_VERSION)
      set(PACKAGE_VERSION_EXACT TRUE)
  endif()
endif()
