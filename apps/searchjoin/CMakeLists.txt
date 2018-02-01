# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/searchjoin/CMakeLists.txt
#
# CMakeLists.txt file for searchjoin.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_searchjoin CXX)
message (STATUS "Configuring apps/searchjoin")

set (SEQAN_APP_VERSION "0.5.8")

if (SEQAN_TRAVIS_BUILD)
    message (STATUS "  Skipping build and test of searchjoin on Travis CI.")
    return ()
endif (SEQAN_TRAVIS_BUILD)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS} -DSEARCHJOIN_HUGEDB=TRUE)

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (s4_search search.cpp finder.h db.h verifier.h writer.h)
add_executable (s4_join join.cpp finder.h db.h verifier.h writer.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (s4_search ${SEQAN_LIBRARIES})
target_link_libraries (s4_join ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (searchjoin)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install searchjoin in ${PREFIX}/bin directory
install (TARGETS s4_search s4_join
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/searchjoin for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (searchjoin)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:searchjoin")
  set (CPACK_PACKAGE_NAME "searchjoin")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "searchjoin")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")

  seqan_configure_cpack_app (searchjoin "searchjoin")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:searchjoin")

