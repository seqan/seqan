# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/seqcons2/CMakeLists.txt
#
# CMakeLists.txt file for seqcons2.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_seqcons2)
message (STATUS "Configuring extras/apps/seqcons2")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Use C++11.
find_package (CXX11)
if (NOT CXX11_FOUND)
  message (STATUS "  C++11 not found, skipping.")
  return ()
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (seqcons2
                seqcons.cpp
                seqcons_app.h
                seqcons_app.cpp
                seqcons_options.h
                seqcons_options.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (seqcons2 ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (seqcons2)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install seqcons2 in ${PREFIX}/bin directory
install (TARGETS seqcons2
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/seqcons2 for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (seqcons2)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqcons2")
  set (CPACK_PACKAGE_NAME "seqcons2")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "seqcons2")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")

  seqan_configure_cpack_app (seqcons2 "seqcons2")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqcons2")

