# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/sgip/CMakeLists.txt
#
# CMakeLists.txt file for SGIP.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_sgip)
message (STATUS "Configuring extras/apps/sgip")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# The SGIP program can only be built with MSVC++ 10 or later.
find_package(CXX11 REQUIRED)

# Disable SGIP for MinGW and VS8 and earlier.
if (NOT CXX11_FOUND)
    message (STATUS "  SGIP cannot be built with MinGW or with MS Visual C++ 8 and earlier.")
    return ()
endif (NOT CXX11_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (sgip sgip.cpp sgip.h sgip_base.h sgip_output.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (sgip ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (sgip)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install sgip in ${PREFIX}/bin directory
install (TARGETS sgip
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/sgip for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/Iso_Data/iso_m2D_m196.A01
               example/r01/iso_r01_m200.A00
               example/r01/iso_r01_m200.A01
               example/r01/iso_r01_m200.B00
               example/r01/iso_r01_m200.B01
               example/srg/latin-4
               example/srg/lattice-4
               example/srg/sts-7
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (sgip)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:sgip")
  set (CPACK_PACKAGE_NAME "sgip")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SGIP - Solution to the Graph Isomorphism Problem")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (sgip "sgip")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:sgip")

