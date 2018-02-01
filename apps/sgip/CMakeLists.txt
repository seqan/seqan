# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/sgip/CMakeLists.txt
#
# CMakeLists.txt file for SGIP.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_sgip CXX)
message (STATUS "Configuring apps/sgip")

set (SEQAN_APP_VERSION "1.4.8")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

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
         DESTINATION ${CMAKE_INSTALL_BINDIR})

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

