# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: apps/fiona/CMakeLists.txt
#
# CMakeLists.txt file for fiona.
# ===========================================================================

cmake_minimum_required (VERSION 3.12)
project (seqan_apps_fiona CXX)
message (STATUS "Configuring apps/fiona")

set (SEQAN_APP_VERSION "0.2.13")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP COMPONENTS CXX)
    find_package (ZLIB)
    find_package (Boost)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# Stop here if we cannot find Boost or OpenMP.
if (NOT Boost_FOUND)
    message (STATUS "  Skipping fiona, missing Boost")
    return ()
endif ()

if (NOT OPENMP_FOUND)
    message (STATUS "  Skipping fiona, missing OpenMP")
    return ()
endif ()

if (NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message (STATUS "  Skipping fiona on 32-bit architectures")
    return ()
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})
#add_definitions("-DQGRAM_LENGTH=5")

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable (compute_gain compute_gain.cpp)
add_executable (fiona fiona.cpp)
add_executable (fiona_illumina fiona.cpp)

set_target_properties (fiona_illumina PROPERTIES COMPILE_DEFINITIONS "FIONA_ILLUMINA")

# Boost 1.74.0 is not compatible with cpp20
set_property(TARGET fiona PROPERTY CXX_STANDARD 17)
set_property(TARGET fiona_illumina PROPERTY CXX_STANDARD 17)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (compute_gain ${SEQAN_LIBRARIES})
target_link_libraries (fiona ${SEQAN_LIBRARIES})
target_link_libraries (fiona_illumina ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (fiona)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install fiona in ${PREFIX}/bin directory
install (TARGETS fiona fiona_illumina compute_gain DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/fiona for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/reads.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------
if (NOT MSVC)
    seqan_add_app_test (fiona)
else ()
    message (STATUS "  Skipping fiona tests on Visual Studio")
endif ()

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable fiona in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} fiona fiona_illumina compute_gain CACHE INTERNAL "")

#Add test to check correct ctd generation.
seqan_add_ctd_test (fiona)
seqan_add_ctd_test (fiona_illumina)
seqan_add_ctd_test (compute_gain)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:fiona")
  set (CPACK_PACKAGE_NAME "fiona")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "fiona")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "David Weese <david.weese@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (fiona "fiona")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:fiona")
