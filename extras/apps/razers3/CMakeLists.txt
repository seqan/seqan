# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/razers3/CMakeLists.txt
#
# CMakeLists.txt file for razers3.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (extras_apps_razers3)
message (STATUS "Configuring extras/apps/razers3")

# Disable RazerS3 for MinGW and VS8 and earlier.
if (MINGW OR (MSVC AND MSVC_VERSION LESS 1500))
    message (STATUS "  RazerS 3 cannot be built with MinGW or with MS Visual C++ 8 and earlier.")
    return ()
endif (MINGW OR (MSVC AND MSVC_VERSION LESS 1500))

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (razers3 razers.cpp
                        razers.h
                        job_queue.h
                        outputFormat.h
                        parallel_job_queue.h
                        parallel_misc.h
                        parallel_store.h
                        paramChooser.h
                        param_tabs.cpp
                        param_tabs.inc
                        param_tabs.h
                        profile_timeline.h
                        razers_match_filter.h
                        razers_paired_match_filter.h
                        razers_matepairs.h
                        razers_matepairs_parallel.h
                        razers_parallel.h
                        razers_window.h
                        readSimulator.h)

#add_executable (razers3_compPHSens compPHSens.cpp razers.h)

add_executable (razers3_simulate_reads simulate_reads.cpp)

add_executable (razers3_quality2prob quality2prob.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (razers3 ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (razers3)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install razers3 in ${PREFIX}/bin directory
install (TARGETS razers3
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/razers3 for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/genome.fa
               example/reads.fa
               example/reads2.fa
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# The parallel code paths are only available when compiling with OpenMP.
if (OPENMP_FOUND)
  seqan_add_app_test (razers3)
endif (OPENMP_FOUND)

# However, the sequential tests are always executed.
seqan_add_app_test (razers3 _sequential)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable razers in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} razers3 CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers3")
  set (CPACK_PACKAGE_NAME "razers3")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "razers3 - Faster, fully sensitive read mapping.")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (razers3 "razers3")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:razers3")

