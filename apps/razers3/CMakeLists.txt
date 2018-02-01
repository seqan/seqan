# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/razers3/CMakeLists.txt
#
# CMakeLists.txt file for razers3.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_razers3 CXX)
message (STATUS "Configuring apps/razers3")

set (SEQAN_APP_VERSION "3.5.8")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Disable RazerS3 on Travis CI because of computational resource limits.
if (SEQAN_TRAVIS_BUILD)
    message (STATUS "  RazerS 3 is not built on Travis CI.")
    return ()
endif ()

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (BZip2)
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
target_link_libraries (razers3_simulate_reads ${SEQAN_LIBRARIES})
target_link_libraries (razers3_quality2prob ${SEQAN_LIBRARIES})
#target_link_libraries (razers3_compPHSens ${SEQAN_LIBRARIES})

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
         DESTINATION ${CMAKE_INSTALL_BINDIR})

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

# Both tests are always executed; without OpenMP the first one is serial, too
seqan_add_app_test (razers3)
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

