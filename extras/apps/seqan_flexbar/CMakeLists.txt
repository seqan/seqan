# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /extras/apps/seqan_flexbar/CMakeLists.txt
#
# CMakeLists.txt file for seqan_flexbar.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_extras_apps_seqan_flexbar)
message (STATUS "Configuring extras/apps/seqan_flexbar")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES OpenMP ZLIB BZip2)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
#add_executable (seqan_flexbar seqan_flexbar.cpp
#                seqan_flexbar.h
#                demultiplex.h
#			    readTrimming.h
#			    adapterTrimming.h
#			    generalProcessing.h)
add_executable(demultiplex_test demultiplex_test.cpp demultiplex.h)
add_executable(trimming_test trimming_test.cpp readTrimming.h)
add_executable(adapter_test adapter_test.cpp adapterTrimming.h)
add_executable(generalProcessing_test generalProcessing_test.cpp generalProcessing.h)

add_library (seqan_flexbar_lib
             demultiplex.h
             readTrimming.h
             adapterTrimming.h
             generalProcessing.h
             seqan_flexbar.h
             seqan_flexbar.cpp)
target_link_libraries (seqan_flexbar_lib ${SEQAN_LIBRARIES})

add_executable (seqFlexQC seqan_flexbar_quality_control.cpp)
add_executable (seqFlexAR seqan_flexbar_adapter_removal.cpp)
add_executable (seqFlexDMulti seqan_flexbar_demultiplexing.cpp)
add_executable (seqFlexFilter seqan_flexbar_filtering.cpp)
add_executable (seqan_flexbar seqan_flexbar_all.cpp)

# Add dependencies found by find_package (SeqAn).
foreach (EXEC trimming_test
         demultiplex_test adapter_test generalProcessing_test
         seqan_flexbar
         seqFlexQC
         seqFlexAR
         seqFlexDMulti
         seqFlexFilter)
    target_link_libraries (${EXEC} seqan_flexbar_lib ${SEQAN_LIBRARIES})
endforeach ()

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (seqan_flexbar)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install seqan_flexbar in ${PREFIX}/bin directory
install (TARGETS seqan_flexbar
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/seqan_flexbar for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (seqan_flexbar)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqanFlexbar in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} seqan_flexbar CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_flexbar")
  set (CPACK_PACKAGE_NAME "seqan_flexbar")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "seqan_flexbar")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (seqan_flexbar "seqan_flexbar")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_flexbar")

