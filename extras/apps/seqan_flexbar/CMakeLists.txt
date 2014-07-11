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

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")

# Update the list of file names below if you add source files to your application.
add_executable(test_sflex_demultiplex        test_demultiplex.cpp demultiplex.h)
add_executable(test_sflex_trimming           test_trimming.cpp read_trimming.h)
add_executable(test_sflex_adapter            test_adapter.cpp adapter_trimming.h)
add_executable(test_sflex_general_processing test_general_processing.cpp general_processing.h)

foreach (TEST demultiplex trimming adapter general_processing)
    target_link_libraries (test_sflex_${TEST} ${SEQAN_LIBRARIES})
endforeach ()

add_library (seqan_flexlib
             demultiplex.h
             read_trimming.h
             adapter_trimming.h
             general_processing.h
             seqan_flexbar.h
             seqan_flexlib.cpp)
target_link_libraries (seqan_flexlib ${SEQAN_LIBRARIES})

set (KNIME_FLEX_NODES sflexQC sflexAR sflexDMulti sflexFilter)

foreach (TARGET ${KNIME_FLEX_NODES} seqan_flexbar)
    add_executable (${TARGET} seqan_flexbar.cpp)
    target_link_libraries (${TARGET} seqan_flexlib ${SEQAN_LIBRARIES})
endforeach ()

set_target_properties (sflexQC PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=QUALITY_CONTROL")
set_target_properties (sflexAR PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=ADAPTER_REMOVAL")
set_target_properties (sflexDMulti PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=DEMULTIPLEXING")
set_target_properties (sflexFilter PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FILTERING")

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

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (seqan_flexbar)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqanFlexbar in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} ${KNIME_FLEX_NODES} CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_flexbar")
  set (CPACK_PACKAGE_NAME "seqan_flexbar")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn::Flexbar - Flexible barcode and adapter removal for sequencing platforms in SeqAn")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Jochen Singer <jochen.singer@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (seqan_flexbar "seqan_flexbar")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:seqan_flexbar")

