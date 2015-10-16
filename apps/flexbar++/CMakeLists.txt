# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/flexbar++/CMakeLists.txt
#
# CMakeLists.txt file for Flexbar++.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_apps_flexbar++)
message (STATUS "Configuring apps/flexbar++")

set (SEQAN_APP_VERSION "1.0.0")
add_definitions ("-D_MULTITHREADED_IO")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package(CXX11 REQUIRED)

# Disable flexbar if no CXX11 available.
if (NOT CXX11_FOUND)
    message (STATUS "  C++11 compiler required to build flexbar++")
    return ()
endif (NOT CXX11_FOUND)


find_package(CXX14 REQUIRED)

# Disable flexbar if no CXX11 available.
if (NOT CXX14_FOUND)
    message (STATUS "  C++14 compiler required to build flexbar++")
    return ()
endif (NOT CXX14_FOUND)

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package (SeqAn REQUIRED)

find_package(Boost 1.36.0)
if(NOT Boost_FOUND)
    message (STATUS "  flexbar cannot be built without Boost.")
    return ()
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")
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

add_library (flexlib
             read.h
			 read_reader.h
			 read_writer.h
			 read_processor.h
			 semaphore.h
             demultiplex.h
             read_trimming.h
             adapter_trimming.h
             general_processing.h
             helper_functions.h
			 general_stats.h
             flexbar++.h
             flexlib.cpp)
target_link_libraries (flexlib ${SEQAN_LIBRARIES})

set (KNIME_FLEX_NODES sflexQC sflexAR sflexDMulti sflexFilter)

foreach (TARGET ${KNIME_FLEX_NODES} flexbar++)
    add_executable (${TARGET} flexbar++.cpp)
    target_link_libraries (${TARGET} flexlib ${SEQAN_LIBRARIES})
endforeach ()

set_target_properties (sflexQC PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=QUALITY_CONTROL")
set_target_properties (sflexAR PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=ADAPTER_REMOVAL")
set_target_properties (sflexDMulti PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=DEMULTIPLEXING")
set_target_properties (sflexFilter PROPERTIES COMPILE_DEFINITIONS "FLEX_PROG=FILTERING")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${CXX14_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (flexbar++)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install flexbar++ in ${PREFIX}/bin directory
install (TARGETS flexbar++
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/flexbar++ for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (flexbar++)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable seqanFlexbar in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} ${KNIME_FLEX_NODES} CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:flexbar++")
  set (CPACK_PACKAGE_NAME "flexbar++")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn::Flexbar++ - Flexible barcode and adapter removal for sequencing platforms in SeqAn")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Benjamin Menkuec <benjamin@menkuec.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (flexbar++ "flexbar++")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:flexbar++")

