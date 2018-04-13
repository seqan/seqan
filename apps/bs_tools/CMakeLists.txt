# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/bs_tools/CMakeLists.txt
#
# CMakeLists.txt file for bs_tools.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_bs_tools CXX)
message (STATUS "Configuring apps/bs_tools")

set (SEQAN_APP_VERSION "0.1.10")

if ((${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") AND ("$ENV{MODEL}" STREQUAL "Nightly"))
    message (STATUS "bs_tools skipped on FreeBSD because math.h rounding errors.")
    return ()
endif ()

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (ZLIB)
    find_package (Boost)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

if (NOT Boost_FOUND OR NOT ZLIB_FOUND)
    message (STATUS "  Boost or zlib not found: Not building bs_tools.")
    return ()
endif (NOT Boost_FOUND OR NOT ZLIB_FOUND)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})
include_directories (SYSTEM ${Boost_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Update the list of file names below if you add source files to your application.
add_executable (bisar bisar.cpp
                    bisar.h
                    bisar_base.h
                    bisar_score.h
                    bisar_score_data.h)

add_executable (casbar casbar.cpp
                    casbar.h
                    casbar_calling.h
                    casbar_realignment.h
                    casbar_consensus_realign.h
                    casbar_score.h
                    casbar_score_data.h
                    casbar_util.h
                    casbar_alphabets.h)

add_executable (four2three four2three.cpp
                    four2three.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (bisar ${SEQAN_LIBRARIES})
target_link_libraries (casbar ${Boost_LIBRARIES} ${SEQAN_LIBRARIES})
target_link_libraries (four2three ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include bs_tools executables in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES}
        bisar
        casbar
        four2three
        CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (bs_tools)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install bs_tools in ${PREFIX}/bin directory
install (TARGETS bisar casbar four2three
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/bs_tools for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "DEVELOP")
    include("SeqAnConsistentFPFlags")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CONSISTENT_FP_FLAGS}")
endif()

# (weese:) temporarily disabled app test
seqan_add_app_test (bs_tools)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:bs_tools")
  set (CPACK_PACKAGE_NAME "bs_tools")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "bs_tools")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Sabrina Krakau <sabrina.krakau@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (bs_tools "bs_tools")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:bs_tools")
