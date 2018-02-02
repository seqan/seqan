# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/insegt/CMakeLists.txt
#
# CMakeLists.txt file for insegt.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_insegt CXX)
message (STATUS "Configuring apps/insegt")

set (SEQAN_APP_VERSION "1.1.10")

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
add_executable (insegt
                insegt.cpp
                base.h
                create_gff.h
                fusion.h
                overlap_module.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (insegt ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (insegt)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install insegt in ${PREFIX}/bin directory
install (TARGETS insegt
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/insegt for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
install (FILES example/annoOutput.gff
               example/annotations.gff
               example/readOutput.gff
               example/tupleOutput.gff
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

# Add app tests if Python interpreter could be found.
if(PYTHONINTERP_FOUND)
  add_test(NAME app_test_insegt COMMAND ${PYTHON_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR})
endif(PYTHONINTERP_FOUND)

seqan_add_app_test (insegt)

# ---------------------------------------------------------------------------- 
# Setup Common Tool Description for Generic Workflow Nodes 
# ---------------------------------------------------------------------------- 
 
# Include executable razers in CTD structure. 
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES} insegt CACHE INTERNAL "") 

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:insegt")
  set (CPACK_PACKAGE_NAME "insegt")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "insegt - Alignment Free Sequence Comparison")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")

  seqan_configure_cpack_app (insegt "insegt")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:insegt")

