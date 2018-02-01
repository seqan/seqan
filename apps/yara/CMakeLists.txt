# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/yara/CMakeLists.txt
#
# CMakeLists.txt file for yara.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_yara CXX)
message (STATUS "Configuring apps/yara")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

if (COMPILER_MSVC)
    message (STATUS "Yara is currently not build on Microsoft Visual because it requires OpenMP 3.0.")
    return ()
endif ()

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  OpenMP 3.0 not found: building yara without multi-threading.")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Zlib not found: building yara without gzipped input and BAM output.")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (STATUS "  BZip2 not found: building yara without bzipped input.")
endif (NOT BZIP2_FOUND)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "0.9.11")

option (YARA_LARGE_CONTIGS "Set to OFF to disable support for more than 32k contigs or contigs longer than 4Gbp." ON)
if (YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)
  add_definitions (-DYARA_LARGE_CONTIGS=1)
endif (YARA_LARGE_CONTIGS AND NOT SEQAN_TRAVIS_BUILD)

# ----------------------------------------------------------------------------
# Disable Compiler Warnings
# ----------------------------------------------------------------------------

# Disables warning about conversion from 64 bit to 32 bit integrals.
# Note that there is a runtime exception thrown if the type is not sufficient.
# Unfortunately the VC does not detect it.
if (MSVC)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4244")  # When compiling on 32 bit windows.
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4267")  # When compiling on 64 bit windows.
endif (MSVC)

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
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(yara_indexer indexer.cpp
                            misc_tags.h
                            misc_timer.h
                            misc_options.h
                            misc_types.h
                            store_seqs.h
                            bits_matches.h
                            index_fm.h)

add_executable(yara_mapper    mapper.cpp
                              mapper.h
                              mapper_collector.h
                              mapper_classifier.h
                              mapper_ranker.h
                              mapper_filter.h
                              mapper_extender.h
                              mapper_verifier.h
                              mapper_aligner.h
                              mapper_writer.h
                              basic_alphabet.h
                              file_pair.h
                              file_prefetched.h
                              store_seqs.h
                              misc_tags.h
                              misc_timer.h
                              misc_options.h
                              misc_types.h
                              bits_reads.h
                              bits_hits.h
                              bits_matches.h
                              bits_context.h
                              bits_seeds.h
                              bits_bucket.h
                              find_extender.h
                              find_verifier.h
                              index_fm.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (yara_indexer ${SEQAN_LIBRARIES})
target_link_libraries (yara_mapper ${SEQAN_LIBRARIES})

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (yara)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install yara in ${PREFIX}/bin directory
install (TARGETS yara_indexer yara_mapper
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/yara for SeqAn release builds.
install (FILES LICENSE README.rst
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (yara)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable yara_mapper in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES}
                           yara_indexer yara_mapper CACHE INTERNAL "")

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
  set (CPACK_PACKAGE_NAME "yara")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Yara - Yet Another Read Aligner")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
  set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")

  seqan_configure_cpack_app (yara "yara")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
