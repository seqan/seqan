# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /CMakeLists.txt
#
# Main CMakeLists.txt file.
# ===========================================================================

project (seqan)

# ===========================================================================
# We require CMake 2.8.2.  Since it is available in the currently stable
# Debian, it should be available everywhere else, too.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)

# ===========================================================================
# Iff the file INST or SUA exists in the source root
# turn on SeqAn instrumentation for experimental data collection.
# ===========================================================================

if (EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)
  set (SEQAN_INSTRUMENTATION YES)
else (EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)
  set (SEQAN_INSTRUMENTATION NO)
endif (EXISTS ${CMAKE_SOURCE_DIR}/INST OR EXISTS ${CMAKE_SOURCE_DIR}/SUA OR EXISTS ${CMAKE_SOURCE_DIR}/../SUA)


# ===========================================================================
# Warn People About the Perils Of In-Source Builds
# ===========================================================================

if (EXISTS "$CMAKE_BINARY_DIR/CMakeLists.txt")
  if (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
    set (IN_SOURCE_BUILD YES)
  else (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
    message ("CAUTION: In-source builds are not allowed by default!")
    message ("")
    message ("(If you want to use Eclipse and have read the 'Getting Started,")
    message ("with Eclipse Tutorial' then see the instructions below.)")
    message ("")
    message ("Instead of doing in source builds, go to into the directory build")
    message ("and call CMake from there.")
    message ("We recommend you to create an additional directory level, this")
    message ("allows you to have multiple build configurations (e.g. Debug,")
    message ("Release) or project types.")
    message ("")
    message ("  cd build")
    message ("  mkdir Debug")
    message ("  cd Debug")
    message ("  cmake ../.. -DCMAKE_BUILD_TYPE=Debug")
    message ("")
    message ("NOTE: CMake has created some files in this directory that you have")
    message ("to remove before executing the above commands.")
    message ("")
    message ("  [Linux, Mac Os X]  rm -rf CMakeCache.txt CMakeFiles")
    message ("  [Windows]          del CMakeCache.txt ; rmdir /S CMakeFiles")
    message ("")
    message ("If you want to use Eclipse then you might want to have in-source")
    message ("builds because of limitations of Eclipse CDT4.  In this case")
    message ("Create a file in this directory called I_WANT_IN_SOURCE_BUILDS:")
    message ("")
    message ("  [Linux, Mac Os X]  touch I_WANT_IN_SOURCE_BUILDS")
    message ("  [Windows]          type nul >> I_WANT_IN_SOURCE_BUILDS")
    message ("")
    message (FATAL_ERROR "Please follow the instructions above.")
  endif (EXISTS "$CMAKE_BINARY_DIR/I_WANT_IN_SOURCE_BUILDS")
endif (EXISTS "$CMAKE_BINARY_DIR/CMakeLists.txt")

# ===========================================================================
# Setup Variables
# ===========================================================================

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo ..."
         FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Clear CTD executable list, will be set on reconfiguration.
set (SEQAN_CTD_EXECUTABLES CACHE INTERNAL "Global list of executables for workflow integration.")

# ===========================================================================
# Setup Modules and Find Packages.
# ===========================================================================

# Note that the find modules in util/make take precedence
# over the built-in ones from CMake.  This is used in
# FindBZip2, for example which we slightly modify.
set (CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/FindTBB
    ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake
    ${CMAKE_MODULE_PATH})

#message (STATUS "CMAKE_MODULE_PATH is ${CMAKE_MODULE_PATH}")

# ===========================================================================
# Add contrib to CMAKE_FIND_ROOT_PATH.
# ===========================================================================

include (SeqAnContribs)

# ===========================================================================
# Build Mode Branches
# ===========================================================================

message (STATUS "Initializing SeqAn Build System...")

set (SEQAN_USE_SEQAN_BUILD_SYSTEM TRUE CACHE INTERNAL "Use SeqAn build system." FORCE)

#set (SEQAN_PREFIX_SHARE .)
include (SeqAnBuildSystem)
set (SEQAN_ROOT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Root source directory." FORCE)
seqan_build_system_init ()
seqan_get_repository_info ()

# ===========================================================================
# Include Subdirectories
# ===========================================================================

message (STATUS "Configuring core")
add_subdirectory (core)
message (STATUS "Configuring extras")
add_subdirectory (extras)
message (STATUS "Configuring docs")
add_subdirectory (docs)
message (STATUS "Configuring sandbox")
add_subdirectory (sandbox)

# ===========================================================================
# Include CPack
# ===========================================================================

include (package)

# ===========================================================================
# Check Optional Libraries, Print Usage
# ===========================================================================

#include (FeatureSummary)

# ===========================================================================
# Setup Common Tool Description for Generic Workflow Nodes
# ===========================================================================

include (SeqAnCtdSetup)
