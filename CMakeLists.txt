# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /CMakeLists.txt
#
# Main CMakeLists.txt file.
# ===========================================================================

project (seqan CXX)

# ===========================================================================
# We require CMake 3.0.0. Since it is available in the currently stable
# Debian, it should be available everywhere else, too.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

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
# Build System Configuration Variables
# ===========================================================================

option (SEQAN_STATIC_APPS "Build static apps." OFF)
option (SEQAN_ARCH_SSE4 "Build SSE optimized binaries." OFF)
option (SEQAN_ARCH_AVX2 "Build AVX2 optimized binaries." OFF)
option (SEQAN_ARCH_AVX512_KNL "Build AVX512 optimized binaries for Knights Landing." OFF)
option (SEQAN_ARCH_AVX512_SKX "Build AVX512 optimized binaries for Skylake." OFF)
option (SEQAN_ARCH_AVX512_CNL "Build AVX512 optimized binaries for Cannon Lake." OFF)
option (SEQAN_ARCH_NATIVE "Build natively optimized binaries." OFF)

# ===========================================================================
# Version Check Variables.
# ===========================================================================

option (SEQAN_VERSION_CHECK_OPT_IN "Enable version check but activate as opt-in." OFF)
option (SEQAN_DISABLE_VERSION_CHECK "Disable version check." OFF) 

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
seqan_get_version ()
seqan_build_system_init ()
seqan_get_repository_info ()

# ===========================================================================
# Include Subdirectories
# ===========================================================================

# Tests, demos and manual are only built when building in DEVELOP mode.
if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    message (STATUS "Configuring tests")
    add_subdirectory (tests)

    message (STATUS "Configuring demos")
    add_subdirectory (demos)

    if (NOT SEQAN_NO_DOX)
        message (STATUS "Configuring manual")
        add_subdirectory (manual)
    endif ()
endif ()

if ((("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP") OR
     ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY")) AND
     (NOT SEQAN_NO_DOX))

    message (STATUS "Configuring dox")
    add_subdirectory (dox)

    if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
        message (STATUS "Configuring util/py_lib")
        add_subdirectory (util/py_lib)
    endif ()
endif ()

if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY")
    message (STATUS "Configuring apps")
    add_subdirectory (apps)
endif ()

# ===========================================================================
# Register Libraries (for GUI project generators)
# ===========================================================================

seqan_setup_library ()

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
