# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2012, Knut Reinert, FU Berlin
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Knut Reinert or the FU Berlin nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
# ============================================================================
#
# This CMake module will try to find SeqAn and its dependencies.  You can use
# it the same way you would use any other CMake module.
#
#   find_package(SeqAn [REQUIRED] ...)
#
# You can control the exact behaviour by setting the following variables.  The
# defaults are given after the variable name.
#
#   SEQAN_FIND_DEPENDENCIES   -- DEFAULT
#   SEQAN_FIND_ENABLE_TESTING -- TRUE if ${CMAKE_BUILD_TYPE} == "Debug", FALSE
#                                otherwise.
#
# For example:
#
#   set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
#   find_package (SeqAn)
#
# The first variable is either "ALL", "DEFAULT" or a list of dependency names
# and gives the names of the dependencies to search for.  The other two
# variables can be used to forcibly enabling/disabling the debug and testing
# mode.
#
# Valid dependencies are:
#
#   ALL     -- Forcibly enable all dependencies.
#   DEFAULT -- Enable default dependencies (zlib, OpenMP if available)
#   NONE    -- Disable all dependencies.
#
#   ZLIB    -- zlib compression library
#   BZip2   -- libbz2 compression library
#   OpenMP  -- OpenMP language extensions to C/C++
#   CUDA    -- CUDA language extensions to C/C++
#
#
# Once the search has been performed, the following variables will be set.
#
#  SEQAN_FOUND           -- Indicate whether SeqAn was found.
#
# These variables are flags that indicate whether the various dependencies
# of the SeqAn library were found.
#
#  SEQAN_HAS_ZLIB
#  SEQAN_HAS_BZIP2
#  SEQAN_HAS_OPENMP
#  SEQAN_HAS_CUDA
#
# These variables give lists that are to be passed to the
# include_directories(), target_link_libraries(), and add_definitions()
# functions.
#
#  SEQAN_INCLUDE_DIRS
#  SEQAN_LIBRARIES
#  SEQAN_DEFINITIONS
#
# Additionally, the following two variables are set.  The first contains
# the include paths for SeqAn, the second for dependencies.  This allows to
# include the dependency headers using include_directories (SYSTEM ...),
# such that warnings from these headers do not appear in the nightly builds.
#
#  SEQAN_INCLUDE_DIRS_MAIN
#  SEQAN_INCLUDE_DIRS_DEPS
#
# The C++ compiler flags to set.
#
#  SEQAN_CXX_FLAGS
#
# The following variables give the version of the SeqAn library, its
# major, minor, and the patch version part of the version string.
#
#  SEQAN_VERSION_STRING
#  SEQAN_VERSION_MAJOR
#  SEQAN_VERSION_MINOR
#  SEQAN_VERSION_PATCH
#
# ============================================================================

include(FindPackageMessage)
include(CheckIncludeFiles)

# ----------------------------------------------------------------------------
# Define Constants.
# ----------------------------------------------------------------------------

set(_SEQAN_DEFAULT_LIBRARIES ZLIB OpenMP)
set(_SEQAN_ALL_LIBRARIES     ZLIB BZip2 OpenMP CUDA)

# ----------------------------------------------------------------------------
# Set variables SEQAN_FIND_* to their default unless they have been set.
# ----------------------------------------------------------------------------

# SEQAN_FIND_DEPENDENCIES
if (NOT SEQAN_FIND_DEPENDENCIES)
  set(SEQAN_FIND_DEPENDENCIES "DEFAULT")
endif ()
if (SEQAN_FIND_DEPENDENCIES STREQUAL "DEFAULT")
  set(SEQAN_FIND_DEPENDENCIES ${_SEQAN_DEFAULT_LIBRARIES})
elseif (SEQAN_FIND_DEPENDENCIES STREQUAL "ALL")
  set(SEQAN_FIND_DEPENDENCIES ${_SEQAN_ALL_LIBRARIES})
elseif (SEQAN_FIND_DEPENDENCIES STREQUAL "NONE")
  set(SEQAN_FIND_DEPENDENCIES)
endif ()

# SEQAN_FIND_ENABLE_TESTING
if (NOT SEQAN_FIND_ENABLE_TESTING)
  set(SEQAN_FIND_ENABLE_TESTING "FALSE")
endif ()

# ----------------------------------------------------------------------------
# Compile-specific settings and workarounds around missing CMake features.
# ----------------------------------------------------------------------------

# Recognize Clang compiler.

set (COMPILER_IS_CLANG FALSE)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set (COMPILER_IS_CLANG TRUE)
endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
		
# Fix CMAKE_COMPILER_IS_GNUCXX for MinGW.

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set (CMAKE_COMPILER_IS_GNUCXX TRUE)
endif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")

# GCC Setup

if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
  # Tune warnings for GCC.
  set (CMAKE_CXX_WARNING_LEVEL 4)
  # NOTE: First location to set SEQAN_CXX_FLAGS at the moment.  If you write
  # to the variable for the first time earlier, update this line to append to
  # the variable instead of overwriting.
  set (SEQAN_CXX_FLAGS "-W -Wall -Wno-long-long -fstrict-aliasing -Wstrict-aliasing")
  set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64)

  # Determine GCC version.
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
               ARGS --version
               OUTPUT_VARIABLE __GCC_VERSION)
  # Remove all but first line.
  STRING(REGEX REPLACE "([^\n]+).*" "\\1" __GCC_VERSION ${__GCC_VERSION})
  # Find out version (3 or 2 components).
  STRING(REGEX REPLACE ".*([0-9])\\.([0-9])\\.([0-9]).*" "\\1\\2\\3"
         __GCC_VERSION ${__GCC_VERSION})
  STRING(REGEX REPLACE ".*([0-9])\\.([0-9]).*" "\\1\\20"
         _GCC_VERSION ${__GCC_VERSION})

  # Add -Wno-longlong if the GCC version is < 4.0.0.  Add -pedantic flag but
  # disable warnings for variadic macros with GCC >= 4.0.0.  Earlier versions
  # warn because of anonymous variadic macros in pedantic mode but do not have
  # a flag to disable these warnings.
  if (400 GREATER _GCC_VERSION)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-long-long")
  else (400 GREATER _GCC_VERSION)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -pedantic -Wno-variadic-macros")
  endif (400 GREATER _GCC_VERSION)
  
  # Force GCC to keep the frame pointer when debugging is enabled.  This is
  # mainly important for 64 bit but does not get into the way on 32 bit either
  # at minimal performance impact.
  if (CMAKE_BUILD_TYPE STREQUAL Debug)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer")
  elseif (CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} ${SEQAN_CXX_FLAGS_RELEASE} -g -fno-omit-frame-pointer")
  endif ()
endif ()

# Windows Setup

if (WIN32)
  # Always set NOMINMAX such that <Windows.h> does not define min/max as
  # macros.
  add_definitions (-DNOMINMAX)
endif (WIN32)

# Visual Studio Setup
if (MSVC)
  # Enable intrinics (e.g. _interlockedIncrease)
  set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /EHsc /Oi")
  # Warning level 3 for MSVC is disabled for now to see how much really bad warnings there are.
  #set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /W3)

  # TODO(holtgrew): This rather belongs into the SeqAn build system and notso much into FindSeqAn.cmake.

  # Force to always compile with W2.
  add_definitions (/W2)

  # Disable warnings about unsecure (although standard) functions.
  add_definitions (-D_SCL_SECURE_NO_WARNINGS)
endif (MSVC)

# ----------------------------------------------------------------------------
# Search for directory seqan.
# ----------------------------------------------------------------------------

option (SEQAN_USE_SEQAN_BUILD_SYSTEM "Whether or not to expect the SeqAn build system." OFF)

if (SEQAN_USE_SEQAN_BUILD_SYSTEM)
  # When using the SeqAn build system, we scan all entries in
  # CMAKE_INCLUDE_PATH for a subdirectory seqan and add all paths to the
  # variable SEQAN_INCLUDE_DIRS_MAIN.
  set (_SEQAN_INCLUDE_DIRS "")
  foreach (_SEQAN_BASEDIR ${CMAKE_INCLUDE_PATH})
    if (EXISTS ${_SEQAN_BASEDIR}/seqan)
      get_filename_component(_SEQAN_BASEDIR "${_SEQAN_BASEDIR}" ABSOLUTE)
      set(_SEQAN_INCLUDE_DIRS ${_SEQAN_INCLUDE_DIRS} ${_SEQAN_BASEDIR})
    endif (EXISTS ${_SEQAN_BASEDIR}/seqan)
  endforeach (_SEQAN_BASEDIR ${CMAKE_INCLUDE_PATH})

  if (_SEQAN_INCLUDE_DIRS)
    set(SEQAN_FOUND        TRUE)
    set(SEQAN_INCLUDE_DIRS_MAIN ${SEQAN_INCLUDE_DIRS_MAIN} ${_SEQAN_INCLUDE_DIRS})
  else (_SEQAN_INCLUDE_DIRS)
    set(SEQAN_FOUND        FALSE)
  endif (_SEQAN_INCLUDE_DIRS)
else (SEQAN_USE_SEQAN_BUILD_SYSTEM)
  # When NOT using the SeqAn build system then we only look for one directory
  # with subdirectory seqan and thus only one library.
  find_path(_SEQAN_BASEDIR "seqan" PATHS ${SEQAN_INCLUDE_PATH})
  mark_as_advanced(_SEQAN_BASEDIR)
  if (_SEQAN_BASEDIR)
    set(SEQAN_FOUND        TRUE)
    set(SEQAN_INCLUDE_DIRS_MAIN ${SEQAN_INCLUDE_DIRS_MAIN} ${_SEQAN_BASEDIR})
  else ()
    set(SEQAN_FOUND        FALSE)
  endif ()
endif (SEQAN_USE_SEQAN_BUILD_SYSTEM)

# ----------------------------------------------------------------------------
# Set defines for debug and testing.
# ----------------------------------------------------------------------------

if (SEQAN_FIND_ENABLE_TESTING)
  set(SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_ENABLE_TESTING=1)
else ()
  set(SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_ENABLE_TESTING=0)
endif ()

# ----------------------------------------------------------------------------
# Search for dependencies.
# ----------------------------------------------------------------------------

# librt, libpthread -- implicit, on Linux only

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set (SEQAN_LIBRARIES ${SEQAN_LIBRARIES} rt pthread)
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
  set (SEQAN_LIBRARIES ${SEQAN_LIBRARIES} pthread)
  set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} "-D_GLIBCXX_USE_C99=1")
endif ()

# libexecinfo -- implicit

check_include_files(execinfo.h _SEQAN_HAVE_EXECINFO)
mark_as_advanced(_SEQAN_HAVE_EXECINFO)
if (_SEQAN_HAVE_EXECINFO)
  set(SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} "-DSEQAN_HAS_EXECINFO=1")
  if (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
    set (SEQAN_LIBRARIES ${SEQAN_LIBRARIES} execinfo)
  endif (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
endif (_SEQAN_HAVE_EXECINFO)


# libstdc++ -- implicit, Mac only (clang seems not to do this automatically)

if (APPLE)
  set (SEQAN_LIBRARIES ${SEQAN_LIBRARIES} stdc++)
endif (APPLE)

# ZLIB

list(FIND SEQAN_FIND_DEPENDENCIES "ZLIB" _SEQAN_FIND_ZLIB)
mark_as_advanced(_SEQAN_FIND_ZLIB)

set (SEQAN_HAS_ZLIB FALSE)
if (NOT _SEQAN_FIND_ZLIB EQUAL -1)
  find_package(ZLIB QUIET)
  if (ZLIB_FOUND)
    set (SEQAN_HAS_ZLIB     TRUE)
    set (SEQAN_LIBRARIES         ${SEQAN_LIBRARIES}         ${ZLIB_LIBRARIES})
    set (SEQAN_INCLUDE_DIRS_DEPS ${SEQAN_INCLUDE_DIRS_DEPS} ${ZLIB_INCLUDE_DIRS})
    set (SEQAN_DEFINITIONS       ${SEQAN_DEFINITIONS}       "-DSEQAN_HAS_ZLIB=1")
  endif ()
endif ()

# BZip2

list(FIND SEQAN_FIND_DEPENDENCIES "BZip2" _SEQAN_FIND_BZIP2)
mark_as_advanced(_SEQAN_FIND_BZIP2)

set (SEQAN_HAS_BZIP2 FALSE)
if (NOT _SEQAN_FIND_BZIP2 EQUAL -1)
  find_package(BZip2 QUIET)
  if (BZIP2_FOUND)
    set (SEQAN_HAS_BZIP2    TRUE)
    set (SEQAN_LIBRARIES         ${SEQAN_LIBRARIES}         ${BZIP2_LIBRARIES})
    set (SEQAN_INCLUDE_DIRS_DEPS ${SEQAN_INCLUDE_DIRS_DEPS} ${BZIP2_INCLUDE_DIRS})
    set (SEQAN_DEFINITIONS       ${SEQAN_DEFINITIONS}       "-DSEQAN_HAS_BZIP2=1")
  endif ()
endif()

# OpenMP

list(FIND SEQAN_FIND_DEPENDENCIES "OpenMP" _SEQAN_FIND_OPENMP)
mark_as_advanced(_SEQAN_FIND_OPENMP)

set (SEQAN_HAS_OPENMP FALSE)
if (NOT _SEQAN_FIND_OPENMP EQUAL -1)
  find_package(OpenMP QUIET)
  # Note that in the following, we do not check for OPENMP_FOUND since this is
  # only true if both C and C++ compiler support OpenMP.  This is not the case
  # if the user has a compiler without OpenMP support by default and overrides
  # only the C++ compiler (e.g. on winter 2013's Mac Os X).
  if (OpenMP_CXX_FLAGS)
    set (SEQAN_HAS_OPENMP   TRUE)
    set (SEQAN_LIBRARIES         ${SEQAN_LIBRARIES}         ${OpenMP_LIBRARIES})
    set (SEQAN_INCLUDE_DIRS_DEPS ${SEQAN_INCLUDE_DIRS_DEPS} ${OpenMP_INCLUDE_DIRS})
    set (SEQAN_DEFINITIONS       ${SEQAN_DEFINITIONS}       "-DSEQAN_HAS_OPENMP=1")
    set (SEQAN_CXX_FLAGS        "${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif ()
endif ()

# CUDA

list(FIND SEQAN_FIND_DEPENDENCIES "CUDA" _SEQAN_FIND_CUDA)
mark_as_advanced(_SEQAN_FIND_CUDA)

set (SEQAN_HAS_CUDA FALSE)
if (SEQAN_ENABLE_CUDA AND NOT _SEQAN_FIND_CUDA EQUAL -1)
  find_package(CUDA QUIET)
  if (CUDA_FOUND)
    set (SEQAN_HAS_CUDA TRUE)
  endif ()
endif (SEQAN_ENABLE_CUDA AND NOT _SEQAN_FIND_CUDA EQUAL -1)

# Build SEQAN_INCLUDE_DIRS from SEQAN_INCLUDE_DIRS_MAIN and SEQAN_INCLUDE_DIRS_DEPS

set (SEQAN_INCLUDE_DIRS ${SEQAN_INCLUDE_DIRS_MAIN} ${SEQAN_INCLUDE_DIRS_DEPS})

# ----------------------------------------------------------------------------
# Determine and set SEQAN_VERSION_* variables.
# ----------------------------------------------------------------------------

if (NOT DEFINED SEQAN_VERSION_STRING)
  if (NOT CMAKE_CURRENT_LIST_DIR)  # CMAKE_CURRENT_LIST_DIR only from cmake 2.8.3.
    get_filename_component (CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
  endif (NOT CMAKE_CURRENT_LIST_DIR)

  try_run (_SEQAN_RUN_RESULT
           _SEQAN_COMPILE_RESULT
           ${CMAKE_BINARY_DIR}/CMakeFiles/SeqAnVersion
           ${CMAKE_CURRENT_LIST_DIR}/SeqAnVersion.cpp
           CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${SEQAN_INCLUDE_DIRS_MAIN}"
           COMPILE_OUTPUT_VARIABLE _COMPILE_OUTPUT
           RUN_OUTPUT_VARIABLE _RUN_OUTPUT)
  if (NOT _RUN_OUTPUT)
    message ("")
    message ("ERROR: Could not determine SeqAn version.")
    message ("COMPILE OUTPUT:")
    message (${_COMPILE_OUTPUT})
  endif (NOT _RUN_OUTPUT)

  string (REGEX REPLACE ".*SEQAN_VERSION_MAJOR:([0-9a-zA-Z]+).*" "\\1" _SEQAN_VERSION_MAJOR ${_RUN_OUTPUT})
  string (REGEX REPLACE ".*SEQAN_VERSION_MINOR:([0-9a-zA-Z]+).*" "\\1" _SEQAN_VERSION_MINOR ${_RUN_OUTPUT})
  string (REGEX REPLACE ".*SEQAN_VERSION_PATCH:([0-9a-zA-Z]+).*" "\\1" _SEQAN_VERSION_PATCH ${_RUN_OUTPUT})
  string (REGEX REPLACE ".*SEQAN_VERSION_PRE_RELEASE:([0-9a-zA-Z]+).*" "\\1" _SEQAN_VERSION_PRE_RELEASE ${_RUN_OUTPUT})

  if (SEQAN_VERSION_PRE_RELEASE EQUAL 1)
    set (_SEQAN_VERSION_DEVELOPMENT "TRUE")
  else ()
    set (_SEQAN_VERSION_DEVELOPMENT "FALSE")
  endif ()

  set (_SEQAN_VERSION_STRING "${_SEQAN_VERSION_MAJOR}.${_SEQAN_VERSION_MINOR}.${_SEQAN_VERSION_PATCH}")
  if (_SEQAN_VERSION_DEVELOPMENT)
    set (_SEQAN_VERSION_STRING "${_SEQAN_VERSION_STRING}_dev")
  endif ()

  # Cache results.
  set (SEQAN_VERSION_MAJOR "${_SEQAN_VERSION_MAJOR}" CACHE INTERNAL "SeqAn major version.")
  set (SEQAN_VERSION_MINOR "${_SEQAN_VERSION_MINOR}" CACHE INTERNAL "SeqAn minor version.")
  set (SEQAN_VERSION_PATCH "${_SEQAN_VERSION_PATCH}" CACHE INTERNAL "SeqAn patch version.")
  set (SEQAN_VERSION_PRE_RELEASE "${_SEQAN_VERSION_PRE_RELEASE}" CACHE INTERNAL "Whether version is a pre-release version version.")
  set (SEQAN_VERSION_STRING "${_SEQAN_VERSION_STRING}" CACHE INTERNAL "SeqAn version string.")

  message (STATUS "  Determined version is ${SEQAN_VERSION_STRING}")
endif (NOT DEFINED SEQAN_VERSION_STRING)

# ----------------------------------------------------------------------------
# Print Variables
# ----------------------------------------------------------------------------

if (SEQAN_FIND_DEBUG)
  message("Result for ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  message("")
  message("  CMAKE_BUILD_TYPE           ${CMAKE_BUILD_TYPE}")
  message("  CMAKE_SOURCE_DIR           ${CMAKE_SOURCE_DIR}")
  message("  CMAKE_INCLUDE_PATH         ${CMAKE_INCLUDE_PATH}")
  message("  _SEQAN_BASEDIR             ${_SEQAN_BASEDIR}")
  message("")
  message("  SEQAN_FOUND                ${SEQAN_FOUND}")
  message("  SEQAN_HAS_ZLIB             ${SEQAN_HAS_ZLIB}")
  message("  SEQAN_HAS_BZIP2            ${SEQAN_HAS_BZIP2}")
  message("  SEQAN_HAS_OPENMP           ${SEQAN_HAS_OPENMP}")
  message("  SEQAN_HAS_CUDA             ${SEQAN_HAS_CUDA}")
  message("")
  message("  SEQAN_INCLUDE_DIRS         ${SEQAN_INCLUDE_DIRS}")
  message("  SEQAN_INCLUDE_DIRS_DEPS    ${SEQAN_INCLUDE_DIRS_DEPS}")
  message("  SEQAN_INCLUDE_DIRS_MAIN    ${SEQAN_INCLUDE_DIRS_MAIN}")
  message("  SEQAN_LIBRARIES            ${SEQAN_LIBRARIES}")
  message("  SEQAN_DEFINITIONS          ${SEQAN_DEFINITIONS}")
  message("  SEQAN_CXX_FLAGS            ${SEQAN_CXX_FLAGS}")
  message("")
  message("  SEQAN_VERSION_STRING       ${SEQAN_VERSION_STRING}")
  message("  SEQAN_VERSION_MAJOR        ${SEQAN_VERSION_MAJOR}")
  message("  SEQAN_VERSION_MINORG       ${SEQAN_VERSION_MINOR}")
  message("  SEQAN_VERSION_PATCH        ${SEQAN_VERSION_PATCH}")
endif ()
