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
# Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
# ============================================================================
# Support for the SeqAn legacy build system.
# ============================================================================

# ---------------------------------------------------------------------------
# Initialize System
#
# We set the variable SEQAN_BUILD_SYSTEM_LEGACY to OFF by default.  If legacy
# functions are used then this is set to ON.  In this case, we will print a
# warning message in seqan_legacy_warning () which is called by the master
# CMakeLists.txt to inform users that the support for the legacy build system
# will eventually go away.
# ---------------------------------------------------------------------------

if (NOT SEQAN_BUILD_SYSTEM_LEGACY)
  set (SEQAN_BUILD_SYSTEM_LEGACY OFF CACHE INTERNAL "Used legacy build system.")
endif (NOT SEQAN_BUILD_SYSTEM_LEGACY)

# ---------------------------------------------------------------------------
# Macro seqan_setup_includes ()
#
# Called in /sandbox/${NAME}/CMakeLists.txt.  Translate to the content of the
# new-style CMakeLists.txt there.
# ---------------------------------------------------------------------------

macro (seqan_setup_includes PATH TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")

  string (REPLACE "SeqAn" "" _PART_NAME "${TARGET_NAME}")

  # Add the paths core/include and extras/include to the paths that CMake
  # searches for libraries.
  set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include
                                                ${CMAKE_CURRENT_SOURCE_DIR}/../../core/include
                                                ${CMAKE_CURRENT_SOURCE_DIR}/../../extras/include)

  # Setup the library modules from extras.
  seqan_setup_library ("${_PART_NAME}" extras core)

  # Add subdirectory for apps.
  add_subdirectory (apps)

  # Demos are required when doing a Whole SeqAn Release (copy demos) or when
  # developing (build demos).
  if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP"))
    add_subdirectory (demos)
  endif ()

  # Tests are only built when building in DEVLOP mode.
  if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    add_subdirectory (tests)
  endif ()
endmacro (seqan_setup_includes PATH TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_make_seqan_available (TARGET_NAME)
#
# Used to register a part of the library.  Does nothing now.
# ---------------------------------------------------------------------------

macro (seqan_make_seqan_available TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
endmacro (seqan_make_seqan_available TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_tests (TARGET_NAME)
#
# Forward to seqan_register_tests().
# ---------------------------------------------------------------------------

macro (seqan_setup_tests)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
  seqan_register_tests ()
endmacro (seqan_setup_tests)

# ---------------------------------------------------------------------------
# Macro seqan_add_test_executable (NAME [SOURCES])
# ---------------------------------------------------------------------------

macro (seqan_add_test_executable NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")

  # ----------------------------------------------------------------------------
  # Dependencies
  # ----------------------------------------------------------------------------

  # Search SeqAn and select dependencies.
  set (SEQAN_FIND_DEPENDENCIES NONE)
  find_package (SeqAn REQUIRED)

  # ----------------------------------------------------------------------------
  # Build Setup
  # ----------------------------------------------------------------------------

  # Add include directories.
  include_directories (${SEQAN_INCLUDE_DIRS})

  # Add definitions set by find_package(SeqAn).
  add_definitions (${SEQAN_DEFINITIONS})

  # Update the list of file names below if you add source files to your test.
  add_executable(${NAME} ${ARGN})

  # Add dependencies found by find_package(SeqAn).
  target_link_libraries (${NAME} ${SEQAN_LIBRARIES})

  # Add CXX flags found by find_package(SeqAn).
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS})

  # ----------------------------------------------------------------------------
  # Register with CTest
  # ----------------------------------------------------------------------------

  add_test(NAME test_${NAME} COMMAND $<TARGET_FILE:${NAME}>)
endmacro (seqan_add_test_executable NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_apps (TARGET_NAME)
#
# Used to setup the demo directory.  Now only forwards to the macro
# seqan_register_apps ().
# ---------------------------------------------------------------------------

macro (seqan_setup_apps TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")

  seqan_register_apps ()
endmacro (seqan_setup_apps TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_demos (TARGET_NAME)
#
# Used to setup the demo directory.  Now only forwards to the macro
# seqan_register_demos ().
# ---------------------------------------------------------------------------

macro (seqan_setup_demos TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")

  seqan_register_demos ("${TARGET_NAME}_")
endmacro (seqan_setup_demos TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_executables (TARGET_NAME)
# ---------------------------------------------------------------------------

macro (seqan_add_all_executables TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
endmacro (seqan_add_all_executables TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_executable (TARGET_NAME)
#
# Used to be used to add executables.  Was only used directly for apps.
#
# Replaced by most of the contents of an apps CMakeLists.txt.
# ---------------------------------------------------------------------------

macro (seqan_add_executable TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
  if (NOT SEQAN_FIND_CALLED)
    set (SEQAN_FIND_CALLED ON)

    # ----------------------------------------------------------------------------
    # Dependencies
    # ----------------------------------------------------------------------------

    # Search SeqAn and select dependencies.
    set (SEQAN_FIND_DEPENDENCIES NONE)
    find_package (SeqAn REQUIRED)

    # ----------------------------------------------------------------------------
    # Build Setup
    # ----------------------------------------------------------------------------

    # Add include directories.
    include_directories (${SEQAN_INCLUDE_DIRS})

    # Add definitions set by find_package(SeqAn).
    add_definitions (${SEQAN_DEFINITIONS})

    # Add CXX flags found by find_package(SeqAn).
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
  endif (NOT SEQAN_FIND_CALLED)

  # Update the list of file names below if you add source files to your application.
  add_executable(${TARGET_NAME} ${ARGN})

  # Add dependencies found by find_package(SeqAn).
  target_link_libraries (${TARGET_NAME} ${SEQAN_LIBRARIES})
endmacro (seqan_add_executable TARGET_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_add_cuda_executable ()
#
# Used to add a CUDA executable.  Was only used directly for tests.
#
# Replaced by most of the contents of a CUDA test CMakeLists.txt.
# ---------------------------------------------------------------------------

macro (seqan_add_cuda_executable TARGET_NAME)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
  if (NOT SEQAN_FIND_CALLED)
    set (SEQAN_FIND_CALLED ON)

    # ----------------------------------------------------------------------------
    # Dependencies
    # ----------------------------------------------------------------------------

    # Search SeqAn and select dependencies.
    set (SEQAN_FIND_DEPENDENCIES NONE)
    find_package (SeqAn REQUIRED)

    # ----------------------------------------------------------------------------
    # Build Setup
    # ----------------------------------------------------------------------------

    # Add include directories.
    include_directories (${SEQAN_INCLUDE_DIRS})

    # Add definitions set by find_package(SeqAn).
    add_definitions (${SEQAN_DEFINITIONS})

    # Add CXX flags found by find_package(SeqAn).
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
  endif (NOT SEQAN_FIND_CALLED)

  if (NOT CUDA_FIND_CALLED)
    set (CUDA_FIND_CALLED ON)
    # Search for CUDA.
    find_package (CUDA)

    # Set CUDA options.
    set (CUDA_PROPAGATE_HOST_FLAGS OFF)
    set (CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
    # Remove -pedantic flag.
    string (REGEX REPLACE "\\-pedantic" ""
            CUDA_CXX_FLAGS ${CUDA_NVCC_FLAGS} ${CMAKE_CXX_FLAGS})

    # Enable .cu as a CXX source file extension for linking.
    list (APPEND CMAKE_CXX_SOURCE_FILE_EXTENSIONS "cu")
    # Add CUT include directories for CUDA.
    cuda_include_directories(${CUDA_CUT_INCLUDE_DIR})
  endif (NOT CUDA_FIND_CALLED)

  if (NOT CUDA_FOUND)
    message (STATUS "CUDA not found.  Not building ${TARGET_NAME}")
  else (NOT CUDA_FOUND)
    # Update the list of file names below if you add source files to your application.
    cuda_add_executable(${TARGET_NAME} ${ARGN})

    # Add dependencies found by find_package(SeqAn).
    target_link_libraries (${TARGET_NAME} ${SEQAN_LIBRARIES})
  endif (NOT CUDA_FOUND)
endmacro (seqan_add_cuda_executable TARGET_NAME)


# ---------------------------------------------------------------------------
# Macro seqan_add_cuda_test_executable ()
#
# Used to add a CUDA test executable.  Was only used directly for tests.
#
# Replaced by most of the contents of a CUDA test CMakeLists.txt.
# ---------------------------------------------------------------------------

macro (seqan_add_cuda_test_executable)
  seqan_add_cuda_executable (${ARGV})

  if (CUDA_FOUND)
    # Add just added CUDA executable as a test.
    add_test (test_${TARGET_NAME} ${TARGET_NAME})
  endif (CUDA_FOUND)
endmacro (seqan_add_cuda_test_executable)

# ---------------------------------------------------------------------------
# Macro seqan_add_all_subdirectories ()
#
# Nothing.
# ---------------------------------------------------------------------------

macro (seqan_add_all_subdirectories)
  set (SEQAN_BUILD_SYSTEM_LEGACY ON CACHE INTERNAL "Used legacy build system")
endmacro (seqan_add_all_subdirectories)

# ---------------------------------------------------------------------------
# macro seqan_legacy_warning ()
#
# Print warning when using legacy build system if it has been used.
# ---------------------------------------------------------------------------

function (seqan_legacy_warning)
  if (NOT SEQAN_BUILD_SYSTEM_LEGACY)
    return ()  # Do not print if not used.
  endif ()

  message ("  ***************************************************************************")
  message ("  * /!\\ NOTICE: YOU HAVE USED THE THE OLD SEQAN BUILD SYSTEMS FUNCTIONS /!\\ *")
  message ("  ***************************************************************************")
  message ("")
  message ("    For now, we still support the old build system such that you can")
  message ("    migrate to the new system on a slower upgrade path.  However, please")
  message ("    do so in the near future since support for the new build system will")
  message ("    stop before the next (1.4) release in early 2013.")
  message ("")
  message ("  ***************************************************************************")
  message ("  * /!\\ NOTICE: YOU HAVE USED THE THE OLD SEQAN BUILD SYSTEMS FUNCTIONS /!\\*")
  message ("  ***************************************************************************")
endfunction (seqan_legacy_warning)
