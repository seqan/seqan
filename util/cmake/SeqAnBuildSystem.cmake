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
# This CMake file defines the necessary macros for the SeqAn build system.
#
# Note that while the SeqAn build system uses the FindSeqAn.cmake module,
# the FindSeqAn.cmake module itself can be used independently from the SeqAn
# build system.
# ============================================================================

# Valid values for SEQAN_BUILD_SYSTEM:
#
# DEVELOP
# SEQAN_RELEASE
# SEQAN_RELEASE_LIBRARY
# APP:${app_name}

include (SeqAnUsabilityAnalyzer)

# ---------------------------------------------------------------------------
# Normalize CMAKE_CXX_FLAGS to be a string.
#
# If we do not do this then setting the environment variable CXXFLAGS will
# cause CMAKE_CXX_FLAGS to become a list and this will generate compiler
# command lines including the list item separator semicolon ";".  This makes
# the compiler command fail.
# ---------------------------------------------------------------------------

#if (CMAKE_CXX_FLAGS)
#  foreach (_FLAG ${CMAKE_CXX_FLAGS})
#    set (_FLAGS "${_FLAGS} ${_FLAG}")
#  endforeach (_FLAG ${CMAKE_CXX_FLAGS})
#  set (CMAKE_CXX_FLAGS "${_FLAGS}")
#endif (CMAKE_CXX_FLAGS)

# ---------------------------------------------------------------------------
# Enable /bigobj flag on Windows.
# ---------------------------------------------------------------------------

# We need the /bigobj switch on windows (for 64 bit builds only actually).
# Set target system to be Windows Vista and later.
if (MSVC)
  add_definitions (/bigobj /D_WIN32_WINNT=0x0600 /DWINVER=0x0600)
elseif (MINGW)
  add_definitions (-D_WIN32_WINNT=0x0600 -DWINVER=0x0600)
endif (MSVC)

# ---------------------------------------------------------------------------
# Set architecture for MinGW.
#
# If we do not set i586 as the architecture for MinGW then generating atomic
# expressions will fail.
# ---------------------------------------------------------------------------

if (MINGW)
	if ("${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
	    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=i586")
	endif ("${CMAKE_SIZEOF_VOID_P}" EQUAL "4")
endif (MINGW)

# ---------------------------------------------------------------------------
# Disable false positive terminal detection in Xcode
# ---------------------------------------------------------------------------

if (CMAKE_GENERATOR STREQUAL Xcode)
  add_definitions (-DSEQAN_NO_TERMINAL)
endif (CMAKE_GENERATOR STREQUAL Xcode)

# ---------------------------------------------------------------------------
# Function add_executable (name [WIN32] [MACOSX_BUNDLE] [EXCLUDE_FROM_ALL]
#                          source1 source2 ... sourceN)
#
# Add an executable with the given name and sources.
#
# We overwrite the built-in function add_executable to automatically add
# behaviour that is required in the SeqAn build system.  This includes:
#
# * Adding dependencies to the SeqAn library headers.
# * Enabling the SeqAn Usability Analyzer (SUA).
#
# Note that it is not possible to overwrite the same function two times.
# ---------------------------------------------------------------------------

function (add_executable NAME)
    # Call overwritten _add_executable.
    _add_executable(${ARGV})

    # Add dependencies on the SeqAn library.
    add_dependencies(${NAME} seqan_library)

    # Add dependency on the SUA target.
    seqan_add_sua_dependency (${NAME})
endfunction (add_executable)

# ---------------------------------------------------------------------------
# Macro seqan_add_app_subdirectory (APP_NAME)
# ---------------------------------------------------------------------------

macro (seqan_add_app_subdirectory APP_NAME)
    if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
         "${SEQAN_BUILD_SYSTEM}" STREQUAL "APP:${APP_NAME}")
        add_subdirectory (${APP_NAME})
    endif ()
endmacro (seqan_add_app_subdirectory)

# ---------------------------------------------------------------------------
# Macro seqan_register_apps ()
#
# Register all apps by adding their subdirectories if they are to be built
# (SEQAN_RELEASE and APP:${app} modes).
# ---------------------------------------------------------------------------

macro (seqan_register_apps)
    # Set SeqAn flags.
    set (SEQAN_FIND_ENABLE_TESTING 0)
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

    # enable static linkage for seqan apps
    if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG AND NOT MINGW)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
      set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
    endif ()

    # Get all direct entries of the current source directory into ENTRIES.
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)

    # Add all values from ${ENTRIES} that are subdirectories and have a file
    # CMakeListst.txt.
    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
                    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP") OR
                    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_APPS") OR
                    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "APP:${ENTRY}"))
                    add_subdirectory(${ENTRY})
                endif ()
            endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
        endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_register_apps)

# ---------------------------------------------------------------------------
# Macro seqan_build_system_init ()
#
# Initialize build system.
# ---------------------------------------------------------------------------

macro (seqan_build_system_init)
    # Enable CTest and command add_test().
    enable_testing ()

    if (NOT SEQAN_BUILD_SYSTEM)
        set (SEQAN_BUILD_SYSTEM "DEVELOP")
    endif (NOT SEQAN_BUILD_SYSTEM)
    set (SEQAN_APP_VERSION "0.0.0" CACHE STRING "Version of the application.")
    set (SEQAN_NIGHTLY_RELEASE FALSE CACHE BOOL "Set to TRUE to enable nightly app releases.")

    if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
        ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY"))
        # Install SeqAn LICENSE, README.rst, CHANGELOG.rst files.
        install (FILES LICENSE
                       README.rst
                       CHANGELOG.rst
                 DESTINATION share/doc/seqan)
    endif ()

    set (SEQAN_BUILD_SYSTEM "DEVELOP" CACHE STRING "Build/Release mode to select. One of DEVELOP SEQAN_RELEASE, APP:\${APP_NAME}. Defaults to DEVELOP.")

    set (_CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set (CMAKE_INCLUDE_PATH ${_CMAKE_INCLUDE_PATH} CACHE STRING "")

    SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
         ${PROJECT_BINARY_DIR}/bin)
endmacro (seqan_build_system_init)

# ---------------------------------------------------------------------------
# Macro seqan_add_app_test (APP_NAME SUFFIX)
#
# Add app test invocation.
# ---------------------------------------------------------------------------

# App tests are run using Python.  Search for Python and register test if the
# Python interpreter could be found.

macro (seqan_add_app_test APP_NAME)
    if (MODEL MATCHES ".*MemCheck.*")
        set (_VALGRIND_FLAG --valgrind)
    else ()
        set (_VALGRIND_FLAG)
    endif ()
    find_package (PythonInterp)
    if (PYTHONINTERP_FOUND)
      add_test (NAME app_test_${APP_NAME}${ARGV1}
                COMMAND ${PYTHON_EXECUTABLE}
                        ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests${ARGV1}.py
                        ${_VALGRIND_FLAG}
                        ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
    endif (PYTHONINTERP_FOUND)
endmacro (seqan_add_app_test APP_NAME)

# ---------------------------------------------------------------------------
# Macro seqan_setup_library ()
#
# * Creates install targets for the library.
# * Writes list SeqAn headers to ${_SEQAN_HEADERS}
# ---------------------------------------------------------------------------

macro (seqan_setup_library)
    # Only install the library if the virtual build packages "SEQAN_RELEASE"
    # or "SEQAN_LIBRARY_ONLY" are chosen.
    if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
        ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY"))
        file (GLOB HEADERS
              RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
              include/seqan/[A-z]*/[A-z]/[A-z]*.h
              include/seqan/[A-z]*/[A-z]*.h
              include/seqan/[A-z]*.h)
        foreach (HEADER ${HEADERS})
            get_filename_component (_DESTINATION ${HEADER} PATH)
            install (FILES ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER} DESTINATION ${_DESTINATION})
        endforeach ()
    endif ()

    # Get list of header and super header files.
    file (GLOB SUPER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/[A-z]*.h)
    file (GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/[A-z]*/[A-z]*.h)
    file (GLOB SUB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan/[A-z]*/[A-z]*/[A-z]*.h)

    # Sort headers for Xcode, ...
    if (SUB_HEADERS)
        list (SORT SUB_HEADERS)
    endif (SUB_HEADERS)

    # Sort headers for Xcode, ...
    if (HEADERS)
        list (SORT HEADERS)
    endif (HEADERS)

    # Sort super-headers for Xcode, ...
    if (SUPER_HEADERS)
        list (SORT SUPER_HEADERS)
    endif (SUPER_HEADERS)

    # Create source groups for Visual Studio (and possibly other IDEs).
    foreach (HEADER ${HEADERS})
        file (RELATIVE_PATH HEADER_REL ${CMAKE_CURRENT_SOURCE_DIR}/include/seqan ${HEADER})
        get_filename_component (MODULE ${HEADER_REL} PATH)
        source_group (seqan\\${MODULE} FILES ${HEADER})
    endforeach (HEADER ${HEADERS})
    source_group (seqan FILES ${SUPER_HEADERS})

#    # CMake bug workaround: For Non-IDE generators there is a bug in cmake.
#    # The SOURCE command in add_custom_target is not recognized there.
#    set (NONIDE_GENERATORS "Unix Makefiles" "MinGW Makefiles")
#    list (FIND NONIDE_GENERATORS ${CMAKE_GENERATOR} FOUND)
#    if (FOUND EQUAL -1)
#        set (IDE_SOURCES SOURCES ${HEADERS} ${SUPER_HEADERS})
#    endif (FOUND EQUAL -1)

    # Add pseudo target for the library part.  Note that the IDE_SOURCES
    # variable includes the "SOURCES" argument for add_custom_target when
    # building with a generator for an IDE.
    add_custom_target (seqan_library SOURCES ${SUB_HEADERS} ${HEADERS} ${SUPER_HEADERS})
endmacro (seqan_setup_library)

# ---------------------------------------------------------------------------
# Macro seqan_setup_install_vars (APP_NAME)
#
# Setup variables for install, depending on build mode.
# ---------------------------------------------------------------------------

macro (seqan_setup_install_vars APP_NAME)
    if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "APP:${APP_NAME}")
        set (SEQAN_PREFIX_SHARE ".")
        set (SEQAN_PREFIX_SHARE_DOC ".")
    else ()
        set (SEQAN_PREFIX_SHARE "share/${APP_NAME}")
        set (SEQAN_PREFIX_SHARE_DOC "share/doc/${APP_NAME}")
    endif ()
endmacro (seqan_setup_install_vars)

# ---------------------------------------------------------------------------
# Macro seqan_configure_cpack_app (APP_NAME APP_DIR)
#
# Setup variables for install, depending on build mode.
#
# Sets defaults for CPACK_PACKAGE_DESCRIPTION_FILE and CPACK_RESOURCE_FILE_LICENSE
# ---------------------------------------------------------------------------

macro (seqan_configure_cpack_app APP_NAME APP_DIR)
  # The following include automates the MS Redistributable installer.
  include (InstallRequiredSystemLibraries)

  if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set (CPACK_GENERATOR "ZIP")
  else ()
    set (CPACK_GENERATOR "ZIP;TBZ2")
  endif ()

  # Set defaults for CPACK_PACKAGE_DESCRIPTION_FILE and CPACK_RESOURCE_FILE_LICENSE
  if (NOT CPACK_PACKAGE_DESCRIPTION_FILE)
    set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
  endif ()
  if (NOT CPACK_RESOURCE_FILE_LICENSE)
    set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
  endif ()

  # Automatically deduce system name for CPack.
  include (SetCPackSystemName)

  # Get SEQAN_APP_VERSION_{MAJOR,MINOR,PATCH} from SEQAN_APP_VERSION.
  set (SEQAN_APP_VERSION_MAJOR "0")
  if (SEQAN_APP_VERSION MATCHES "^([0-9]+).*")
    string (REGEX REPLACE "^([0-9]+).*" "\\1" _SEQAN_APP_VERSION_MAJOR "${SEQAN_APP_VERSION}")
  endif ()
  if (_SEQAN_APP_VERSION_MAJOR)
    set(SEQAN_APP_VERSION_MAJOR "${_SEQAN_APP_VERSION_MAJOR}")
  endif ()
  set (SEQAN_APP_VERSION_MINOR "0")
  if (SEQAN_APP_VERSION MATCHES "^[0-9]+\\.([0-9]+).*")
    string (REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" _SEQAN_APP_VERSION_MINOR "${SEQAN_APP_VERSION}")
  endif ()
  if (_SEQAN_APP_VERSION_MINOR)
    set(SEQAN_APP_VERSION_MINOR "${_SEQAN_APP_VERSION_MINOR}")
  endif ()
  set (SEQAN_APP_VERSION_PATCH "0")
  if (SEQAN_APP_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.([0-9]+)$")
    string (REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)$" "\\1" _SEQAN_APP_VERSION_PATCH "${SEQAN_APP_VERSION}")
  endif ()
  if (_SEQAN_APP_VERSION_PATCH)
    set(SEQAN_APP_VERSION_PATCH "${_SEQAN_APP_VERSION_PATCH}")
  endif ()

  # Setup the app version.  SEQAN_APP_VERSION_{MAJOR,MINOR,PATCH} have
  # to be set.  To create nightly releases, set SEQAN_NIGHTLY_RELEASE to
  # TRUE on the command line.
  if (SEQAN_NIGHTLY_RELEASE)
    include (GetCurrentDate)
    set (CPACK_PACKAGE_VERSION "${CURRENT_YEAR}${CURRENT_MONTH}${CURRENT_DAY}")
  else ()
    set (CPACK_PACKAGE_VERSION "${SEQAN_APP_VERSION_MAJOR}.${SEQAN_APP_VERSION_MINOR}.${SEQAN_APP_VERSION_PATCH}")
  endif ()
  set (CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_APP_VERSION_MAJOR}")
  set (CPACK_PACKAGE_VERSION_MINOR "${SEQAN_APP_VERSION_MINOR}")
  set (CPACK_PACKAGE_VERSION_PATCH "${SEQAN_APP_VERSION_PATCH}")

  set (CPACK_PACKAGE_INSTALL_DIRECTORY "${APP_DIR} ${CPACK_PACKAGE_VERSION}")

  include (CPack)
endmacro (seqan_configure_cpack_app)


# ---------------------------------------------------------------------------
# Macro seqan_setup_cuda_vars ([DISABLE_WARNINGS] [DEBUG_DEVICE]
#                              [ARCH sm_xx] [FLAGS flags ...])
#
# Setup CUDA variables.
# ---------------------------------------------------------------------------

macro (seqan_setup_cuda_vars)
  cmake_parse_arguments(_SEQAN_CUDA
                        "DISABLE_WARNINGS;DEBUG_DEVICE"
                        "ARCH"
                        "FLAGS"
                        ${ARGN})
  if (SEQAN_HAS_CUDA)
    # Wrap nvcc to make cudafe output gcc-like.
    find_program (COLOR_NVCC colornvcc PATHS ${CMAKE_SOURCE_DIR}/util NO_DEFAULT_PATH)
    set (CUDA_NVCC_EXECUTABLE ${COLOR_NVCC})

    # Build CUDA targets from the given architecture upwards.
    set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch ${_SEQAN_CUDA_ARCH} ${_SEQAN_CUDA_FLAGS}")

    # Add debug symbols to device code.
    if (_SEQAN_CUDA_DISABLE_WARNINGS)
      set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -G")
    endif ()

    # Add flags for the CUDA compiler.
    list (APPEND CUDA_NVCC_FLAGS_RELEASE "-O3")
    list (APPEND CUDA_NVCC_FLAGS_MINSIZEREL "-O3")
    list (APPEND CUDA_NVCC_FLAGS_RELWITHDEBINFO "-O3 -g -lineinfo")
    list (APPEND CUDA_NVCC_FLAGS_DEBUG "-O0 -g -lineinfo")

    if (_SEQAN_CUDA_DISABLE_WARNINGS)
      # Disable all CUDA warnings.
      set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --disable-warnings")
    else ()
      # Disable only Thrust warnings.
      string (REGEX REPLACE "-Wall" ""
              SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS}")
      string (REGEX REPLACE "-pedantic" ""
              SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS}")
      if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-unused-parameter")
      endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif ()

    # Fix CUDA on OSX.
    if (APPLE AND COMPILER_IS_CLANG)
      # (weese:) I had to deactivate the C compiler override to make it compile again
      # NVCC mistakes /usr/bin/cc as gcc.
      #list (APPEND CUDA_NVCC_FLAGS "-ccbin /usr/bin/clang")
      # NVCC does not support libc++.
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libstdc++")
    endif ()
  endif ()
endmacro (seqan_setup_cuda_vars)


# ---------------------------------------------------------------------------
# Function seqan_get_version()
#
# Sets the variables SEQAN_VERSION, SEQAN_VERSION_MAJOR, SEQAN_VERSION_MINOR,
# SEQAN_VERSION_PATCH, determined from seqan/version.h
# ---------------------------------------------------------------------------

macro (seqan_get_version)
  try_run(_SEQAN_RUN_RESULT
          _SEQAN_COMPILE_RESULT
          ${CMAKE_BINARY_DIR}/CMakeFiles/SeqAnVersion
          ${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/SeqAnVersion.cpp
          CMAKE_FLAGS -DINCLUDE_DIRECTORIES:STRING=${CMAKE_CURRENT_SOURCE_DIR}/include
          COMPILE_OUTPUT_VARIABLE _COMPILE_OUTPUT
          RUN_OUTPUT_VARIABLE _RUN_OUTPUT)
  if (NOT _RUN_OUTPUT)
    message("")
    message("ERROR: Could not determine SeqAn version.")
    message("COMPILE OUTPUT:")
    message(${_COMPILE_OUTPUT})
  endif (NOT _RUN_OUTPUT)
  string(REGEX REPLACE ".*SEQAN_VERSION_MAJOR:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_MAJOR ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_MINOR:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_MINOR ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_PATCH:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_PATCH ${_RUN_OUTPUT})
  string(REGEX REPLACE ".*SEQAN_VERSION_PRE_RELEASE:([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_PRE_RELEASE ${_RUN_OUTPUT})
  set(SEQAN_VERSION "${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}")
#  if (SEQAN_VERSION_PRE_RELEASE STREQUAL 1)
#    set(SEQAN_VERSION "pre${SEQAN_VERSION}")
#  endif (SEQAN_VERSION_PRE_RELEASE STREQUAL 1)
endmacro (seqan_get_version)

# ---------------------------------------------------------------------------
# Function seqan_get_repository_info()
#
# Sets the variables SEQAN_DATE and SEQAN_REVISION determined from git.
# ---------------------------------------------------------------------------

macro (seqan_get_repository_info)
  set (_SEQAN_GIT_DIR "${CMAKE_SOURCE_DIR}/.git")

  # Get Git information.
  if (EXISTS ${_SEQAN_GIT_DIR})
    find_package (GitInfo QUIET)
    if (GIT_FOUND)
      GIT_WC_INFO (${CMAKE_SOURCE_DIR} _SEQAN)
    endif ()
  else ()
    message(STATUS "No revision system found.")
  endif ()

  # Set SeqAn date of last commit.
  if (_SEQAN_WC_LAST_CHANGED_DATE)
    set (SEQAN_DATE "${_SEQAN_WC_LAST_CHANGED_DATE}")
    message (STATUS "  Determined repository date is ${SEQAN_DATE}")
  else ()
    message (STATUS "  Repository date not determined.")
  endif ()

  # Set SeqAn repository revision.
  if (_SEQAN_WC_REVISION)
    set (SEQAN_REVISION "${_SEQAN_WC_REVISION}" CACHE INTERNAL "SeqAn repository revision.")
    message (STATUS "  Determined repository revision is ${SEQAN_REVISION}")
   else ()
    set (SEQAN_REVISION "tarball" CACHE INTERNAL "SeqAn repository revision.")
    message (STATUS "  Repository revision not determined.")
  endif ()
endmacro (seqan_get_repository_info)

# ---------------------------------------------------------------------------
# Macro _seqan_setup_demo_test(cpp_file executable)
#
# When called with the file PATH.cpp, it will check whether PATH.cpp.stdout
# and/or PATH.cpp.stderr exists.  If this is the case then we will add a test
# that runs the demo and compares the standard output/error stream with the
# given file.
#
# Used in seqan_build_demos_develop().
# ---------------------------------------------------------------------------
macro (_seqan_setup_demo_test CPP_FILE EXECUTABLE)
    set (STDOUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${CPP_FILE}.stdout")
    set (STDERR_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${CPP_FILE}.stderr")
    if (EXISTS "${STDOUT_PATH}" OR EXISTS "${STDERR_PATH}")
        # Build the path to the demo_checker.py script.
        set (CHECKER_PATH "${CMAKE_SOURCE_DIR}/util/bin/demo_checker.py")

        # Compose arguments to the demo_checker.py script.
        if (MSVC)
            # Add buildtype path and ".exe" suffix under Windows.
            set (ARGS "--binary-path" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${CMAKE_BUILD_TYPE}/${EXECUTABLE}.exe")
        elseif (WIN32)
          # Add ".exe" suffix for all other Windows compilers, e.g. MinGW.
            set (ARGS "--binary-path" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXECUTABLE}.exe")
        else ()
            set (ARGS "--binary-path" "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${EXECUTABLE}")
        endif ()

        if (EXISTS "${STDOUT_PATH}")
            set (ARGS ${ARGS} "--stdout-path" "${STDOUT_PATH}")
        endif ()
        if (EXISTS "${STDERR_PATH}")
            set (ARGS ${ARGS} "--stderr-path" "${STDERR_PATH}")
        endif()

        # Add the test.
        find_package (PythonInterp)
        if (PYTHONINTERP_FOUND)
          add_test (NAME test_${EXECUTABLE}
                    COMMAND ${PYTHON_EXECUTABLE} ${CHECKER_PATH} ${ARGS})
          #message(STATUS "add_test (NAME test_${EXECUTABLE} COMMAND ${PYTHON_EXECUTABLE} ${CHECKER_PATH} ${ARGS})")
        endif (PYTHONINTERP_FOUND)
    endif ()
endmacro (_seqan_setup_demo_test CPP_FILE)

# ---------------------------------------------------------------------------
# Macro seqan_register_demos([prefix])
#
# Use this in demos directories and subdirectories.
#
# This is only used when doing a Whole SeqAn Release or when developing.
# When doing a SeqAn Release then we copy over the demos, otherwise we build
# them.
# ---------------------------------------------------------------------------

# NOTE that we look with default SeqAn dependencies and also build if some are not found. The demos themselves must contain the appropriate #if preprocessor statements.

# Install all demo source files.
macro (seqan_install_demos_release)
    # Set flags for SeqAn.
    set (SEQAN_FIND_ENABLE_TESTING 0)
    set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DSEQAN_ENABLE_DEBUG=0")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")

    # Get a list of all .cpp and .cu files in the current directory.
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cu)

    # Get path to current source directory, relative from root.  Will be used to install demos in.
    file (RELATIVE_PATH INSTALL_DIR "${SEQAN_ROOT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")

    # Install all demo files into "share/doc/seqan/demos" (demos comes from INSTALL_DIR).
    install (FILES ${ENTRIES} DESTINATION "share/doc/seqan/${INSTALL_DIR}")
endmacro (seqan_install_demos_release)

macro (seqan_build_demos_develop PREFIX)
    # Get a list of all .cpp and .cu files in the current directory.
    file (GLOB_RECURSE ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cu)

    # Find SeqAn with all dependencies.
    set (SEQAN_FIND_DEPENDENCIES ALL)
    find_package (SeqAn REQUIRED)
    find_package (CXX11)
    find_package (OpenMP)
    if (OPENMP_FOUND AND CMAKE_COMPILER_IS_GNUCXX)
        set(SEQAN_LIBRARIES ${SEQAN_LIBRARIES} -lgomp)
    endif()
    # Setup include directories and definitions for SeqAn; flags follow below.
    include_directories (${SEQAN_INCLUDE_DIRS})
    add_definitions (${SEQAN_DEFINITIONS})

    # Supress unused parameter warnings for demos.
    if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-unused-parameter")
    endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

    # Setup flags for CUDA demos.
    seqan_setup_cuda_vars(ARCH sm_20 DEBUG_DEVICE DISABLE_WARNINGS)

    # Add SeqAn flags to CXX and NVCC flags.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${CXX11_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #NOTE(h-2): it is not clear, why the flags have to be added to the defs, but they are not included otherwise
    add_definitions(${CMAKE_CXX_FLAGS})

    # Add all demos with found flags in SeqAn.
    foreach (ENTRY ${ENTRIES})
        string (REPLACE "/" "_" BIN_NAME "${ENTRY}")
        string (REPLACE "\\" "_" BIN_NAME "${BIN_NAME}")
        get_filename_component (BIN_NAME "${BIN_NAME}" NAME_WE)

        get_filename_component (FILE_NAME "${ENTRY}" NAME)
        if ("${FILE_NAME}" MATCHES "\\.cu$")
            if (SEQAN_HAS_CUDA)
                cuda_add_executable(${PREFIX}${BIN_NAME} ${ENTRY})
                target_link_libraries (${PREFIX}${BIN_NAME} ${SEQAN_LIBRARIES})
                if (APPLE AND COMPILER_IS_CLANG)
                    set_target_properties (${PREFIX}${BIN_NAME} PROPERTIES LINK_FLAGS -stdlib=libstdc++)
                endif ()
                _seqan_setup_demo_test (${ENTRY} ${PREFIX}${BIN_NAME})
            endif ()
        else ()
            add_executable(${PREFIX}${BIN_NAME} ${ENTRY})
            target_link_libraries (${PREFIX}${BIN_NAME} ${SEQAN_LIBRARIES})
            _seqan_setup_demo_test (${ENTRY} ${PREFIX}${BIN_NAME})
        endif ()
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_build_demos_develop)

function (seqan_register_demos)
    # Set optional parameter with index 0 into variable PREFIX.
    if (${ARGC} GREATER 0)
        set (PREFIX ${ARGV0})
    else (${ARGC} GREATER 0)
        set (PREFIX "")
    endif (${ARGC} GREATER 0)

    # Install demo source files when releasing and build demos when developing.
    if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE")
        seqan_install_demos_release ()
    elseif ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
        seqan_build_demos_develop ("${PREFIX}")
    endif ()
endfunction (seqan_register_demos)

# ---------------------------------------------------------------------------
# Macro seqan_register_tests ()
# ---------------------------------------------------------------------------

# Switch to testing mode and include all subdirectories with a CMakeLists.txt
# file inside them.  This function should be called in the CMakeLists.txt in
# the tests directories before including subdirectories.
#
# The following will happen:
#
# * Setting definitions SEQAN_ENABLE_DEBUG=1 and SEQAN_ENABLE_TESTING=1.
# * If the ${MODEL} variable is NightlyCoverage OR ExperimentalCoverage,
#   and the compiler is GCC C++ then symbols for test coverate are added.
# * All subdirectories with a CMakeLists.txt file inside will be added.

macro (seqan_register_tests)
    # Setup flags for tests.
    set (SEQAN_FIND_ENABLE_DEBUG TRUE)
    set (SEQAN_FIND_ENABLE_TESTING TRUE)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    # Remove NDEBUG definition for tests.
    string (REGEX REPLACE "-DNDEBUG" ""
            CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
    string (REGEX REPLACE "-DNDEBUG" ""
            CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

    # Conditionally enable coverage mode by setting the appropriate flags.
    if (MODEL STREQUAL "NightlyCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
            add_definitions(-DSEQAN_ENABLE_CHECKPOINTS=0)
        endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif (MODEL STREQUAL "NightlyCoverage")
    if (MODEL STREQUAL "ExperimentalCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
            add_definitions(-DSEQAN_ENABLE_CHECKPOINTS=0)
        endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif (MODEL STREQUAL "ExperimentalCoverage")

    # Add all subdirectories that have a CMakeLists.txt inside them.
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)
    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory(${ENTRY})
            endif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
        endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
    endforeach (ENTRY ${ENTRIES})
endmacro (seqan_register_tests)

