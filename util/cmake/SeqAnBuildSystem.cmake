# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
# SEQAN_RELEASE_APPS
# SEQAN_RELEASE_LIBRARY
# APP:${app_name}

include (SeqAnUsabilityAnalyzer)
include (CheckCXXCompilerFlag)

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
# Is it a 32 bit platform?
# ---------------------------------------------------------------------------
if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(SEQAN_32BIT_TARGET_PLATFORM 1)
    set(SEQAN_64BIT_TARGET_PLATFORM 0)
else()
    set(SEQAN_32BIT_TARGET_PLATFORM 0)
    set(SEQAN_64BIT_TARGET_PLATFORM 1)
endif()

# ---------------------------------------------------------------------------

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
# Macro seqan_register_apps ()
#
# Register all apps by adding their subdirectories if they are to be built
# (SEQAN_RELEASE and APP:${app} modes).
# ---------------------------------------------------------------------------

macro (seqan_register_apps)

    # enable static linkage for seqan apps
    if (SEQAN_STATIC_APPS AND (NOT CMAKE_SYSTEM_NAME MATCHES "Windows"))
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        if (APPLE)
            # static build not supported on apple, but at least we can include gcc libs
            if (CMAKE_COMPILER_IS_GNUCXX)
                set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++")
            endif (CMAKE_COMPILER_IS_GNUCXX)
        else (APPLE)
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

            # make sure -rdynamic isn't added automatically
            set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
            # make sure -fPIC isn't added automatically
            set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)

            # for ENTIRELY UNKOWN reasons find_package returns libs for static
            # linking (.a) when building as SEQAN_RELEASE_APPS or APP:*  but
            # .so when building DEVELOP. In the latter case it also encloses the
            # static libs with -Bdynamic which turns static off for system libs.
            # Here we remove these (so statics works), but only for NON-DEVELOP
            if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
                # make sure -Wl,-Bdynamic isn't added automatically
                set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
            endif (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
        endif (APPLE)
    endif (SEQAN_STATIC_APPS AND (NOT CMAKE_SYSTEM_NAME MATCHES "Windows"))

    # Get all direct entries of the current source directory into ENTRIES.
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)

    # Add all values from ${ENTRIES} that are subdirectories and have a file
    # CMakeListst.txt.
    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP") OR
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

    # GENERAL SETUP
    set (_CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set (CMAKE_INCLUDE_PATH ${_CMAKE_INCLUDE_PATH} CACHE STRING "")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1")
#     set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DSEQAN_ENABLE_DEBUG=1" PARENT_SCOPE)
    # Enable global exception handler for all "official" stuff
    set (SEQAN_DEFINITIONS "${SEQAN_DEFINITIONS} -DSEQAN_GLOBAL_EXCEPTION_HANDLER=1")
    set (CMAKE_RUNTIME_OUTPUT_DIRECTORY
         ${PROJECT_BINARY_DIR}/bin)

    # Set Warnings
    if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        # Disable warnings about unsecure (although standard) functions.
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /W2 /D_SCL_SECURE_NO_WARNINGS") # TODO(h-2): raise this to W4
    else ()
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -W -Wall -pedantic -fstrict-aliasing -Wstrict-aliasing")
        set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64)
        if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            # disable some warnings on ICC
            set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -wd3373,2102")
        endif ()
    endif ()

    if (NOT SEQAN_BUILD_SYSTEM)
        set (SEQAN_BUILD_SYSTEM "DEVELOP" CACHE STRING "Build/Release mode to select. One of DEVELOP SEQAN_RELEASE, APP:\${APP_NAME}. Defaults to DEVELOP.")
    endif (NOT SEQAN_BUILD_SYSTEM)
    set (SEQAN_APP_VERSION "0.0.0" CACHE STRING "Version of the application.")
    set (SEQAN_NIGHTLY_RELEASE FALSE CACHE BOOL "Set to TRUE to enable nightly app releases.")

    ## options
    if (NOT SEQAN_64BIT_TARGET_PLATFORM)
        set (SEQAN_ARCH_SSE4 FALSE)
    endif ()

    if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set (SEQAN_STATIC_APPS FALSE)
        set (SEQAN_ARCH_SSE4 FALSE)
        set (SEQAN_ARCH_NATIVE FALSE)
    endif ()

    if (${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD")
        set (SEQAN_STATIC_APPS FALSE)
        if (SEQAN_ARCH_NATIVE)
            set (SEQAN_ARCH_NATIVE FALSE)
            set (SEQAN_ARCH_SSE4 TRUE)
            message (STATUS "OpenBSD does not support native, but SSE4 was activated instead.")
        endif ()
    endif ()

    if (SEQAN_STATIC_APPS)
        message (STATUS "Building static apps.")
        # implementation in seqan_register_apps()
    endif ()

    # machine specific optimizations
    if (SEQAN_ARCH_NATIVE)
        message (STATUS "Building binaries optimized for this specific CPU. They might not work elsewhere.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=native")
        if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ipo -no-prec-div -fp-model fast=2 -xHOST")
        endif ()
    elseif (SEQAN_ARCH_SSE4)
        message (STATUS "Building optimized binaries up to SSE4 and POPCNT.")
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mmmx -msse -msse2 -msse3 -mssse3 -msse4")
        if (NOT ${CMAKE_SYSTEM_NAME} STREQUAL "OpenBSD")
            set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -mpopcnt")
        endif ()
        if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
            set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -ipo -no-prec-div -fp-model fast=2")
        endif ()
    endif ()
    # TODO(h-2): for icc on windows, replace the " -" in SEQAN_CXX_FLAGS with " /"
    #            find out whether clang/c2 takes - or / options

    # automatic c++ standard detection/selection
    if (NOT MSVC)
        # find the highest c++ standard and select it
        check_cxx_compiler_flag("-std=c++11" CXX11_DETECTED)
        check_cxx_compiler_flag("-std=c++14" CXX14_DETECTED)
    endif ()

    set (CXX11_FOUND ${CXX11_DETECTED} CACHE INTERNAL "Availability of c++11")
    set (CXX14_FOUND ${CXX14_DETECTED} CACHE INTERNAL "Availability of c++14")

    if (CXX14_FOUND)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    elseif (CXX11_FOUND)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif ()
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
    if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY"))

        # Install SeqAn LICENSE, README.rst, CHANGELOG.rst files.
        install (FILES LICENSE
                       README.rst
                       CHANGELOG.rst
                 DESTINATION share/doc/seqan)
        # Install pkg-config file, except on Windows.
        if (NOT CMAKE_SYSTEM_NAME MATCHES Windows)
            configure_file("util/pkgconfig/seqan.pc.in" "${CMAKE_BINARY_DIR}/util/pkgconfig/seqan-${SEQAN_VERSION_MAJOR}.pc" @ONLY)
            install(FILES "${CMAKE_BINARY_DIR}/util/pkgconfig/seqan-${SEQAN_VERSION_MAJOR}.pc" DESTINATION share/pkgconfig)
        endif (NOT CMAKE_SYSTEM_NAME MATCHES Windows)
        # Install FindSeqAn TODO(h-2) rename FindSeqAn.cmake to FindSeqAn${SEQAN_VERSION_MAJOR}.cmake after 2.x cycle
        install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/util/cmake/FindSeqAn.cmake" DESTINATION share/cmake/Modules)

        # Install headers
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
  elseif (CMAKE_VERSION VERSION_LESS "3.1") # TXZ support since 3.1
    set (CPACK_GENERATOR "ZIP;TBZ2")
  else ()
    set (CPACK_GENERATOR "ZIP;TXZ")
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
# Function seqan_get_version()
#
# Sets the variables SEQAN_VERSION, SEQAN_VERSION_MAJOR, SEQAN_VERSION_MINOR,
# SEQAN_VERSION_PATCH, determined from seqan/version.h
# ---------------------------------------------------------------------------

macro (seqan_get_version)
  # Read from CMAKE_SOURCE_DIR the /include/seqan/version.h
  get_filename_component(_SEQAN_VERSION_H "${CMAKE_SOURCE_DIR}/include/seqan/version.h" ABSOLUTE)
  # If file wasn't found seqan version is set to 0.0.0
  set (_SEQAN_VERSION_IDS MAJOR MINOR PATCH PRE_RELEASE)
  foreach (_ID ${_SEQAN_VERSION_IDS})
    set(_SEQAN_VERSION_${_ID} "0")
  endforeach(_ID ${_SEQAN_VERSION_IDS})

  # Error log if version.h not found, otherwise read version from
  # version.h and cache it.
  if (NOT EXISTS ${_SEQAN_VERSION_H})
    message ("")
    message ("ERROR: Could not determine SeqAn version.")
    message ("Could not find file: ${_SEQAN_VERSION_H}")
  else ()
    foreach (_ID ${_SEQAN_VERSION_IDS})
      file (STRINGS ${_SEQAN_VERSION_H} _VERSION_${_ID} REGEX ".*SEQAN_VERSION_${_ID}.*")
      string (REGEX REPLACE ".*SEQAN_VERSION_${_ID}[ |\t]+([0-9a-zA-Z]+).*" "\\1" SEQAN_VERSION_${_ID} ${_VERSION_${_ID}})
    endforeach(_ID ${_SEQAN_VERSION_IDS})
  endif ()
  set (SEQAN_VERSION_STRING "${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}")
endmacro (seqan_get_version)

# ---------------------------------------------------------------------------
# Function seqan_get_repository_info()
#
# Sets the variables SEQAN_DATE and SEQAN_REVISION determined from git.
# ---------------------------------------------------------------------------

macro (seqan_get_repository_info)
  set (_SEQAN_GIT_DIR "${CMAKE_SOURCE_DIR}/.git")
  message (STATUS "  Selected repository dir: ${CMAKE_SOURCE_DIR}")
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
    # icc doesn't cope with spaces..
    string(REPLACE " " "_" SEQAN_DATE "${SEQAN_DATE}")
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

function (seqan_register_demos PREFIX)
    # Get a list of all .cpp and .cu files in the current directory.
    file (GLOB_RECURSE ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*.cu)

    # Find SeqAn with all dependencies.
    find_package (OpenMP)
    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn REQUIRED)

    # Supress unused parameter warnings for demos.
    if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
        set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -Wno-unused-parameter")
    endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)

    # Add SeqAn flags to CXX and NVCC flags.
    # Set to PARENT_SCOPE since this macro is executed from within a function which declares it's own scope.
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}" PARENT_SCOPE)
    # Setup include directories and definitions for SeqAn; flags follow below.
    include_directories (${SEQAN_INCLUDE_DIRS})
    add_definitions (${SEQAN_DEFINITIONS})

    # Add all demos with found flags in SeqAn.
    foreach (ENTRY ${ENTRIES})
        set (SKIP FALSE)
        # workaround a bug in llvm35 on FreeBSD
        if ((ENTRY MATCHES "zip") AND
            (${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD") AND
            (COMPILER_IS_CLANG) AND
            (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5.0) AND
            (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6.0))
            set (SKIP TRUE)
        # bug in visual studio
        elseif ((ENTRY MATCHES "queue_example.cpp") AND
                COMPILER_IS_MSVC)
            set (SKIP TRUE)
        endif ()

        if (SKIP)
            message(STATUS "${ENTRY} skipped on this platform." )
        else (SKIP)
            string (REPLACE "/" "_" BIN_NAME "${ENTRY}")
            string (REPLACE "\\" "_" BIN_NAME "${BIN_NAME}")
            get_filename_component (BIN_NAME "${BIN_NAME}" NAME_WE)

            get_filename_component (FILE_NAME "${ENTRY}" NAME)
            add_executable(${PREFIX}${BIN_NAME} ${ENTRY})
            target_link_libraries (${PREFIX}${BIN_NAME} ${SEQAN_LIBRARIES})
            _seqan_setup_demo_test (${ENTRY} ${PREFIX}${BIN_NAME})
        endif (SKIP)
    endforeach (ENTRY ${ENTRIES})
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
    set (SEQAN_DEFINITIONS ${SEQAN_DEFINITIONS} -DSEQAN_ENABLE_TESTING=1)

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
        endif (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
    endif (MODEL STREQUAL "NightlyCoverage")
    if (MODEL STREQUAL "ExperimentalCoverage")
        if (CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
            set (LDFLAGS "${LDFLAGS} -fprofile-arcs -ftest-coverage")
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

