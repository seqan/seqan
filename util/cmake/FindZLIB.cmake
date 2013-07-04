# - Find zlib
# Find the native ZLIB includes and library.
# Once done this will define
#
#  ZLIB_INCLUDE_DIRS   - where to find zlib.h, etc.
#  ZLIB_LIBRARIES      - List of libraries when using zlib.
#  ZLIB_FOUND          - True if zlib found.
#
#  ZLIB_VERSION_STRING - The version of zlib found (x.y.z)
#  ZLIB_VERSION_MAJOR  - The major version of zlib
#  ZLIB_VERSION_MINOR  - The minor version of zlib
#  ZLIB_VERSION_PATCH  - The patch version of zlib
#  ZLIB_VERSION_TWEAK  - The tweak version of zlib
#
# The following variable are provided for backward compatibility
#
#  ZLIB_MAJOR_VERSION  - The major version of zlib
#  ZLIB_MINOR_VERSION  - The minor version of zlib
#  ZLIB_PATCH_VERSION  - The patch version of zlib

#=============================================================================
# Copyright 2001-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# This is based on the original FinZLIB.cmake file from the CMake
# distribution.  However, we also look for zlib_d.lib etc. such
# that we can have both debug and release builds.

FIND_PATH(ZLIB_INCLUDE_DIR zlib.h
    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Zlib;InstallPath]/include"
)

SET(ZLIB_NAMES z zlib zlib1 zlibd zlibd1)
if (MSVC)
	SET(ZLIB_NAMES_D z_d zlib_d zlib1_d zlibd_d zlibd1_d)
endif (MSVC)

FIND_LIBRARY(ZLIB_LIBRARY
    NAMES
        ${ZLIB_NAMES}
    PATHS
        "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Zlib;InstallPath]/lib"
    PATH_SUFFIXES
        "" lib lib32
)

if (MSVC)
	FIND_LIBRARY(ZLIB_LIBRARY_D
		NAMES
			${ZLIB_NAMES_D}
		PATHS
			"[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Zlib;InstallPath]/lib"
	)
endif (MSVC)

MARK_AS_ADVANCED(ZLIB_LIBRARY ZLIB_INCLUDE_DIR)

IF(ZLIB_INCLUDE_DIR AND EXISTS "${ZLIB_INCLUDE_DIR}/zlib.h")
    FILE(STRINGS "${ZLIB_INCLUDE_DIR}/zlib.h" ZLIB_H REGEX "^#define ZLIB_VERSION \"[^\"]*\"$")

    STRING(REGEX REPLACE "^.*ZLIB_VERSION \"([0-9]+).*$" "\\1" ZLIB_VERSION_MAJOR "${ZLIB_H}")
    STRING(REGEX REPLACE "^.*ZLIB_VERSION \"[0-9]+\\.([0-9]+).*$" "\\1" ZLIB_VERSION_MINOR  "${ZLIB_H}")
    STRING(REGEX REPLACE "^.*ZLIB_VERSION \"[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1" ZLIB_VERSION_PATCH "${ZLIB_H}")
    SET(ZLIB_VERSION_STRING "${ZLIB_VERSION_MAJOR}.${ZLIB_VERSION_MINOR}.${ZLIB_VERSION_PATCH}")

    # only append a TWEAK version if it exists:
    SET(ZLIB_VERSION_TWEAK "")
    IF( "${ZLIB_H}" MATCHES "^.*ZLIB_VERSION \"[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*$")
        SET(ZLIB_VERSION_TWEAK "${CMAKE_MATCH_1}")
        SET(ZLIB_VERSION_STRING "${ZLIB_VERSION_STRING}.${ZLIB_VERSION_TWEAK}")
    ENDIF( "${ZLIB_H}" MATCHES "^.*ZLIB_VERSION \"[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*$")

    SET(ZLIB_MAJOR_VERSION "${ZLIB_VERSION_MAJOR}")
    SET(ZLIB_MINOR_VERSION "${ZLIB_VERSION_MINOR}")
    SET(ZLIB_PATCH_VERSION "${ZLIB_VERSION_PATCH}")
ENDIF()

# handle the QUIETLY and REQUIRED arguments and set ZLIB_FOUND to TRUE if 
# all listed variables are TRUE
if (NOT CMAKE_CURRENT_LIST_DIR)  # CMAKE_CURRENT_LIST_DIR only from cmake 2.8.3.
    get_filename_component(CMAKE_CURRENT_LIST_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
endif (NOT CMAKE_CURRENT_LIST_DIR)
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ZLIB REQUIRED_VARS ZLIB_LIBRARY ZLIB_INCLUDE_DIR
                                       VERSION_VAR ZLIB_VERSION_STRING)

IF(ZLIB_FOUND)
	#message(STATUS "Zlib found at ${ZLIB_LIBRARY} ${ZLIB_LIBRARY_D}")
    SET(ZLIB_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR})
	if (MSVC)
		SET(ZLIB_LIBRARIES debug ${ZLIB_LIBRARY_D} optimized ${ZLIB_LIBRARY})
	else (MSVC)
		SET(ZLIB_LIBRARIES ${ZLIB_LIBRARY})
	endif (MSVC)
ELSE()
	message(STATUS "Zlib could not be found!")
ENDIF()

