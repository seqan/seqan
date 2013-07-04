# - Try to find BZip2
# Once done this will define
#
#  BZIP2_FOUND - system has BZip2
#  BZIP2_INCLUDE_DIR - the BZip2 include directory
#  BZIP2_LIBRARIES - Link these to use BZip2
#  BZIP2_NEED_PREFIX - this is set if the functions are prefixed with BZ2_

#=============================================================================
# Copyright 2006-2009 Kitware, Inc.
# Copyright 2006 Alexander Neundorf <neundorf@kde.org>
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

# This is based on the original FindZLIB.cmake file from the CMake
# distribution.  However, we also look for zlib_d.lib etc. such
# that we can have both debug and release builds.  Furthermore,
# we also look for libbz2.lib since that's the name of the bz2 library
# on windows and the lib is not automatically inferred as on Unices.

FIND_PATH(BZIP2_INCLUDE_DIR bzlib.h)

if (MSVC)
	FIND_LIBRARY(BZIP2_LIBRARY   NAMES libbz2 bz2 bzip2)
	FIND_LIBRARY(BZIP2_LIBRARY_D NAMES libbz2_d bz2_d bzip2_d)
else (MSVC)
    FIND_LIBRARY(BZIP2_LIBRARY   NAMES libbz2 bz2 bzip2 PATH_SUFFIXES "" lib lib32)
    FIND_LIBRARY(BZIP2_LIBRARY_D NAMES libbz2_d bz2_d bzip2_d PATH_SUFFIXES "" lib lib32)
endif (MSVC)

mark_as_advanced(BZIP_LIBRARY BZIP_LIBRARY_D)

# handle the QUIETLY and REQUIRED arguments and set BZip2_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BZip2 DEFAULT_MSG BZIP2_LIBRARY BZIP2_INCLUDE_DIR)

IF (BZIP2_FOUND)
    #message(STATUS "BZip2 libraries found at ${BZIP2_LIBRARY} ${BZIP2_LIBRARY_D}")
    INCLUDE(CheckLibraryExists)
    CHECK_LIBRARY_EXISTS(${BZIP2_LIBRARY} BZ2_bzCompressInit "" BZIP2_NEED_PREFIX)
	IF (MSVC)
		SET(BZIP2_LIBRARIES debug ${BZIP2_LIBRARY_D} optimized ${BZIP2_LIBRARY})
	ELSE (MSVC)
		SET(BZIP2_LIBRARIES ${BZIP2_LIBRARY})
	ENDIF (MSVC)
ELSE (BZIP2_FOUND)
	message(STATUS "BZip2 libraries could not be found!")
ENDIF (BZIP2_FOUND)


MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR BZIP2_LIBRARIES)

