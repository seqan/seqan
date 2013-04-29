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
# CMake file that is included from the root CMakeLists.txt.
#
# It sets variables for configuring CPack and then invokes CPack such that the
# "make package" command is available.
# ============================================================================

if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_APPS"))
    include (InstallRequiredSystemLibraries)
endif ()

if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_APPS"))
    include (SetCPackSystemName)

    # NOTE that you have to run "make docs" before running cpack.  The reason
    # is that we cannot add dependencies to the install target at the moment.
    # See: http://public.kitware.com/Bug/view.php?id=8438

    # ===========================================================================
    # Archive Packages (.tar & .tar.bz2)
    # ===========================================================================

    if (WIN32)
        SET (CPACK_GENERATOR "ZIP")
    else ()
        SET (CPACK_GENERATOR "ZIP;TBZ2")
    endif ()
    if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE")
      SET(CPACK_PACKAGE_NAME "seqan")
    elseif ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY")
      SET(CPACK_PACKAGE_NAME "seqan-library")
    elseif ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_APPS")
      SET(CPACK_PACKAGE_NAME "seqan-apps")
    endif ()
    SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "SeqAn - The C++ library for sequence analysis.")
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>")
    SET(CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
    SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
    SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")

    seqan_get_version()

    if (SEQAN_NIGHTLY_RELEASE)
      include (GetCurrentDate)
      set (CPACK_PACKAGE_VERSION "${CURRENT_YEAR}${CURRENT_MONTH}${CURRENT_DAY}")
      set (CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION}")
    else ()
      set (CPACK_PACKAGE_VERSION "${SEQAN_VERSION_MAJOR}.${SEQAN_VERSION_MINOR}.${SEQAN_VERSION_PATCH}")
      set (CPACK_PACKAGE_VERSION "${SEQAN_VERSION}")
    endif (SEQAN_NIGHTLY_RELEASE)
    SET(CPACK_PACKAGE_VERSION_MAJOR "${SEQAN_VERSION_MAJOR}")
    SET(CPACK_PACKAGE_VERSION_MINOR "${SEQAN_VERSION_MINOR}")
    SET(CPACK_PACKAGE_VERSION_PATCH "${SEQAN_VERSION_PATCH}")
    SET(CPACK_PACKAGE_INSTALL_DIRECTORY "SeqAn ${CPACK_PACKAGE_VERSION}")

    if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY")
        set (CPACK_PACKAGE_FILE_NAME "seqan-library-${CPACK_PACKAGE_VERSION}")
    endif ("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE_LIBRARY")

    # Should be the last include.
    INCLUDE(CPack)
endif ()
