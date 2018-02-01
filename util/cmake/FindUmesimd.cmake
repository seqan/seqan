# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
# Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
# Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
# ============================================================================
#
#.rst:
# FindUmesimd
# -----------
#
# Try to find UME::SIMD includes
#
# Once done this will define
#
# ::
#
#   UMESIMD_FOUND - system has UME::SIMD
#   UMESIMD_INCLUDE_DIR - the UME::SIMD include directory
#   UMESIMD_VERSION_STRING - the version of UME::SIMD found

find_path(UMESIMD_CHILD_INCLUDE_DIR UMESimd.h PATH_SUFFIXES umesimd)

if (UMESIMD_CHILD_INCLUDE_DIR AND EXISTS "${UMESIMD_CHILD_INCLUDE_DIR}/UMESimd.h")
    get_filename_component(UMESIMD_INCLUDE_DIR ${UMESIMD_CHILD_INCLUDE_DIR} DIRECTORY)

    file(STRINGS "${UMESIMD_CHILD_INCLUDE_DIR}/UMESimd.h" UMESIMD_H REGEX "#define UME_SIMD_VERSION_(MAJOR|MINOR|PATCH)")
    string(REGEX REPLACE "#define UME_SIMD_VERSION_(MAJOR|MINOR|PATCH) " "" UMESIMD_VERSION_STRING "${UMESIMD_H}")
    string(REGEX REPLACE ";" "." UMESIMD_VERSION_STRING "${UMESIMD_VERSION_STRING}")
endif ()

# set UMESIMD_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(
    Umesimd
    REQUIRED_VARS UMESIMD_INCLUDE_DIR
    VERSION_VAR UMESIMD_VERSION_STRING)

mark_as_advanced(UMESIMD_INCLUDE_DIR)
