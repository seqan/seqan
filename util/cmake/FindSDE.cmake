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
# FindSDE
# ---------
#
# Try to find Intel® Software Development Emulator
#
# Once done this will define
#
# ::
#
#   SDE_FOUND - system has Intel® SDE
#   SDE_EXECUTABLE - the Intel® SDE executable (full path)
#   SDE_VERSION_STRING - the version of Intel® SDE found

set(SDE_EXECUTABLE)
set(SDE_VERSION_STRING)

# first try 64bit and then 32bit version of sde
foreach(sde_exec "sde64" "sde")
    find_program(SDE_EXECUTABLE
        NAMES "${sde_exec}"
        PATHS
            /usr/bin
            /usr/local/bin
            /opt/local/bin
        DOC "Intel Software Development Emulator"
    )

    if (SDE_EXECUTABLE)
      execute_process(COMMAND "${SDE_EXECUTABLE}" "--version" OUTPUT_VARIABLE SDE_OUTPUT_VERSION)
      if (SDE_OUTPUT_VERSION MATCHES "Software Development Emulator")
          string(REGEX REPLACE ".*Software Development Emulator.*Version:[ ]*([0-9]+\\.[^ ]+).*" "\\1" SDE_VERSION_STRING "${SDE_OUTPUT_VERSION}")
          break()
      endif()
    endif()
endforeach()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
    SDE
    FOUND_VAR SDE_FOUND
    REQUIRED_VARS SDE_EXECUTABLE SDE_VERSION_STRING
    VERSION_VAR SDE_VERSION_STRING
    FAIL_MESSAGE "Could NOT find SDE: Intel(R) Software Development Emulator (sde)"
)

mark_as_advanced(SDE_EXECUTABLE)
