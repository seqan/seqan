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
# Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
# ============================================================================
# Autotection of the the system name.
# ============================================================================

# Setting CMAKE_SYSTEM_PROCESSOR from the command line does not work, so we
# set it here.  We need this to make the naming of the package file automatic.
if (SEQAN_SYSTEM_PROCESSOR)
  set(CMAKE_SYSTEM_PROCESSOR "${SEQAN_SYSTEM_PROCESSOR}")
endif ()

# some platforms (e.g. FreeBSD) use different names here
if (CMAKE_SYSTEM_PROCESSOR STREQUAL "amd64")
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
  set(CPACK_RPM_PACKAGE_ARCHITECTURE "amd64")
  set(CMAKE_SYSTEM_PROCESSOR "x86_64")
endif ()

if (CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" AND (SEQAN_ARCH_SSE4))
  if (SEQAN_ARCH_AVX2)
    set (CMAKE_SYSTEM_PROCESSOR "x86_64_avx2") 
  else ()
    set (CMAKE_SYSTEM_PROCESSOR "x86_64_sse4") 
  endif ()
endif()

if (NOT DEFINED CPACK_SYSTEM_NAME)
  set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif (NOT DEFINED CPACK_SYSTEM_NAME)

if (${CPACK_SYSTEM_NAME} MATCHES Windows)
  set(CPACK_SYSTEM_NAME Windows-${CMAKE_SYSTEM_PROCESSOR})
elseif (${CPACK_SYSTEM_NAME} MATCHES Darwin)
  set(CPACK_SYSTEM_NAME Mac-${CMAKE_SYSTEM_PROCESSOR})
endif (${CPACK_SYSTEM_NAME} MATCHES Windows)

