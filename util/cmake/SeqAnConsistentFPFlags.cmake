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
# An attempt to produce the same floating point output across all platforms
# for some of our not-so-well written test cases
# ============================================================================

if (COMPILER_LINTEL AND CMAKE_CXX_COMPILER_ID VERSION_LESS 17.0.0)
  # https://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler
  set(SEQAN_CONSISTENT_FP_FLAGS "-fp-model precise -fp-model source -fimf-arch-consistency=true -no-fma")
elseif (COMPILER_WINTEL AND CMAKE_CXX_COMPILER_ID VERSION_LESS 17.0.0)
  # note: bs_tools didn't produce the same output with these settings
  set(SEQAN_CONSISTENT_FP_FLAGS "/fp:precise /fp:source /Qimf-arch-consistency:true /Qfma-")
elseif (COMPILER_LINTEL)
  # at least version 17.0.0
  set(SEQAN_CONSISTENT_FP_FLAGS "-fp-model consistent")
elseif (COMPILER_WINTEL)
  # at least version 17.0.0
  set(SEQAN_CONSISTENT_FP_FLAGS "/fp:consistent")
elseif(COMPILER_GCC)
  # https://randomascii.wordpress.com/2013/07/16/floating-point-determinism/
  # cc1plus: sorry, unimplemented: -fexcess-precision=standard for C++
  set(SEQAN_CONSISTENT_FP_FLAGS "-ffloat-store -ffp-contract=off")
elseif(COMPILER_CLANG)
  set(SEQAN_CONSISTENT_FP_FLAGS "-ffp-contract=off")
elseif(COMPILER_MSVC)
  # https://msdn.microsoft.com/en-us/library/e7s85ffb.aspx
  # note: bs_tools didn't produce the same output with this setting
  set(SEQAN_CONSISTENT_FP_FLAGS "/fp:precise")
endif()
