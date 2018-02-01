// ==========================================================================
//                           pair_align_global.cpp
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: René Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Source file for pairwise global alignment.
// ==========================================================================

#include "pair_align_lib.h"

// Set default values if preprocessor symbols are not set.
#ifndef SUFFIX_GAP_TOP
#define SUFFIX_GAP_TOP 0
#endif  // SUFFIX_GAP_TOP

#ifndef SUFFIX_GAP_LEFT
#define SUFFIX_GAP_LEFT 0
#endif  // SUFFIX_GAP_LEFT

#ifndef SUFFIX_GAP_RIGHT
#define SUFFIX_GAP_RIGHT 0
#endif  // SUFFIX_GAP_RIGHT

#ifndef SUFFIX_GAP_BOTTOM
#define SUFFIX_GAP_BOTTOM 0
#endif  // SUFFIX_GAP_BOTTOM

// Defines the suffix for every generated object file.
#define PAIR_ALIGN_GLOBAL_FUNC_SUFFIX SEQAN_JOIN(SEQAN_JOIN(SUFFIX_GAP_TOP, SUFFIX_GAP_LEFT),   \
SEQAN_JOIN(SUFFIX_GAP_RIGHT, SUFFIX_GAP_BOTTOM))

// This macro expands to the function resulting from the set preprocessor values SUFFIX_GAP_*.
// Note that the cmake script generates for every AlignConfig combination
// a new object file where the preprocessor definitions are set accordingly.
#define PAIR_ALIGN_DEFINE_GLOBAL_FUNC(name)                                                                       \
void name(Options const & options)                                                                                 \
{                                                                                                                  \
    AlignConfig<SUFFIX_GAP_TOP, SUFFIX_GAP_LEFT, SUFFIX_GAP_RIGHT, SUFFIX_GAP_BOTTOM> config;                      \
    if (options.method == "nw")                                                                                    \
        pairAlignConfig(options, NeedlemanWunsch(), config);                                                       \
    else                                                                                                           \
        pairAlignConfig(options, Gotoh(), config);                                                                 \
}

PAIR_ALIGN_DEFINE_GLOBAL_FUNC(SEQAN_JOIN(pairAlignGlobal_, PAIR_ALIGN_GLOBAL_FUNC_SUFFIX))
