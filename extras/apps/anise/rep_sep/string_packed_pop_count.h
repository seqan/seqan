// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_STRING_PACKED_POP_COUNT_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_STRING_PACKED_POP_COUNT_H_

#include <seqan/sequence.h>
#include <seqan/misc/misc_bit_twiddling.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function popCount()
// ----------------------------------------------------------------------------

template <typename THostSpec>
inline size_t popCount(String<bool, Packed<THostSpec> > const & str)
{
    typedef String<bool, Packed<THostSpec> > TBitVector;
    typedef typename Host<TBitVector>::Type THost;
    typedef typename Iterator<THost const, Standard>::Type TConstIterator;
    typedef typename Value<THost>::Type THostValue;
    
    if (empty(host(str)))
        return 0;
        
    TConstIterator itBegin = begin(host(str), seqan::Standard());
    TConstIterator itEnd = end(host(str), seqan::Standard());

    unsigned lastBits = length(host(str)) % BitsPerValue<THostValue>::VALUE;
    size_t result = 0;
    if (lastBits)
        --itEnd;
                           
    for (; itBegin != itEnd; ++itBegin)
        result += popCount(itBegin->i);

    if (lastBits)
    {
        __uint64 mask = (1 << (lastBits - 1)) - 1 | 1;
        result += popCount(itEnd->i & mask);
    }

    return result;
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_STRING_PACKED_POP_COUNT_H_
