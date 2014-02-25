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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Thread-safe / lock-free sequence operations.
// ==========================================================================

#ifndef SEQAN_PARALLEL_SEQUENCE_H_
#define SEQAN_PARALLEL_SEQUENCE_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _incLength()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TParallel>
inline typename Size<String<TValue, Alloc<TSpec> > >::Type
_incLength(String<TValue, Alloc<TSpec> > & me, Tag<TParallel> const & tag)
{
    return atomicInc(me.data_end, tag) - begin(me, Standard());
}

// ----------------------------------------------------------------------------
// Function appendValue(Insist, Parallel); Atomic
// ----------------------------------------------------------------------------

template <typename TTargetValue, typename TTargetSpec, typename TValue>
inline void
appendValue(String<TTargetValue, TTargetSpec> & me, TValue const & _value, Insist, Parallel)
{
    valueConstruct(begin(me, Standard()) + _incLength(me, Parallel()) - 1, _value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_SEQUENCE_H_
