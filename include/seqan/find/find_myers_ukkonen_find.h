// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_FIND_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_FIND_H

#include <seqan/find/find_myers_ukkonen_base.h>
#include <seqan/find/find_myers_ukkonen_pattern.h>

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// find

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern)
{
    typedef typename Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >::TPatternState TPatternState;
    return find(finder, pattern, static_cast<TPatternState&>(pattern));
}

template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
                  int const minScore)
{
    setScoreLimit(state, minScore);
    return find(finder, pattern, state);
}

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  TNeedle const & needle,
                  PatternState_<TNeedle2, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
                  int const minScore)
{
    setScoreLimit(state, minScore);
    return find(finder, needle, state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
                  int const minScore)
{
    return find(finder, pattern, pattern, minScore); //static cast
}

} // namespace seqan

#endif // #ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_FIND_H
