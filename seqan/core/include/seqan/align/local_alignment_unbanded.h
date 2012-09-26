// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Interface functions for unbanded local alignment.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_

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
// Function localAlignment()                                  [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);

    LocalAlignmentFinder<TScoreValue> laFinder;
    return _smithWaterman(row(align, 0), row(align, 1), laFinder, scoringScheme, 0);
}

// ----------------------------------------------------------------------------
// Function localAlignment()                                   [unbanded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                           Gaps<TSequenceV, TGapsSpecV> & gapsV,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    LocalAlignmentFinder<TScoreValue> laFinder;
    return _smithWaterman(gapsH, gapsV, laFinder, scoringScheme, 0);
}

// ----------------------------------------------------------------------------
// Function localAlignment()                     [unbanded, Graph<Alignment<>>]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    return _localAlignment(alignmentGraph,
                           value(stringSet(alignmentGraph), 0),
                           value(stringSet(alignmentGraph), 1),
                           positionToId(stringSet(alignmentGraph), 0),
                           positionToId(stringSet(alignmentGraph), 1),
                           scoringScheme, SmithWaterman());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                    [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    // TODO(holtgrew): This is currently not supported.
    // TODO(holtgrew): Write tests!
    (void)fragmentString;
    (void)strings;
    (void)scoringScheme;
    SEQAN_FAIL("Write me!");
    return 0;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
