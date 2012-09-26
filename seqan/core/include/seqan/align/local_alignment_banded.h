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
// Interface functions for banded local alignment.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_H_

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
// Function localAlignment()                                    [banded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           int lowerDiag,
                           int upperDiag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);

	clearGaps(row(align, 0));
	clearGaps(row(align, 1));
    clearClipping(row(align, 0));
    clearClipping(row(align, 1));

	LocalAlignmentFinder<TScoreValue> finder;
	_initLocalAlignmentFinder(source(row(align, 0)), source(row(align, 1)), finder, BandedWatermanEggert(), lowerDiag, upperDiag);
	finder.needReinit = false;

	return _localAlignment(finder, row(align, 0), row(align, 1), scoringScheme, 0, lowerDiag, upperDiag, BandedWatermanEggert());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                                     [banded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                           Gaps<TSequenceV, TGapsSpecV> & gapsV,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           int lowerDiag,
                           int upperDiag)
{
	clearGaps(gapsH);
	clearGaps(gapsV);
    clearClipping(gapsH);
    clearClipping(gapsV);

	LocalAlignmentFinder<TScoreValue> finder;
	_initLocalAlignmentFinder(source(gapsH), source(gapsV), finder, BandedWatermanEggert(), lowerDiag, upperDiag);
	finder.needReinit = false;

	return _localAlignment(finder, gapsH, gapsV, scoringScheme, 0, lowerDiag, upperDiag, BandedWatermanEggert());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                      [banded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           int lowerDiag,
                           int upperDiag)
{
    // TODO(holtgrew): This is currently not supported.
    (void)fragmentString;
    (void)strings;
    (void)scoringScheme;
    (void)lowerDiag;
    (void)upperDiag;
    SEQAN_FAIL("Write me!");
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_BANDED_H_
