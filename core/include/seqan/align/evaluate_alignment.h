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
// Compute alignment score given a pairwise alignment.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_
#define CORE_INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AlignmentStats
// ----------------------------------------------------------------------------

/*!
 * @class AlignmentStats
 * @headerfile <seqan/align.h>
 * @brief Statistics about a tabular alignment.
 *
 * The default constructor initializes all members to 0.
 *
 * @var unsigned AlignmentStats::numGapOpens;
 * @brief Number of gap open events.
 *
 * @var unsigned AlignmentStats::numGapExtensions;
 * @brief Number of gap extension events.
 *
 * @var unsigned AlignmentStats::numMatches;
 * @brief Number of match (identity) events.
 *
 * @var unsigned AlignmentStats::numMismatches;
 * @brief Number of mismatch (not identity) events.
 *
 * @var unsigned AlignmentStats::numPositiveScores;
 * @brief Number of residues aligned with positive score (0 is counted as positive).
 *
 * @var unsigned AlignmentStats::numNegativeScores;
 * @brief Number of residues aligned with negative score.
 *
 * @var int AlignmentStats::alignmentScore;
 * @brief The resulting alignment score.
 */

struct AlignmentStats
{
    // Number of gap opens/gap extensions.
    unsigned numGapOpens;
    unsigned numGapExtensions;
    // Number of matches, mismatches.
    unsigned numMatches;
    unsigned numMismatches;
    // Number of aligned residues with positive/negative scores.
    unsigned numPositiveScores;
    unsigned numNegativeScores;

    // The alignment score.
    int alignmentScore;

    AlignmentStats() : numGapOpens(0), numGapExtensions(0), numMatches(0), numMismatches(0),
                       numPositiveScores(0), numNegativeScores(0), alignmentScore(0)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn AlignmentStats#clear
 * @brief Clear AlignmentStats object.
 *
 * @signature void clear(stats);
 *
 * @param[in,out] stats AlignmentStats object to clear.
 */

inline
void clear(AlignmentStats & stats)
{
    stats.numGapOpens = 0;
    stats.numGapExtensions = 0;
    stats.numMatches = 0;
    stats.numMismatches = 0;
    stats.numPositiveScores = 0;
    stats.numNegativeScores = 0;
    stats.alignmentScore = 0;
}

// ----------------------------------------------------------------------------
// Function computeAlignmentStats()
// ----------------------------------------------------------------------------

/*!
 * @fn computeAlignmentStats
 * @headerfile <seqan/align.h>
 * @brief Compute alignment statistics.
 *
 * @signature TScoreVal computeAlignmentStats([stats, ]align, scoringScheme);
 *
 * @param[out] stats The @link AlignmentStats @endlink object to store alignment statistics in.
 * @param[in]  align The @link Align @endlink object to score.
 * @param[in]  score The @link Score @endlink object to use for the scoring scheme.
 *
 * @see AlignmentStats
 *
 * @section Examples
 *
 * @include demos/align/compute_alignment_stats.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/compute_alignment_stats.cpp.stdout
 */

template <typename TSource, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
TScoreVal computeAlignmentStats(AlignmentStats & stats,
                                Align<TSource, TAlignSpec> const & align,
                                Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ_MSG(length(rows(align)), 2u, "Only works with pairwise alignments.");
    SEQAN_ASSERT_EQ_MSG(length(row(align, 0)), length(row(align, 1)), "Invalid alignment!");
    clear(stats);

    typedef Align<TSource, TAlignSpec> const TAlign;
    typedef typename Row<TAlign>::Type TGaps;
    typedef typename Iterator<TGaps, Standard>::Type TGapsIter;
    typedef typename Value<typename Source<TGaps>::Type>::Type TAlphabet;

    // Get iterators.
    TGapsIter it0 = begin(row(align, 0));
    TGapsIter itEnd0 = end(row(align, 0));
    TGapsIter it1 = begin(row(align, 1));
    TGapsIter itEnd1 = end(row(align, 1));

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;

    for (; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1)
    {
        if (isGap(it0))
        {
            if (isGapOpen0)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            isGapOpen0 = true;
        }
        else
        {
            isGapOpen0 = false;
        }

        if (isGap(it1))
        {
            if (!isGapOpen1)
            {
                stats.numGapOpens += 1;
                stats.alignmentScore += scoreGapOpen(scoringScheme);
            }
            else
            {
                stats.numGapExtensions += 1;
                stats.alignmentScore += scoreGapExtend(scoringScheme);
            }
            isGapOpen1 = true;
        }
        else
        {
            isGapOpen1 = false;
        }

        if (!isGap(it0) && !isGap(it1))
        {
            // Compute the alignment score and register in stats.
            TAlphabet c0 = *it0, c1 = *it1;
            TScoreVal scoreVal = score(scoringScheme, c0, c1);
            stats.alignmentScore += scoreVal;
            // Register other statistics.
            bool isMatch = (c0 == c1);
            bool isPositive = (scoreVal >= 0);
            stats.numMatches += isMatch;
            stats.numMismatches += !isMatch;
            stats.numPositiveScores += isPositive;
            stats.numNegativeScores += !isPositive;
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);

    return stats.alignmentScore;
}

template <typename TGaps, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
TScoreVal computeAlignmentStats(Align<TGaps, TAlignSpec> const & align,
                                Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    AlignmentStats stats;
    (void)stats;
    return computeAlignmentStats(stats, align, scoringScheme);
}

// ----------------------------------------------------------------------------
// Function alignmentIdentityString()
// ----------------------------------------------------------------------------

/*!
 * @fn Align#alignmentIdentityString()
 * @headerfile <seqan/align.h>
 * @brief Returns a CharString of the similarity of positions in an alignment
 *
 * @signature CharString alignmentIdentityString();
 *
 * @param align The @link Align @endlink object to score for similarity.
 * @param match The character to use for similarity at match positions
 * @param mismatch The character to use for similarity at mismatch positions
 * @param gap The character to use for similarity at gap positions
 *
 * @signature CharString A string of characters representing alignment position similarity
 *
 * @section Examples
 *
 * @include demos/align/compute_match_string.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/compute_match_string.cpp.stdout
 */


template <typename TGaps>
String<char> alignmentIdentityString(TGaps & queryRow,
                                     TGaps & referenceRow,
                                     char const matchChar,
                                     char const mismatchChar,
                                     char const gapChar)
{
    // Define short-hand names for various types to be used
    typedef typename Iterator<TGaps, Standard>::Type TGapsIter;
    typedef typename Value<typename Source<TGaps>::Type>::Type TAlphabet;

    // Define and size the char string we will return
    String<char> matchString;
    unsigned alignLength = std::min(length(queryRow), length(referenceRow));
    resize(matchString, alignLength, Exact());

    // Get iterators.
    TGapsIter it0 = begin(queryRow);
    TGapsIter itEnd0 = end(queryRow);
    TGapsIter it1 = begin(referenceRow);
    TGapsIter itEnd1 = end(referenceRow);

    // Counter for the position in the similarity string
    unsigned pos = 0;

    for (; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1, ++pos)
    {
        if (isGap(it0) || isGap(it1))
        {
            matchString[pos] = gapChar;
        }
        else
        {
            // Compute the alignment score and register in stats.
            TAlphabet c0 = *it0, c1 = *it1;
            if (c0 == c1)
            {
                matchString[pos] = matchChar;
            }
            else
            {
                matchString[pos] = mismatchChar;
            }
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);

    return matchString;
}

template <typename TGaps>
String<char> alignmentIdentityString(TGaps const & queryRow,
                                     TGaps const & referenceRow)
{
    return alignmentIdentityString(queryRow, referenceRow, '|', '*', ' ');
}

template <typename TAlphabet, typename TAlignSpec>
String<char> alignmentIdentityString(Align<TAlphabet, TAlignSpec> const & align,
                                     char const matchChar,
                                     char const mismatchChar,
                                     char const gapChar)
{
    // Check that the supplied alignment is pairwise
    SEQAN_ASSERT_EQ_MSG(length(rows(align)), 2u, "Only works with pairwise alignments.");

    typedef Align<TAlphabet, TAlignSpec> const TAlign;
    typedef typename Row<TAlign>::Type TGaps;

    TGaps queryRow = row(align, 0);
    TGaps referenceRow = row(align, 1);

    return alignmentIdentityString(queryRow, referenceRow, matchChar, mismatchChar, gapChar);
}

template <typename TSource, typename TAlignSpec>
String<char> alignmentIdentityString(Align<TSource, TAlignSpec> const & align)
{
    return alignmentIdentityString(align, '|', '*', ' ');
}


// ----------------------------------------------------------------------------
// Function alignmentSimilarityString()
// ----------------------------------------------------------------------------

/*!
 * @fn Align#alignmentSimilarityString()
 * @headerfile <seqan/align.h>
 * @brief Returns a CharString of the similarity of positions in an alignment
 *
 * @signature CharString computeAlignmentSimilarityString();
 *
 * @param align The @link Align @endlink object to score for similarity.
 * @param score The @link Score @endlink object to use for the scoring scheme.
 * @param match The character to use for similarity at match positions
 * @param positive The character to use for similarity at positive-scoring non-match positions
 * @param mismatch The character to use for similarity at negative-scoring non-match positions
 * @param gap The character to use for similarity at gap positions
 *
 * @signature CharString A string of characters representing alignment position similarity
 *
 * @section Examples
 *
 * @include demos/align/compute_match_string.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/compute_match_string.cpp.stdout
 */

template <typename TGaps, typename TScoreVal, typename TScoreSpec>
String<char> alignmentSimilarityString(TGaps & queryRow,
                                       TGaps & referenceRow,
                                       Score<TScoreVal, TScoreSpec> const & scoringScheme,
                                       char const matchChar,
                                       char const positiveChar,
                                       char const mismatchChar,
                                       char const gapChar)
{
    // Define short-hand names for various types to be used
    typedef typename Iterator<TGaps, Standard>::Type TGapsIter;
    typedef typename Value<typename Source<TGaps>::Type>::Type TAlphabet;

    // Define and size the char string we will return
    String<char> matchString;
    unsigned alignLength = std::min(length(queryRow), length(referenceRow));
    resize(matchString, alignLength, Exact());

    // Get iterators.
    TGapsIter it0 = begin(queryRow);
    TGapsIter itEnd0 = it0 + alignLength;
    TGapsIter it1 = begin(referenceRow);
    TGapsIter itEnd1 = it1 + alignLength;

    // Counter for the position in the similarity string
    unsigned pos = 0;

    for (; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1, ++pos)
    {
        if (isGap(it0) || isGap(it1))
        {
            matchString[pos] = gapChar;
        }
        else
        {
            TAlphabet c0 = *it0, c1 = *it1;
            if (c0 == c1)
            {
                matchString[pos] = matchChar;
            }
            else if (score(scoringScheme, c0, c1) > 0)
            {
                matchString[pos] = positiveChar;
            }
            else 
            {
                matchString[pos] = mismatchChar;
            }
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);

    return matchString;
}


template <typename TGaps, typename TScoreVal, typename TScoreSpec>
String<char> alignmentSimilarityString(TGaps const & queryRow,
                                       TGaps const & referenceRow,
                                       Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    return alignmentSimilarityString(queryRow, referenceRow, scoringScheme, 
                                     '|', ':', '.', ' ');
}

template <typename TAlphabet, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
String<char> alignmentSimilarityString(Align<TAlphabet, TAlignSpec> const & align,
                                       Score<TScoreVal, TScoreSpec> const & scoringScheme,
                                       char const matchChar,
                                       char const positiveChar,
                                       char const mismatchChar,
                                       char const gapChar)
{
    // Check that the supplied alignment is pairwise
    SEQAN_ASSERT_EQ_MSG(length(rows(align)), 2u, "Only works with pairwise alignments.");

    typedef Align<TAlphabet, TAlignSpec> const TAlign;
    typedef typename Row<TAlign>::Type TGaps;

    TGaps queryRow = row(align, 0);
    TGaps referenceRow = row(align, 1);

    return alignmentSimilarityString(queryRow, referenceRow, scoringScheme, 
                                     matchChar, positiveChar, mismatchChar, gapChar);
}

template <typename TSource, typename TAlignSpec, typename TScoreVal, typename TScoreSpec>
String<char> alignmentSimilarityString(Align<TSource, TAlignSpec> const & align,
                                       Score<TScoreVal, TScoreSpec> const & scoringScheme)
{
    return alignmentSimilarityString(align, scoringScheme, '|', ':', '.', ' ');
}


}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_ALIGN_EVALUATE_ALIGNMENT_H_
