// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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
// Global alignment interface for the unbanded Needleman-Wunsch and Gotoh
// algorithms.
//
// We define the interface functions pretty explicitly (versus just TAlign,
// TFragments etc.) so the candidates the compiler gives when resolution to
// the globalFunction() fails is actually meaningful.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_

namespace seqan2 {

// ============================================================================
// Forwards
// ============================================================================

template <typename TScoreValue, typename TSpec>
class Score;
template <typename TSpec>
class Graph;
template <typename TStringSet, typename TCargo, typename TGraphSpec>
struct Alignment;
template <typename TSize, typename TFragmentSpec>
class Fragment;

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
// Function globalAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn globalAlignment
 * @headerfile <seqan/align.h>
 * @headerfile <seqan/align_parallel.h>
 * @brief Computes the best global pairwise alignment.
 *
 * @signature TScoreCollection globalAlignment([exec,] alignCollection,                  scoringScheme, [alignConfig,] [lowerDiag, upperDiag]);
 * @signature TScoreCollection globalAlignment([exec,] gapsHCollection, gapsVCollection, scoringScheme, [alignConfig,] [lowerDiag, upperDiag]);
 * @signature TScoreVal globalAlignment(align,          scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(gapsH, gapsV,   scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(frags, strings, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(alignGraph,     scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 *
 * @param[in]     exec             The @link ExecutionPolicy @endlink used for the alignment algorithm.
 * @param[in,out] alignCollection  A collection of @link Align @endlink objects containing the sequences to compute the alignment for.
 * @param[in,out] gapsHCollection  A collection of @link Gaps @endlink objects containing the first sequences to compute the alignment for.
 * @param[in,out] gapsVCollection  A collection of @link Gaps @endlink objects containing the second sequences to compute the alignment for.
 * @param[in,out] align        The @link Align @endlink object to use for storing the pairwise alignment.
 * @param[in,out] gapsH        The @link Gaps @endlink object for the first row (horizontal in the DP matrix).
 * @param[in,out] gapsV        The @link Gaps @endlink object for the second row (vertical in the DP matrix).
 * @param[in,out] frags        String of @link Fragment @endlink objects to store alignment in.
 * @param[in]     strings      StringSet of length two with the strings to align.
 * @param[in,out] alignGraph   Alignment Graph for the resulting alignment.  Must be initialized with two strings.
 * @param[in]     scoringScheme The @link Score scoring scheme @endlink to use for the alignment.  Note that
 *                              the user is responsible for ensuring that the scoring scheme is compatible with <tt>algorithmTag</tt>.
 * @param[in]     alignConfig  @link AlignConfig @endlink instance to use for the alignment configuration.
 * @param[in]     lowerDiag    Optional lower diagonal (<tt>int</tt>).
 * @param[in]     upperDiag    Optional upper diagonal (<tt>int</tt>).
 * @param[in]     algorithmTag Tag to select the alignment algorithm (see @link AlignmentAlgorithmTags @endlink).
 *
 * @return TScoreVal   Score value of the resulting alignment  (Metafunction: @link Score#Value @endlink of
 *                     the type of <tt>scoringScheme</tt>).
 *
 * @return TScoreCollection A collection of computed scores for every aligned sequence pair.  The value type of
 *                          this score is @link String @endlink over the score type of the passed scoring scheme.
 *                          (Metafunction: @link Score#Value @endlink of the type of <tt>scoringScheme</tt>).
 *
 * There exist multiple overloads for this function with four configuration dimensions.
 *
 * First, you can select whether begin and end gaps are free in either sequence using <tt>alignConfig</tt>.
 *
 * Second, you can select the type of the target storing the alignment. This can be either an @link Align @endlink
 * object, two @link Gaps @endlink objects, a @link AlignmentGraph @endlink, or a string of @link Fragment @endlink
 * objects. @link Align @endlink objects provide an interface to tabular alignments with the restriction of all rows
 * having the same type. Using two @link Gaps @endlink objects has the advantage that you an align sequences with
 * different types, for example @link DnaString @endlink and @link Dna5String @endlink. @link AlignmentGraph Alignment
 * Graphs @endlink provide a graph-based representation of segment-based colinear alignments. Using @link Fragment
 * @endlink strings is useful for collecting many pairwise alignments, for example in the construction of @link
 * AlignmentGraph Alignment Graphs @endlink for multiple-sequence alignments (MSA).
 *
 * Third, you can optionally give a band for the alignment using <tt>lowerDiag</tt> and <tt>upperDiag</tt>. The center
 * diagonal has index <tt>0</tt>, the <tt>i</tt>th diagonal below has index <tt>-i</tt>, the <tt>i</tt>th above has
 * index <tt>i</tt>.
 *
 * Fourth, you can select the algorithm to use with <tt>algorithmTag</tt>.  This can be one of @link
 * AlignmentAlgorithmTags#NeedlemanWunsch @endlink and @link AlignmentAlgorithmTags#Gotoh @endlink.  The
 * Needleman-Wunsch algorithm supports scoring schemes with linear gap costs only while Gotoh's algorithm also allows
 * affine gap costs.
 *
 * The available alignment algorithms all have some restrictions.  Gotoh's algorithm can handle arbitrary substitution
 * and affine gap scores.  Needleman-Wunsch is limited to linear gap scores.  The implementation of Hirschberg's
 * algorithm is further limited that it does not support <tt>alignConfig</tt> objects or banding.  The implementation of
 * the Myers-Hirschberg algorithm further limits this to only support edit distance (as scores, matches are scored with
 * 0, mismatches are scored with -1).
 *
 * The examples below show some common use cases.
 *
 * @section Examples
 *
 * Global alignment of two sequences using an @link Align @endlink object and
 * the Needleman-Wunsch algorithm.
 *
 * @include demos/dox/align/global_alignment_unbanded.cpp
 *
 * The output is as follows:
 *
 * @include demos/dox/align/global_alignment_unbanded.cpp.stdout
 *
 * Global banded alignment of two sequences using two @link Gaps @endlink objects and the Gotoh algorithm.
 *
 * @include demos/dox/align/global_alignment_banded.cpp
 *
 * The output is as follows:
 *
 * @include demos/dox/align/global_alignment_banded.cpp.stdout
 *
 * https://seqan.readthedocs.io/en/main/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
 *
 * @section Execution policies
 *
 * SeqAn supports parallel and vectorised execution of pairwise alignments.
 * This additional interface takes not just a single pair of sequences but a collection of sequence pairs
 * either in form of a collection over @link Align @endlink objects or as two collections of @link Gaps @endlink
 * where both collections must have the same size. The collection based interface allows to additionally specify an
 * @link ExecutionPolicy @endlink. This execution policy can be used to select the multi-threaded or vectorised
 * implementation or the combination of thereof for the alignment algorithm. SeqAn implements an inter-sequence
 * vectorisation scheme which means that <tt>x</tt> alignments are computed in parallel in one
 * <a href="https://en.wikipedia.org/wiki/SIMD">SIMD</a> vector where <tt>x</tt> is the number of elements a vector
 * can compute in parallel. This depends on the architecture's supported SIMD vector width (128 bit, 256 bit or 512 bit)
 * and the selected score type, e.g. <tt>int16_t</tt>. For example on a CPU architecture that supports SSE4 and a score
 * type of <tt>int16_t</tt>, <tt>128/16 = 8</tt> alignments can be computed in parallel on a single core.
 * In addition, the execution policy can be configured for multi-threaded execution, such that chunks of sequence
 * pairs from the initial collection are spawned and executed on different threads. The number of threads can be
 * set via the execution policy.
 *
 * @note In order to get the performance advantage of vectorised execution, one has to compile the code with the
 *       respective CPU flags. For most Intel based CPUs the compiler flag <tt>-msse4</tt> can be used for gcc and
 *       clang to enable vectorisation for 128 bit wide registers. For CPUs that support wider registers please read
 *       the respective documentation on how to select the correct compilation flags.
 *
 * @warning There are currently some limitations in the use of the execution policy.
 *          The @link WavefrontExecutionPolicy @endlink cannot yet compute the traceback and is only allowed for the
 *          @link globalAlignmentScore @endlink and @link localAlignmentScore @endlink interface.
 *          The banded version only works for collections if all sequences within one collection have the
 *          same size.
 *
 * The following example shows an example for a multi-threaded and vectorised execution of global alignments
 * for two collections of sequences:
 *
 * @include demos/dox/align/global_alignment_unbanded_execution_policy.cpp
 *
 * @section References
 *
 * <ul>
 *   <li>Needleman SB, Wunsch CD: A general method applicable to the search for similarities in the amino acid sequence
 *       of two proteins. J Mol Biol 1970, 48(3): 443-53.</li>
 *   <li>Gotoh O: An improved algorithm for matching biological sequences. J Mol Biol 1982, 162(3):705-8</li>
 * </ul>
 *
 * @see localAlignment
 * @see globalAlignmentScore
 * @see AlignmentAlgorithmTags
 */

// ----------------------------------------------------------------------------
// Function globalAlignment()                                 [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec, typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            TAlgoTag const & /*algoTag*/)
{
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    String<TTraceSegment> trace;
    TScoreValue res;
    DPScoutState_<Default> dpScoutState;
    res  = _setUpAndRunAlignment(trace, dpScoutState, source(row(align, 0)), source(row(align, 1)), scoringScheme,
                                 TAlignConfig2(), TGapModel());

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), trace);
    return res;
}

// Interface without AlignConfig<>.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(align, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                                  [unbanded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            TAlgoTag const & /*algoTag*/)
{
    typedef typename Size<TSequenceH>::Type TSize;
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    String<TTraceSegment> traceSegments;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(traceSegments, dpScoutState, source(gapsH), source(gapsV), scoringScheme,
                                            TAlignConfig2(), TGapModel());
    _adaptTraceSegmentsTo(gapsH, gapsV, traceSegments);
    return res;
}

// Interface without AlignConfig<>.
template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                   [unbanded, Graph<Alignment<> >]
// ----------------------------------------------------------------------------

// Full interface.
template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            TAlgoTag const & /*algoTag*/)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Position<TGraph>::Type TPosition;
    typedef typename Size<TGraph>::Type TSize;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    String<TTraceSegment> traceSegments;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(traceSegments, dpScoutState, value(stringSet(alignmentGraph), 0),
                                            value(stringSet(alignmentGraph), 1), scoringScheme, TAlignConfig2(),
                                            TGapModel());

    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.
template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                   [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.
template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                            TAlgoTag const & /*algoTag*/)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    String<TTraceSegment> traceSegments;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(traceSegments, dpScoutState, value(strings, 0), value(strings, 1),
                                            scoringScheme, TAlignConfig2(), TGapModel());

    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.
template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()
// ----------------------------------------------------------------------------

/*!
 * @fn globalAlignmentScore
 * @headerfile <seqan/align.h>
 * @headerfile <seqan/align_parallel.h>
 * @brief Computes the best global pairwise alignment score.
 *
 * @signature TScoreCollection globalAlignmentScore([exec,] seqHCollection, seqVCollection, scoringScheme[, alignConfig][, lowerDiag, upperDiag]);
 * @signature TScoreCollection globalAlignmentScore([exec,] seqH, seqVCollection, scoringScheme[, alignConfig][, lowerDiag, upperDiag]);
 * @signature TScoreVal globalAlignmentScore(strings,    scoringScheme[, alignConfig][, lowerDiag, upperDiag][, algorithmTag]);
 * @signature TScoreVal globalAlignmentScore(seqH, seqV, {MyersBitVector | MyersHirschberg});
 * @signature TScoreVal globalAlignmentScore(strings,    {MyersBitVector | MyersHirschberg});
 *
 * @param[in] exec           The @link ExecutionPolicy @endlink used for the alignment algorithm.
 * @param[in]seqHCollection  A collection of sequences aligned pairwise against the respective sequence in
 *                           <tt>seqVCollection</tt>.
 * @param[in]seqVCollection  A collection of sequences aligned pairwise against the respective sequence in
 *                           <tt>seqHCollection</tt> or with <tt>seqH</tt>.
 * @param[in] strings        A @link StringSet @endlink containing two sequences.  Type: StringSet.
 * @param[in] seqH           A single sequence to be aligned against <tt>seqV</tt> or <tt>seqVCollection</tt>.
 * @param[in] seqV           A single sequence to be aligned against <tt>seqH</tt>.
 * @param[in] alignConfig    The @link AlignConfig @endlink to use for the alignment.  Type: AlignConfig
 * @param[in] scoringScheme  The scoring scheme to use for the alignment.  Note that the user is responsible for ensuring
 *                           that the scoring scheme is compatible with <tt>algorithmTag</tt>.
 *                           Type: @link Score @endlink.
 * @param[in] lowerDiag      Optional lower diagonal.  Types: <tt>int</tt>
 * @param[in] upperDiag      Optional upper diagonal.  Types: <tt>int</tt>
 * @param[in] algorithmTag   The Tag for picking the alignment algorithm. Types: @link PairwiseLocalAlignmentAlgorithms
 *                           @endlink.
 *
 * @return TScoreVal   Score value of the resulting alignment  (Metafunction: @link Score#Value @endlink of
 *                     the type of <tt>scoringScheme</tt>).
 *
 * @return TScoreCollection A collection of computed scores for every aligned sequence pair.  The value type of
 *                          this score is @link String @endlink over the score type of the passed scoring scheme.
 *                          (Metafunction: @link Score#Value @endlink of the type of <tt>scoringScheme</tt>).
 *
 * This function does not perform the (linear time) traceback step after the (mostly quadratic time) dynamic programming
 * step.  Note that Myers' bit-vector algorithm does not compute an alignment (only in the Myers-Hirschberg variant) but
 * scores can be computed using <tt>globalAlignmentScore</tt>.
 * Global alignment score can be either used with two sequences or two sets of sequences of equal size.
 *
 * The same limitations to algorithms as in @link globalAlignment @endlink apply.  Furthermore, the
 * <tt>MyersBitVector</tt> and <tt>MyersHirschberg</tt> variants can only be used without any other parameter.
 *
 * @section Execution policies
 *
 * SeqAn supports parallel and vectorised execution of pairwise alignments.
 * This additional interface takes not just a single pair of sequences but a collection of sequences where both
 * collections must have the same size.
 * The collection based interface allows to additionally specify an
 * @link ExecutionPolicy @endlink. This execution policy can be used to select the multi-threaded or vectorised
 * implementation or the combination thereof for the alignment algorithm. SeqAn implements an inter-sequence
 * vectorisation scheme which means that <tt>x</tt> alignments are computed in parallel in one
 * <a href="https://en.wikipedia.org/wiki/SIMD">SIMD</a> vector where <tt>x</tt> is the number of elements a vector
 * can compute in parallel. This depends on the architecture's supported SIMD vector width (128 bit, 256 bit or 512 bit)
 * and the selected score type, e.g. <tt>int16_t</tt>. For example on a CPU architecture that supports SSE4 and a score
 * type of <tt>int16_t</tt>, <tt>128/16 = 8</tt> alignments can be computed in parallel on a single core.
 *
 * In addition, the execution policy can be configured for multi-threaded execution, such that either chunks of sequence
 * pairs from the initial collection are spawned and executed on different threads or an intra-sequence parallelization
 * is used to parallelize a single alignment. In total the following execution modes are possible:
 * <i>sequential</i>, <i>parallel</i>, <i>wave-front</i>, <i>vectorized</i>, <i>parallel+vectorized</i> and
 * <i>wave-front+vectorized</i>.
 *
 * The wave-front execution can be selected via the @link WavefrontExecutionPolicy @endlink, which can also be combined
 * with a vectorized execution. In addition the wave-front execution parallelizes a single pairwise alignment, while the
 * standard @link ParallelismTags#Parallel @endlink specialization does only parallelizes the sequence set via chunking.
 * Note,
 *
 * @note In order to get the performance advantage of vectorised execution one has to compile the code with the
 *       respective CPU flags. For most Intel based CPUs the compiler flag <tt>-msse4</tt> can be used for gcc and
 *       clang to enable vectorisation for 128 bit wide registers. For CPUs that support wider register please read
 *       the respective documentation on how to select the correct compilation flags.
 *
 * @warning There are currently some limitations in the use of the execution policy.
 *          The banded version is at the moment only supported for the following execution modes: <i>sequential</i>,
 *          <i>parallel</i>, <i>vectorized</i> and <i>parallel+vectorized</i>.
 *          The banded version only works for collections if all sequences within one collection have the
 *          same size.
 *
 * The following example shows an example for a multi-threaded and vectorised execution of global alignments
 * for two collections of sequences:
 *
 * @include demos/dox/align/global_alignment_unbanded_score_execution_policy_parallel.cpp
 *
 * The following example shows an example for a wavefront and vectorised execution of global alignments
 * for two collections of large sequences:
 *
 * @include demos/dox/align/global_alignment_unbanded_score_execution_policy_wavefront.cpp
 *
 * @see https://seqan.readthedocs.io/en/main/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
 * @see globalAlignment
 *
 * @datarace thread-safe. No shared state is modified during the execution and concurrent invocations of this function
 * on the same data does not cause any race conditions.
 */

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [unbanded, 2 Strings]
// ----------------------------------------------------------------------------

template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
SEQAN_FUNC_DISABLE_IF(And<And<Is<ContainerConcept<TSequenceH>>, Is<ContainerConcept<typename Value<TSequenceH>::Type>>>,
                          And<Is<ContainerConcept<TSequenceV>>, Is<ContainerConcept<typename Value<TSequenceV>::Type>>>
                         >,
                      TScoreValue)
globalAlignmentScore(TSequenceH const & seqH,
                     TSequenceV const & seqV,
                     Score<TScoreValue, TScoreSpec> const & scoringScheme,
                     AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                     TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    DPScoutState_<Default> dpScoutState;
    String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments.
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, TAlignConfig2(), TGapModel());
}

// Interface without AlignConfig<>.
template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
auto globalAlignmentScore(TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
auto globalAlignmentScore(TSequenceH const & seqH,
                          TSequenceV const & seqV,
                          Score<TScoreValue, TScoreSpec> const & scoringScheme,
                          AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec>
auto globalAlignmentScore(TSequenceH const & seqH,
                                 TSequenceV const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [unbanded, StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                                 TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    SEQAN_ASSERT_EQ(length(strings), 2u);

    DPScoutState_<Default> dpScoutState;
    String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments.
    return _setUpAndRunAlignment(traceSegments, dpScoutState, strings[0], strings[1], scoringScheme, TAlignConfig2(),
                                 TGapModel());
}

// Interface without AlignConfig<>.
template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 TAlgoTag const & algoTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()       [unbanded, SIMD version, 2x StringSet]
// ----------------------------------------------------------------------------

template <typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSeqH>>, Is<ContainerConcept<typename Value<TSeqH>::Type>>>,
                         And<Is<ContainerConcept<TSeqV>>, Is<ContainerConcept<typename Value<TSeqV>::Type>>>
                        >,
                    String<TScoreValue>)
globalAlignmentScore(TSeqH const & stringsH,
                     TSeqV const & stringsV,
                     Score<TScoreValue, TScoreSpec> const & scoringScheme,
                     AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                     TAlgoTag const & /*algoTag*/)
{
    SEQAN_ASSERT_EQ(length(stringsH), length(stringsV));
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    return _alignWrapper(stringsH, stringsV, scoringScheme, TAlignConfig2(), TGapModel());
}

// --------------------------------------------------------------------------------
// Function globalAlignmentScore()    [unbanded, SIMD version, String vs StringSet]
// --------------------------------------------------------------------------------
template <typename TString1, typename TString2, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
String<TScoreValue> globalAlignmentScore(TString1 const & stringH,
                                         StringSet<TString2, TSpec> const & stringsV,
                                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                         AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                                         TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    return _alignWrapper(stringH, stringsV, scoringScheme, TAlignConfig2(), TGapModel());
}

// Interface without AlignConfig<>.
template <typename TString1, typename TString2, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
String<TScoreValue> globalAlignmentScore(TString1 const & stringH,
                                         StringSet<TString2, TSpec> const & stringsV,
                                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                         TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(stringH, stringsV, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TString1, typename TString2, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
String<TScoreValue> globalAlignmentScore(TString1 const & stringH,
                                         StringSet<TString2, TSpec> const & stringsV,
                                         Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                         AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(stringH, stringsV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignmentScore(stringH, stringsV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TString1, typename TString2, typename TSpec,
          typename TScoreValue, typename TScoreSpec>
String<TScoreValue> globalAlignmentScore(TString1 const & stringH,
                                         StringSet<TString2, TSpec> const & stringsV,
                                         Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(stringH, stringsV, scoringScheme, alignConfig);
}

// ============================================================================
// Many-vs-Many align interfaces.
// ============================================================================

// ----------------------------------------------------------------------------
// Function globalAlignment()            [unbanded, SIMD version, GapsH, GapsV]
// ----------------------------------------------------------------------------

template <typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSeqH>>, Is<ContainerConcept<typename Value<TSeqH>::Type>>>,
                         And<Is<ContainerConcept<TSeqV>>, Is<ContainerConcept<typename Value<TSeqV>::Type>>>
                        >,
                      String<TScoreValue>)
globalAlignment(TSeqH & gapSeqSetH,
                TSeqV & gapSeqSetV,
                Score<TScoreValue, TScoreSpec> const & scoringScheme,
                AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & /*alignConfig*/,
                TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec>              TAlignConfig;
    typedef typename SubstituteAlignConfig_<TAlignConfig>::Type         TFreeEndGaps;
    typedef AlignConfig2<DPGlobal, DPBandConfig<BandOff>, TFreeEndGaps> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type                 TGapModel;

    return _alignWrapper(gapSeqSetH, gapSeqSetV, scoringScheme, TAlignConfig2(), TGapModel());
}

// Interface without algorithm tag.
template <typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSeqH>>, Is<ContainerConcept<typename Value<TSeqH>::Type>>>,
                         And<Is<ContainerConcept<TSeqV>>, Is<ContainerConcept<typename Value<TSeqV>::Type>>>
                        >,
                     String<TScoreValue>)
globalAlignment(TSeqH & gapSeqSetH,
                TSeqV & gapSeqSetV,
                Score<TScoreValue, TScoreSpec> const & scoringScheme,
                AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(gapSeqSetH, gapSeqSetV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(gapSeqSetH, gapSeqSetV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSeqH>>, Is<ContainerConcept<typename Value<TSeqH>::Type>>>,
                         And<Is<ContainerConcept<TSeqV>>, Is<ContainerConcept<typename Value<TSeqV>::Type>>>
                        >,
                     String<TScoreValue>)
globalAlignment(TSeqH & gapSeqSetH,
                TSeqV & gapSeqSetV,
                Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapSeqSetH, gapSeqSetV, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()        [unbanded, SIMD version, StringSet<Align>]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec, typename TAlgoTag>
String<TScoreValue> globalAlignment(StringSet<Align<TSequence, TAlignSpec> > & alignSet,
                                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                    AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                    TAlgoTag const & algoTag)
{
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Row<TAlign>::Type TGapSequence;

    StringSet<TGapSequence, Dependent<> > gapSetH;
    StringSet<TGapSequence, Dependent<> > gapSetV;
    reserve(gapSetH, length(alignSet));
    reserve(gapSetV, length(alignSet));

    for (auto & align : alignSet)
    {
        appendValue(gapSetH, row(align, 0));
        appendValue(gapSetV, row(align, 1));
    }

    return globalAlignment(gapSetH, gapSetV, scoringScheme, alignConfig, algoTag);
}

// Interface without AlignConfig<>.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
String<TScoreValue> globalAlignment(StringSet<Align<TSequence, TAlignSpec> > & align,
                                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                    TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
String<TScoreValue> globalAlignment(StringSet<Align<TSequence, TAlignSpec> > & align,
                                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                    AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(align, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.
template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
String<TScoreValue> globalAlignment(StringSet<Align<TSequence, TAlignSpec> > & align,
                                    Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig);
}

}  // namespace seqan2

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_
