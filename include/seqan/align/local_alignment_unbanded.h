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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Interface functions for unbanded local alignment.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_

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
// Function localAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn localAlignment
 * @headerfile <seqan/align.h>
 * @headerfile <seqan/align_parallel.h>
 * @brief Computes the best pairwise local alignment using the Smith-Waterman algorithm.
 *
 * @signature TScoreCollection localAlignment([exec,] alignCollection,                  scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreCollection localAlignment([exec,] gapsHCollection, gapsVCollection, scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignment(align,          scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignment(gapsH, gapsV,   scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignment(fragmentString, scoringScheme, [lowerDiag, upperDiag]);
 *
 * @param[in]     exec             The @link ExecutionPolicy @endlink used for the alignment algorithm.
 * @param[in,out] alignCollection  A collection of @link Align @endlink objects containing the sequences to compute the alignment for.
 * @param[in,out] gapsHCollection  A collection of @link Gaps @endlink objects containing the first sequences to compute the alignment for.
 * @param[in,out] gapsVCollection  A collection of @link Gaps @endlink objects containing the second sequences to compute the alignment for.
 * @param[in,out] gapsH Horizontal gapped sequence in alignment matrix. Types: @link Gaps @endlink
 * @param[in,out] gapsV Vertical gapped sequence in alignment matrix. Types: @link Gaps @endlink
 * @param[in,out] align An @link Align @endlink object that stores the alignment. The
 *                      number of rows must be 2 and the sequences must have already
 *                      been set. <tt>align[0]</tt> is the horizontal one in the
 *                      alignment matrix alignment, <tt>align[1]</tt> is the vertical
 *                      one.
 * @param[in,out] fragmentString
 *                      String of @link Fragment @endlink objects. The sequence
 *                      with id <tt>0</tt> is the horizontal one, the sequence
 *                      with id <tt>1</tt> is the vertical one.
 * @param[in] scoringScheme
 *                      The @link Score scoring scheme @endlink to use for the alignment.
 * @param[in] lowerDiag Optional lower diagonal (<tt>int</tt>).
 * @param[in] upperDiag Optional upper diagonal (<tt>int</tt>).
 *
 * @return TScoreVal Score value of the resulting alignment  (Metafunction @link Score#Value @endlink of the type of
 *                   <tt>scoringScheme</tt>).
 *
 * @return TScoreCollection A collection of computed scores for every aligned sequence pair.  The value type of
 *                          this score is @link String @endlink over the score type of the passed scoring scheme.
 *                          (Metafunction: @link Score#Value @endlink of the type of <tt>scoringScheme</tt>).
 *
 * The Waterman-Eggert algorithm (local alignment with declumping) is available through the @link
 * LocalAlignmentEnumerator @endlink class.
 *
 * When using @link Gaps @endlink and @link Align @endlink objects, only parts (i.e. one infix) of each sequence will be
 * aligned.  This will be presented to the user by setting the clipping begin and end position of the gaps (the rows in
 * the case of @link Align @endlink objects).  When using @link Fragment @endlink strings, these parts of the sequences
 * will not appear in any fragment.
 *
 * There exist multiple overloads for this function with two configuration dimensions.
 *
 * First, you can select the type of the target storing the alignment. This can be either an @link Align @endlink
 * object, two @link Gaps @endlink objects, or a string of @link Fragment @endlink objects. @link Align @endlink objects
 * provide an interface to tabular alignments with the restriction of all rows having the same type. Using two @link
 * Gaps @endlink objects has the advantage that you an align sequences with different types, for example @link DnaString
 * @endlink and @link Dna5String @endlink. Using @link Fragment @endlink strings is useful for collecting many pairwise
 * alignments, for example in the construction of @link AlignmentGraph Alignment Graphs @endlink for multiple- sequence
 * alignments (MSA).
 *
 * Second, you can optionally give a band for the alignment using <tt>lowerDiag</tt> and <tt>upperDiag</tt>. The center
 * diagonal has index <tt>0</tt>, the <tt>i</tt>th diagonal below has index <tt>-i</tt>, the <tt>i</tt>th above has
 * index <tt>i</tt>.
 *
 * The examples below show some common use cases.
 *
 * @section Examples
 *
 * Local alignment of two sequences using an @link Align @endlink object.
 *
 * @code{.cpp}
 * Dna5String seqH = "CGATT";
 * Dna5String seqV = "CGAAATT";
 *
 * Align<Dna5String> align;
 * resize(rows(align), 2);
 * assignSource(row(align, 0), seqH);
 * assignSource(row(align, 0), seqV);
 * Score<int, Simple> scoringScheme(2, -1, -2);
 *
 * int result = localAlignment(align, scoringScheme);
 * @endcode
 *
 * Local banded alignment of two sequences using two @link Gaps @endlink objects.
 *
 * @code{.cpp}
 * Dna5String seqH = "CGATT";
 * Gaps<Dna5String, ArrayGaps> gapsH(seqH);
 * DnaString seqV = "CGAAATT";
 * Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);
 *
 * Score<int, Simple> scoringScheme(5, -3, -1, -5);
 *
 * int result = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);
 * @endcode
 *
 * @section Execution policies
 *
 * SeqAn supports parallel and vectorised execution of local pairwise alignments.
 * This additional interface takes not just a single pair of sequences but a collection of sequence pairs
 * either in form of a collection over @link Align @endlink objects or as two collections of @link Gaps @endlink
 * where both collections must have the same size. The collection based interface allows to additionally specify an
 * @link ExecutionPolicy @endlink. This execution policy can be used to select the multi-threaded or vectorised
 * implementation or the combination there of of the alignment algorithm. SeqAn implements an inter-sequence
 * vectorisation scheme which means that <tt>x</tt> alignments are computed in parallel in one
 * <a href="https://en.wikipedia.org/wiki/SIMD">SIMD</a> vector where <tt>x</tt> is the number of elements a vector
 * can compute in parallel. This depends on the architecture's supported SIMD vector width (128 bit, 256 bit or 512 bit)
 * and the selected score type, e.g. <tt>int16_t</tt>. For example on a CPU architecture that supports SSE4 and a score
 * type of <tt>int16_t</tt>, <tt>128/16 = 8</tt> alignments can be computed in parallel on a single core.
 * In addition, the execution policy can be configured for multi-threaded execution, such that chunks of sequence
 * pairs from the initial collection are spawned and executed on different threads. The number of threads can be
 * set via the execution policy.
 *
 * @note In order to get the performance advantage of vectorised execution one has to compile the code with the
 *       respective CPU flags. For most Intel based CPUs the compiler flag <tt>-msse4</tt> can be used for gcc and
 *       clang to enable vectorisation for 128 bit wide registers. For CPUs that support wider register please read
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
 * @include demos/dox/align/local_alignment_unbanded_execution_policy.cpp
 *
 * https://seqan.readthedocs.io/en/develop/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
 *
 * @section References
 *
 * <ul>
 *   <li>Smith TF, Waterman, MS: Identification of Common Molecular Subsequences. J Mol Biol 1981, 147(1):195-7.</li>
 * </ul>
 *
 * @see globalAlignment
 * @see LocalAlignmentEnumerator
 * @see PairwiseLocalAlignmentAlgorithms
 */

// ----------------------------------------------------------------------------
// Function localAlignment()                                  [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, source(row(align, 0)), source(row(align, 1)),
                                            scoringScheme, TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), trace);
    return res;
}

 template <typename TSequence, typename TAlignSpec,
           typename TScoreValue, typename TScoreSpec>
 TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
 {
     SEQAN_ASSERT(length(rows(align)) == 2u);
     if (_usesAffineGaps(scoringScheme, source(row(align, 0)), source(row(align, 1))))
         return localAlignment(align, scoringScheme, AffineGaps());
     else
         return localAlignment(align, scoringScheme, LinearGaps());
 }

// ----------------------------------------------------------------------------
// Function localAlignment()                                   [unbanded, Gaps]
// ----------------------------------------------------------------------------

 template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TScoreValue,
           typename TScoreSpec, typename TTag>
 TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TTag const & tag)
 {
     typedef typename Size<TSequenceH>::Type TSize;
     typedef typename Position<TSequenceH>::Type TPosition;
     typedef TraceSegment_<TPosition, TSize> TTraceSegment;
     typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

     String<TTraceSegment> trace;
     DPScoutState_<Default> dpScoutState;
     TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, source(gapsH), source(gapsV), scoringScheme,
                                             TAlignConfig2(), tag);
     _adaptTraceSegmentsTo(gapsH, gapsV, trace);
     return res;
 }

 template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                           Gaps<TSequenceV, TGapsSpecV> & gapsV,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
     if (_usesAffineGaps(scoringScheme, source(gapsH), source(gapsV)))
         return localAlignment(gapsH, gapsV, scoringScheme, AffineGaps());
     else
         return localAlignment(gapsH, gapsV, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                     [unbanded, Graph<Alignment<>>]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Position<TGraph>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, value(stringSet(alignmentGraph), 0),
                                            value(stringSet(alignmentGraph), 1), scoringScheme, TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), trace);
    return res;
}

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT(length(stringSet(alignmentGraph)) == 2u);

    if (_usesAffineGaps(scoringScheme, stringSet(alignmentGraph)[0], stringSet(alignmentGraph)[1]))
        return localAlignment(alignmentGraph, scoringScheme, AffineGaps());
    else
        return localAlignment(alignmentGraph, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                    [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, value(strings, 0), value(strings, 1), scoringScheme,
                                            TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), trace);
    return res;
}

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT(length(strings) == 2u);

    if (_usesAffineGaps(scoringScheme, strings[0], strings[1]))
        return localAlignment(fragmentString, strings, scoringScheme, AffineGaps());
    else
        return localAlignment(fragmentString, strings, scoringScheme, LinearGaps());
}

// ============================================================================
// Many-vs-Many align interfaces.
// ============================================================================

// ----------------------------------------------------------------------------
// Function localAlignment()             [unbanded, SIMD version, GapsH, GapsV]
// ----------------------------------------------------------------------------

template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSequenceH>>, Is<ContainerConcept<typename Value<TSequenceH>::Type>>>,
                         And<Is<ContainerConcept<TSequenceV>>, Is<ContainerConcept<typename Value<TSequenceV>::Type>>>
                        >, String<TScoreValue>)
localAlignment(TSequenceH & gapSeqSetH,
               TSequenceV & gapSeqSetV,
               Score<TScoreValue, TScoreSpec> const & scoringScheme,
               TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> >    TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type                      TGapModel;

    return _alignWrapper(gapSeqSetH, gapSeqSetV, scoringScheme, TAlignConfig2(), TGapModel());
}

template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec>
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSequenceH>>, Is<ContainerConcept<typename Value<TSequenceH>::Type>>>,
                         And<Is<ContainerConcept<TSequenceV>>, Is<ContainerConcept<typename Value<TSequenceV>::Type>>>
                        >, String<TScoreValue>)
localAlignment(TSequenceH & gapSeqSetH,
               TSequenceV & gapSeqSetV,
               Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
   if (_usesAffineGaps(scoringScheme, source(gapSeqSetH[0]), source(gapSeqSetV[0])))
        return localAlignment(gapSeqSetH, gapSeqSetV, scoringScheme, AffineGaps());
   else
        return localAlignment(gapSeqSetH, gapSeqSetV, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignment()         [unbanded, SIMD version, StringSet<Align>]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
inline String<TScoreValue>
localAlignment(StringSet<Align<TSequence, TAlignSpec> > & alignSet,
               Score<TScoreValue, TScoreSpec> const & scoringScheme,
               TAlgoTag const & algoTag)
{
    typedef Align<TSequence, TAlignSpec>    TAlign;
    typedef typename Row<TAlign>::Type      TGapSequence;

    StringSet<TGapSequence, Dependent<> > gapSetH;
    StringSet<TGapSequence, Dependent<> > gapSetV;
    reserve(gapSetH, length(alignSet));
    reserve(gapSetV, length(alignSet));

    for (auto & align : alignSet)
    {
        appendValue(gapSetH, row(align, 0));
        appendValue(gapSetV, row(align, 1));
    }

    return localAlignment(gapSetH, gapSetV, scoringScheme, algoTag);
}

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
String<TScoreValue> localAlignment(StringSet<Align<TSequence, TAlignSpec> > & align,
                                   Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
   if (_usesAffineGaps(scoringScheme, source(row(align[0], 0)), source(row(align[0], 1))))
        return localAlignment(align, scoringScheme, AffineGaps());
   else
        return localAlignment(align, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignmentScore()                         [unbanded, TSequence]
// ----------------------------------------------------------------------------

/*!
 * @fn localAlignmentScore
 * @headerfile <seqan/align.h>
 * @headerfile <seqan/align_parallel.h>
 * @brief Computes the best global pairwise alignment score.
 *
 * @signature TScoreCollection localAlignmentScore([exec,] seqHCollection, seqVCollection, scoringScheme[, alignConfig][, lowerDiag, upperDiag]);
 * @signature TScoreCollection localAlignmentScore([exec,] seqH, seqVCollection, scoringScheme[, alignConfig][, lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignmentScore(seqH, seqV, scoringScheme[, lowerDiag, upperDiag]);
 *
 * @param[in] exec          The @link ExecutionPolicy @endlink used for the alignment algorithm.
 * @param[in]seqHCollection A collection of sequences aligned pairwise against the respective sequence in
 *                          <tt>seqVCollection</tt>.
 * @param[in]seqVCollection A collection of sequences aligned pairwise against the respective sequence in
 *                          <tt>seqHCollection</tt> or with <tt>seqH</tt>.
 * @param[in] seqH          A single sequence to be aligned against <tt>seqV</tt> or <tt>seqVCollection</tt>.
 * @param[in] seqV          A single sequence to be aligned against <tt>seqH</tt>.
 * @param[in] scoringScheme The scoring scheme to use for the alignment.  Note that the user is responsible for ensuring
 *                          that the scoring scheme is compatible with <tt>algorithmTag</tt>.  Type: @link Score @endlink.
 * @param[in] lowerDiag     Optional lower diagonal.  Types: <tt>int</tt>
 * @param[in] upperDiag     Optional upper diagonal.  Types: <tt>int</tt>
 *
 * @return TScoreVal   Score value of the resulting alignment  (Metafunction: @link Score#Value @endlink of
 *                     the type of <tt>scoringScheme</tt>).
 *
 * @return TScoreCollection A collection of computed scores for every aligned sequence pair.  The value type of
 *                          this score is @link String @endlink over the score type of the passed scoring scheme.
 *                          (Metafunction: @link Score#Value @endlink of the type of <tt>scoringScheme</tt>).
 *
 * This function does not perform the (linear time) traceback step after the (mostly quadratic time) dynamic programming
 * step. Local alignment score can be either used with two sequences or two sets of sequences of equal size.
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
 * @include demos/dox/align/local_alignment_unbanded_score_execution_policy_parallel.cpp
 *
 * The following example shows an example for a wavefront and vectorised execution of global alignments
 * for two collections of large sequences:
 *
 * @include demos/dox/align/local_alignment_unbanded_score_execution_policy_wavefront.cpp
 *
 * @see https://seqan.readthedocs.io/en/develop/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
 * @see localAlignment
 *
 * @datarace thread-safe. No shared state is modified during the execution and concurrent invocations of this function
 * on the same data does not cause any race conditions.
 */

template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
SEQAN_FUNC_DISABLE_IF(And<And<Is<ContainerConcept<TSequenceH>>, Is<ContainerConcept<typename Value<TSequenceH>::Type>>>,
                          And<Is<ContainerConcept<TSequenceV>>, Is<ContainerConcept<typename Value<TSequenceV>::Type>>>
                         >, TScoreValue)
localAlignmentScore(TSequenceH const & seqH,
                    TSequenceV const & seqV,
                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                    TAlgoTag const & /*algoTag*/)
{
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<>, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    DPScoutState_<Default> dpScoutState;
    String<TraceSegment_<unsigned, unsigned> > traceSegments;  // Dummy segments.
    return _setUpAndRunAlignment(traceSegments, dpScoutState, seqH, seqV, scoringScheme, TAlignConfig2{}, TGapModel{});
}

// ----------------------------------------------------------------------------
// Function localAlignmentScore()                   [unbanded, Simd, TSequence]
// ----------------------------------------------------------------------------

template <typename TSeqH,
          typename TSeqV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
inline
SEQAN_FUNC_ENABLE_IF(And<And<Is<ContainerConcept<TSeqH>>, Is<ContainerConcept<typename Value<TSeqH>::Type>>>,
                         And<Is<ContainerConcept<TSeqV>>, Is<ContainerConcept<typename Value<TSeqV>::Type>>>
                        >, String<TScoreValue>)
localAlignmentScore(TSeqH const & stringsH,
                    TSeqV const & stringsV,
                    Score<TScoreValue, TScoreSpec> const & scoringScheme,
                    TAlgoTag const & /*algoTag*/)
{
    SEQAN_ASSERT_EQ(length(stringsH), length(stringsV));
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<>, TracebackOff> TAlignConfig2;
    typedef typename SubstituteAlgoTag_<TAlgoTag>::Type TGapModel;

    return _alignWrapper(stringsH, stringsV, scoringScheme, TAlignConfig2(), TGapModel());
}

// Interface without algorithm tag.
template <typename TSequenceH,
          typename TSequenceV,
          typename TScoreValue, typename TScoreSpec>
auto localAlignmentScore(TSequenceH const & seqH,
                         TSequenceV const & seqV,
                         Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return localAlignmentScore(seqH, seqV, scoringScheme, NeedlemanWunsch());
    else
        return localAlignmentScore(seqH, seqV, scoringScheme, Gotoh());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
