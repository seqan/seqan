// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_
#define INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_

namespace seqan
{

#if SEQAN_ALIGN_SIMD_PROFILE
struct AlignSimdProfile_
{
    double preprTimer = 0.0;
    double alignTimer = 0.0;
    double traceTimer = 0.0;

    void clear()
    {
        preprTimer = 0.0;
        alignTimer = 0.0;
        traceTimer = 0.0;
    }
};

    double timer = 0.0;

AlignSimdProfile_ profile;
#endif

// ============================================================================
// Forwards
// ============================================================================

template <unsigned LENGTH>
struct VectorLength_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSimdVector_, typename TSeqH_, typename TSeqV_>
struct SimdAlignVariableLengthTraits
{
    using TSimdVector   = TSimdVector_;
    using TSeqH         = TSeqH_;
    using TSeqV         = TSeqV_;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _createSimdRepImpl()
// ----------------------------------------------------------------------------

template <typename TSimdVecs,
          typename TStrings,
          size_t ...I>
inline void
_createSimdRepImpl(TSimdVecs & vecs,
                   TStrings const & strs,
                   std::index_sequence<I...> const & /*unsued*/)
{
    for (size_t pos = 0; pos < length(vecs); ++pos)
        fillVector(vecs[pos], strs[I][pos]...);
}

template <typename TSimdVecs,
          typename TStrings>
inline void
_createSimdRepImpl(TSimdVecs & simdStr,
                   TStrings const & strings)
{
    using TSimdVec = typename Value<TSimdVecs>::Type;
    constexpr auto length = LENGTH<TSimdVec>::VALUE;
    _createSimdRepImpl(simdStr, strings, std::make_index_sequence<length>());
}

// Actually precompute value if scoring scheme is score matrix and simd version.
template <typename TSeqValue,
typename TScoreValue, typename TScore>
inline SEQAN_FUNC_ENABLE_IF(And<Is<SimdVectorConcept<TSeqValue> >, IsScoreMatrix_<TScore> >, TSeqValue)
_precomputeScoreMatrixOffset(TSeqValue const & seqVal,
                             Score<TScoreValue, ScoreSimdWrapper<TScore> > const & /*score*/)
{
    return createVector<TSeqValue>(TScore::VALUE_SIZE) * seqVal;
}

// ----------------------------------------------------------------------------
// Function _prepareAndRunSimdAlignment()
// ----------------------------------------------------------------------------

template <typename TStringSimdH,
          typename TStringSimdV,
          typename TSequencesH,
          typename TSequencesV>
inline void
_prepareSimdAlignment(TStringSimdH & stringSimdH,
                      TStringSimdV & stringSimdV,
                      TSequencesH const & seqH,
                      TSequencesV const & seqV,
                      DPScoutState_<SimdAlignEqualLength> const & /*unused*/)
{
    resize(stringSimdH, length(seqH[0]));
    resize(stringSimdV, length(seqV[0]));
    _createSimdRepImpl(stringSimdH, seqH);
    _createSimdRepImpl(stringSimdV, seqV);
}

template <typename TStringSimdH,
          typename TStringSimdV,
          typename TSequencesH,
          typename TSequencesV,
          typename TTraits>
inline void
_prepareSimdAlignment(TStringSimdH & stringSimdH,
                      TStringSimdV & stringSimdV,
                      TSequencesH const & seqH,
                      TSequencesV const & seqV,
                      DPScoutState_<SimdAlignVariableLength<TTraits> > & state)
{
    using TSimdVecH SEQAN_TYPEDEF_FOR_DEBUG = typename Value<TStringSimdH>::Type;
    SEQAN_ASSERT_EQ(length(seqH), length(seqV));
    SEQAN_ASSERT_EQ(static_cast<decltype(length(seqH))>(LENGTH<TSimdVecH>::VALUE), length(seqH));

    using TSimd = typename Value<TStringSimdH>::Type;
    using TPadStringH = ModifiedString<typename Value<TSequencesH const>::Type, ModPadding>;
    using TPadStringV = ModifiedString<typename Value<TSequencesV const>::Type, ModPadding>;

    resize(state.lengthsH, length(seqH));
    resize(state.lengthsV, length(seqV));
    resize(state.endsH, length(seqH));
    resize(state.endsV, length(seqV));

    for (unsigned i = 0; i < length(seqH); ++i)
    {
        state.lengthsH[i] = length(seqH[i]);
        state.lengthsV[i] = length(seqV[i]);
        state.endsH[i] = i;
        state.endsV[i] = i;
    }

    setHost(state.sortedEndsH, state.lengthsH);
    setHost(state.sortedEndsV, state.lengthsV);
    setCargo(state.sortedEndsH, state.endsH);
    setCargo(state.sortedEndsV, state.endsV);

    auto maxLengthLambda = [](auto& lengthLhs, auto& lengthRhs) { return lengthLhs < lengthRhs; };
    sort(state.sortedEndsH, maxLengthLambda, Serial());
    sort(state.sortedEndsV, maxLengthLambda, Serial());

    size_t maxH = back(state.sortedEndsH);
    size_t maxV = back(state.sortedEndsV);

    // and we have to prepare the bit masks of the DPScoutState
    resize(state.masks,  maxV + 1, createVector<TSimd>(0));
    resize(state.masksV, maxV + 1, createVector<TSimd>(0));
    resize(state.masksH, maxH + 1, createVector<TSimd>(0));

    // Create Stringset with padded strings.
    StringSet<TPadStringH> paddedH;
    StringSet<TPadStringV> paddedV;
    resize(paddedH, length(seqH));
    resize(paddedV, length(seqV));

    for(unsigned i = 0; i < length(seqH); ++i)
    {
        setHost(paddedH[i], seqH[i]);
        setHost(paddedV[i], seqV[i]);
        expand(paddedH[i], maxH);
        expand(paddedV[i], maxV);

        // mark the original end position of the alignment in the masks (with -1, all bits set)
        assignValue(state.masksH[state.lengthsH[i]], i, -1);
        assignValue(state.masksV[state.lengthsV[i]], i, -1);
    }

    // now create SIMD representation
    resize(stringSimdH, maxH);
    resize(stringSimdV, maxV);
    _createSimdRepImpl(stringSimdH, paddedH);
    _createSimdRepImpl(stringSimdV, paddedV);
}

template <typename TResult,
          typename TTraces,
          typename TSequencesH,
          typename TSequencesV,
          typename TScore,
          typename TAlgo, typename TBand, typename TFreeEndGaps, typename TTraceback,
          typename TGapModel>
inline void
_prepareAndRunSimdAlignment(TResult & results,
                            TTraces & traces,
                            TSequencesH const & seqH,
                            TSequencesV const & seqV,
                            TScore const & scoringScheme,
                            AlignConfig2<TAlgo, TBand, TFreeEndGaps, TTraceback> const & alignConfig,
                            TGapModel const & /*gapModel*/)
{
    auto seqLengthH = length(seqH[0]);
    auto seqLengthV = length(seqV[0]);
    auto zipView = makeZipView(seqH, seqV);
    bool allSameLength = std::all_of(begin(zipView, Standard()), end(zipView, Standard()),
                                     [seqLengthH, seqLengthV](auto param)
                                     {
                                         return (length(std::get<0>(param)) == seqLengthH) &&
                                                (length(std::get<1>(param)) == seqLengthV);
                                     });

    String<TResult, Alloc<OverAligned> > stringSimdH;
    String<TResult, Alloc<OverAligned> > stringSimdV;
    if(allSameLength)
    {
        DPScoutState_<SimdAlignEqualLength> state;
        _prepareSimdAlignment(stringSimdH, stringSimdV, seqH, seqV, state);
        results = _setUpAndRunAlignment(traces, state, stringSimdH, stringSimdV, scoringScheme, alignConfig, TGapModel());
    }
    else
    {
        DPScoutState_<SimdAlignVariableLength<SimdAlignVariableLengthTraits<TResult, TSequencesH, TSequencesV> > > state;
        _prepareSimdAlignment(stringSimdH, stringSimdV, seqH, seqV, state);

        state.dimV = length(stringSimdV) + 1;
        state.isLocalAlignment = IsLocalAlignment_<TAlgo>::VALUE;
        state.right = IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE;
        state.bottom = IsFreeEndGap_<TFreeEndGaps, DPLastRow>::VALUE;

        results = _setUpAndRunAlignment(traces, state, stringSimdH, stringSimdV, scoringScheme, alignConfig, TGapModel());
    }
}

// ----------------------------------------------------------------------------
// Function _alignWrapperSimd(); Score; StringSet vs. StringSet
// ----------------------------------------------------------------------------

template <typename TSetH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSetH>>,
                                   Is<ContainerConcept<typename Value<TSetH>::Type>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<ContainerConcept<typename Value<TSetV>::Type>>>
                               >::VALUE,
                          int> = 0>
inline auto
_alignWrapperSimd(TSetH const & stringsH,
                  TSetV const & stringsV,
                  Score<TScoreValue, TScoreSpec> const & scoringScheme,
                  TAlignConfig const & config,
                  TGapModel const & /*gaps*/)
{
    typedef typename SimdVector<int16_t>::Type TSimdAlign;

    unsigned const numAlignments = length(stringsV);
    unsigned const sizeBatch = LENGTH<TSimdAlign>::VALUE;
    unsigned const fullSize = sizeBatch * ((numAlignments + sizeBatch - 1) / sizeBatch);

    String<TScoreValue> results;
    resize(results, numAlignments);

    StringSet<String<Nothing> > trace;  // We need to declare it, but it will not be used.

    // Create a SIMD scoring scheme.
    Score<TSimdAlign, ScoreSimdWrapper<Score<TScoreValue, TScoreSpec> > > simdScoringScheme(scoringScheme);

    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        TSimdAlign resultsBatch;
        if (SEQAN_UNLIKELY(numAlignments < pos + sizeBatch))
        {
            StringSet<std::remove_const_t<typename Value<TSetH>::Type>, Dependent<> > depSetH;
            StringSet<std::remove_const_t<typename Value<TSetV>::Type>, Dependent<> > depSetV;
            for (unsigned i = pos; i < fullSize; ++i)
            {
                if (i >= numAlignments)
                {
                    appendValue(depSetH, back(stringsH));
                    appendValue(depSetV, back(stringsV));
                }
                else
                {
                    appendValue(depSetH, stringsH[i]);
                    appendValue(depSetV, stringsV[i]);
                }
            }
            SEQAN_ASSERT_EQ(length(depSetH), sizeBatch);
            SEQAN_ASSERT_EQ(length(depSetV), sizeBatch);

            _prepareAndRunSimdAlignment(resultsBatch, trace, depSetH, depSetV, simdScoringScheme, config, TGapModel());
        }
        else
        {
            auto infSetH = infixWithLength(stringsH, pos, sizeBatch);
            auto infSetV = infixWithLength(stringsV, pos, sizeBatch);

            _prepareAndRunSimdAlignment(resultsBatch, trace, infSetH, infSetV, simdScoringScheme, config, TGapModel());
        }

        // TODO(rrahn): Could be parallelized!
        for(auto x = pos; x < pos + sizeBatch && x < numAlignments; ++x)
            results[x] = resultsBatch[x - pos];
    }
    return results;
}

// ----------------------------------------------------------------------------
// Function _alignWrapperSimd(); Score; String vs. StringSet
// ----------------------------------------------------------------------------

template <typename TSeqH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSeqH>>,
                                   Not<Is<ContainerConcept<typename Value<TSeqH>::Type>>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<ContainerConcept<typename Value<TSetV>::Type>>>
                              >::VALUE,
                           int> = 0>
inline auto
_alignWrapperSimd(TSeqH const & stringH,
                  TSetV const & stringsV,
                  Score<TScoreValue, TScoreSpec> const & scoringScheme,
                  TAlignConfig const & config,
                  TGapModel const & /*gaps*/)
{
    typedef typename SimdVector<int16_t>::Type TSimdAlign;

    unsigned const numAlignments = length(stringsV);
    unsigned const sizeBatch = LENGTH<TSimdAlign>::VALUE;
    unsigned const fullSize = sizeBatch * ((numAlignments + sizeBatch - 1) / sizeBatch);

    String<TScoreValue> results;
    resize(results, numAlignments);

    // Prepare strings.
    StringSet<TSeqH, Dependent<> > setH;
    for (auto i = 0u; i < sizeBatch; ++i)
        appendValue(setH, stringH);

    StringSet<String<Nothing> > trace;  // We need to declare it, but it will not be used.

    // Create a SIMD scoring scheme.
    Score<TSimdAlign, ScoreSimdWrapper<Score<TScoreValue, TScoreSpec> > > simdScoringScheme(scoringScheme);

    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        TSimdAlign resultsBatch;
        if (SEQAN_UNLIKELY(numAlignments < pos + sizeBatch))
        {
            StringSet<std::remove_const_t<typename Value<TSetV>::Type>, Dependent<> > depSetV;
            for (unsigned i = pos; i < fullSize; ++i)
            {
                if (i >= numAlignments)
                    appendValue(depSetV, back(stringsV));
                else
                    appendValue(depSetV, stringsV[i]);
            }
            SEQAN_ASSERT_EQ(length(depSetV), sizeBatch);

            _prepareAndRunSimdAlignment(resultsBatch, trace, setH, depSetV, simdScoringScheme, config, TGapModel());
        }
        else
        {
            auto infSetV = infixWithLength(stringsV, pos, sizeBatch);
            _prepareAndRunSimdAlignment(resultsBatch, trace, setH, infSetV, simdScoringScheme, config, TGapModel());
        }
        // TODO(rrahn): Could be parallelized!
        for(auto x = pos; x < pos + sizeBatch && x < numAlignments; ++x)
            results[x] = resultsBatch[x - pos];
    }
    return results;
}

// ----------------------------------------------------------------------------
// Function _alignWrapperSimd(); Gaps
// ----------------------------------------------------------------------------

template <typename TSetH,
          typename TSetV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TGapModel,
          std::enable_if_t<And<And<Is<ContainerConcept<TSetH>>,
                                   Is<AlignedSequenceConcept<typename Value<TSetH>::Type>>>,
                               And<Is<ContainerConcept<TSetV>>,
                                   Is<AlignedSequenceConcept<typename Value<TSetV>::Type>>>
                              >::VALUE,
                          int> = 0>
inline auto
_alignWrapperSimd(TSetH & gapSeqSetH,
                  TSetV & gapSeqSetV,
                  Score<TScoreValue, TScoreSpec> const & scoringScheme,
                  TAlignConfig const & config,
                  TGapModel const & /*gaps*/)
{
    typedef typename Value<TSetH>::Type                                 TGapSequenceH;
    typedef typename Value<TSetV>::Type                                 TGapSequenceV;
    typedef typename Size<TGapSequenceH>::Type                          TSize;
    typedef typename Position<TGapSequenceH>::Type                      TPosition;
    typedef TraceSegment_<TPosition, TSize>                             TTraceSegment;

    typedef typename SimdVector<int16_t>::Type                          TSimdAlign;

#if SEQAN_ALIGN_SIMD_PROFILE
    timer = sysTime();
#endif

    unsigned const numAlignments = length(gapSeqSetH);
    unsigned const sizeBatch = LENGTH<TSimdAlign>::VALUE;
    unsigned const fullSize = sizeBatch * ((numAlignments + sizeBatch - 1) / sizeBatch);

    String<TScoreValue> results;
    resize(results, numAlignments);

    // Create a SIMD scoring scheme.
    Score<TSimdAlign, ScoreSimdWrapper<Score<TScoreValue, TScoreSpec> > > simdScoringScheme(scoringScheme);

    // Prepare string sets with sequences.
    StringSet<std::remove_const_t<typename Source<TGapSequenceH>::Type>, Dependent<> > depSetH;
    StringSet<std::remove_const_t<typename Source<TGapSequenceV>::Type>, Dependent<> > depSetV;
    reserve(depSetH, fullSize);
    reserve(depSetV, fullSize);
    for (unsigned i = 0; i < fullSize; ++i)
    {
        if (i >= numAlignments)
        {
            appendValue(depSetH, source(back(gapSeqSetH)));
            appendValue(depSetV, source(back(gapSeqSetV)));
        }
        else
        {
            appendValue(depSetH, source(gapSeqSetH[i]));
            appendValue(depSetV, source(gapSeqSetV[i]));
        }
    }

    // Run alignments in batches.
    for (auto pos = 0u; pos < fullSize; pos += sizeBatch)
    {
        auto infSetH = infixWithLength(depSetH, pos, sizeBatch);
        auto infSetV = infixWithLength(depSetV, pos, sizeBatch);

        TSimdAlign resultsBatch;

        StringSet<String<TTraceSegment> > trace;
        resize(trace, sizeBatch, Exact());

        _prepareAndRunSimdAlignment(resultsBatch, trace, infSetH, infSetV, simdScoringScheme, config, TGapModel());

        // copy results and finish traceback
        // TODO(rrahn): Could be parallelized!
        // to for_each call
        for(auto x = pos; x < pos + sizeBatch && x < numAlignments; ++x)
        {
            results[x] = resultsBatch[x - pos];
            _adaptTraceSegmentsTo(gapSeqSetH[x], gapSeqSetV[x], trace[x - pos]);
        }
#if SEQAN_ALIGN_SIMD_PROFILE
        profile.traceTimer += sysTime() - timer;
        timer = sysTime();
#endif
    }
    return results;
}

}  // namespace seqan
#endif  // #ifndef INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_
