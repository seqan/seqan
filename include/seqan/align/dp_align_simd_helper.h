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

template <typename TSimdVector_, typename TSeqH_, typename TSeqV_, typename TDPProfile_>
struct SimdAlignVariableLengthTraits
{
    using TSimdVector   = TSimdVector_;
    using TSeqH         = TSeqH_;
    using TSeqV         = TSeqV_;
    using TDPProfile    = TDPProfile_;
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
                      String<size_t> & lengthsH,
                      String<size_t> & lengthsV,
                      DPScoutState_<SimdAlignVariableLength<TTraits> > & state)
{
    SEQAN_ASSERT_EQ(length(seqH), length(seqV));
    SEQAN_ASSERT_EQ(static_cast<decltype(length(seqH))>(LENGTH<typename Value<TStringSimdH>::Type>::VALUE), length(seqH));

    using TSimdVector = typename TTraits::TSimdVector;
    using TSimdValueType = typename Value<TSimdVector>::Type;

    using TPadStringH = ModifiedString<typename Value<TSequencesH const>::Type, ModPadding>;
    using TPadStringV = ModifiedString<typename Value<TSequencesV const>::Type, ModPadding>;

    resize(lengthsH, length(seqH), Exact{});
    resize(lengthsV, length(seqV), Exact{});

    for (unsigned i = 0; i < length(seqH); ++i)
    {
        lengthsH[i] = length(seqH[i]);
        lengthsV[i] = length(seqV[i]);
    }

    // Sort and remove unique elements from length vectors.
    auto maxLengthLambda = [](auto& lengthLhs, auto& lengthRhs) { return lengthLhs < lengthRhs; };
    std::sort(begin(lengthsH, Standard{}), end(lengthsH, Standard{}), maxLengthLambda);
    std::sort(begin(lengthsV, Standard{}), end(lengthsV, Standard{}), maxLengthLambda);

    erase(lengthsH,
          std::unique(begin(lengthsH, Standard{}), end(lengthsH, Standard{})) - begin(lengthsH, Standard{}),
          length(lengthsH));
    erase(lengthsV,
          std::unique(begin(lengthsV, Standard{}), end(lengthsV, Standard{})) - begin(lengthsV, Standard{}),
          length(lengthsV));

    // Initialize iterator to the lengths vectors.
    state.nextEndsH = begin(lengthsH, Rooted{});
    state.nextEndsV = begin(lengthsV, Rooted{});

    size_t maxH = back(lengthsH);
    size_t maxV = back(lengthsV);

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

        // Store the end points as vector in both dimensions.
        assignValue(state.endPosVecH, i, static_cast<TSimdValueType>(length(seqH[i])));
        assignValue(state.endPosVecV, i, static_cast<TSimdValueType>(length(seqV[i])));
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
        using TDPProfile = typename SetupAlignmentProfile_<TAlgo, TFreeEndGaps, TGapModel, TTraceback>::Type;

        DPScoutState_<SimdAlignVariableLength<SimdAlignVariableLengthTraits<TResult,
                                                                            TSequencesH,
                                                                            TSequencesV,
                                                                            TDPProfile>>> state;
        String<size_t> lengthsH;
        String<size_t> lengthsV;
        _prepareSimdAlignment(stringSimdH, stringSimdV, seqH, seqV, lengthsH, lengthsV, state);

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
    typedef typename SimdVector<TScoreValue>::Type TSimdAlign;

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
    typedef typename SimdVector<TScoreValue>::Type TSimdAlign;

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

    typedef typename SimdVector<TScoreValue>::Type                      TSimdAlign;

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
