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

#define SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos, chrPos)    data[strPos + 1][chrPos], data[strPos][chrPos]
#define SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos, chrPos)    SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos + 2, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos, chrPos)    SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos + 4, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos, chrPos)   SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos + 8, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_32(data, strPos, chrPos)   SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos + 16, chrPos), SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos, chrPos)

#define SEQAN_CREATE_SIMD_REP_FILL_IMPL_2(MACRO, data, chrPos) MACRO(data, 0, chrPos)
#define SEQAN_CREATE_SIMD_REP_FILL_IMPL(data, chrPos, SIZE) SEQAN_CREATE_SIMD_REP_FILL_IMPL_2(SEQAN_CREATE_SIMD_REP_IMPL_##SIZE, data, chrPos)

#define SEQAN_CREATE_SIMD_REP_IMPL(SIZE)                                                \
template <typename TSimdVecs, typename TStrings>                                        \
inline void _createSimdRepImpl(TSimdVecs & simdStr,                                     \
                               TStrings const & strings,                                \
                               VectorLength_<SIZE> const & /*size*/)                    \
{                                                                                       \
    auto itB = begin(simdStr, Standard());                                              \
    auto itE = end(simdStr, Standard());                                                \
    for (auto it = itB; it != itE; ++it)                                                \
        fillVector(*it, SEQAN_CREATE_SIMD_REP_FILL_IMPL(strings, it - itB, SIZE));      \
}

SEQAN_CREATE_SIMD_REP_IMPL(2)
SEQAN_CREATE_SIMD_REP_IMPL(4)
SEQAN_CREATE_SIMD_REP_IMPL(8)
SEQAN_CREATE_SIMD_REP_IMPL(16)
SEQAN_CREATE_SIMD_REP_IMPL(32)

template <typename TSimdVecs, typename TStrings>
inline void _createSimdRepImpl(TSimdVecs & simdStr,
                               TStrings const & strings)
{
    _createSimdRepImpl(simdStr, strings, VectorLength_<LENGTH<typename Value<TSimdVecs>::Type>::VALUE>());
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
                            TGapModel const & /*gapModel*/,
                            SimdAlignEqualLength const & /*tag*/)
{
    String<TResult, Alloc<OverAligned> > stringSimdH;
    String<TResult, Alloc<OverAligned> > stringSimdV;

    resize(stringSimdH, length(seqH[0]));
    resize(stringSimdV, length(seqV[0]));
    _createSimdRepImpl(stringSimdH, seqH);
    _createSimdRepImpl(stringSimdV, seqV);

    DPScoutState_<SimdAlignEqualLength> state;
    results = _setUpAndRunAlignment(traces, state, stringSimdH, stringSimdV, scoringScheme, alignConfig, TGapModel());
}

template <typename TResult,
          typename TTraces,
          typename TSequencesH,
          typename TSequencesV,
          typename TScore,
          typename TAlgo, typename TBand, typename TFreeEndGaps, typename TTraceback,
          typename TGapModel,
          typename TTraits>
inline void
_prepareAndRunSimdAlignment(TResult & results,
                            TTraces & traces,
                            TSequencesH const & seqH,
                            TSequencesV const & seqV,
                            TScore const & scoringScheme,
                            AlignConfig2<TAlgo, TBand, TFreeEndGaps, TTraceback> const & alignConfig,
                            TGapModel const & /*gapModel*/,
                            SimdAlignVariableLength<TTraits> const /*tag*/)
{
    SEQAN_ASSERT_EQ(length(seqH), length(seqV));
    SEQAN_ASSERT_EQ(static_cast<decltype(length(seqH))>(LENGTH<TResult>::VALUE), length(seqH));

    using TPadStringH = ModifiedString<typename Value<TSequencesH>::Type, ModPadding>;
    using TPadStringV = ModifiedString<typename Value<TSequencesV>::Type, ModPadding>;

    String<TResult, Alloc<OverAligned> > stringSimdH;
    String<TResult, Alloc<OverAligned> > stringSimdV;

    DPScoutState_<SimdAlignVariableLength<SimdAlignVariableLengthTraits<TResult, TSequencesH, TSequencesV> > > state;

    String<size_t> lengthsH;
    String<size_t> lengthsV;

    resize(lengthsH, length(seqH));
    resize(lengthsV, length(seqV));
    resize(state.endsH, length(seqH));
    resize(state.endsV, length(seqV));

    for (unsigned i = 0; i < length(seqH); ++i)
    {
        lengthsH[i] = length(seqH[i]) - 1;
        lengthsV[i] = length(seqV[i]) - 1;
        state.endsH[i] = i;
        state.endsV[i] = i;
    }

//    std::cout << "lengthsH: " <<std::endl;
//    for (auto p : lengthsH)
//        std::cout << p << " ";
//    std::cout << std::endl;
//
//    std::cout << "lengthsV: " <<std::endl;
//    for (auto p : lengthsV)
//        std::cout << p << " ";
//    std::cout << std::endl;

    setHost(state.sortedEndsH, lengthsH);
    setHost(state.sortedEndsV, lengthsV);
    setCargo(state.sortedEndsH, state.endsH);
    setCargo(state.sortedEndsV, state.endsV);

    auto maxLengthLambda = [](auto& lengthLhs, auto& lengthRhs) { return lengthLhs < lengthRhs; };
    sort(state.sortedEndsH, maxLengthLambda, Serial());
    sort(state.sortedEndsV, maxLengthLambda, Serial());

//    std::cout << "state.sortedEndsH: " <<std::endl;
//    for (auto p : state.sortedEndsH)
//        std::cout << p << " ";
//    std::cout << std::endl;
//
//    std::cout << "state.sortedEndsV: " <<std::endl;
//    for (auto p : state.sortedEndsV)
//        std::cout << p << " ";
//    std::cout << std::endl;

    size_t maxH = back(state.sortedEndsH) + 1;
    size_t maxV = back(state.sortedEndsV) + 1;

    // and we have to prepare the bit masks of the DPScoutState
    resize(state.masks,  maxV, createVector<TResult>(0));
    resize(state.masksV, maxV, createVector<TResult>(0));
    resize(state.masksH, maxH, createVector<TResult>(0));

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
        assignValue(state.masksH[lengthsH[i]], i, -1);
        assignValue(state.masksV[lengthsV[i]], i, -1);
    }

    // now create SIMD representation
    resize(stringSimdH, maxH);
    resize(stringSimdV, maxV);
    _createSimdRepImpl(stringSimdH, paddedH);
    _createSimdRepImpl(stringSimdV, paddedV);

    state.dimV = length(stringSimdV);
    state.isLocalAlignment = IsLocalAlignment_<TAlgo>::VALUE;
    state.right = IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE;
    state.bottom = IsFreeEndGap_<TFreeEndGaps, DPLastRow>::VALUE;

    results = _setUpAndRunAlignment(traces, state, stringSimdH, stringSimdV, scoringScheme, alignConfig, TGapModel());
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
    if(allSameLength)
        _prepareAndRunSimdAlignment(results, traces, seqH, seqV, scoringScheme, alignConfig, TGapModel(), SimdAlignEqualLength());
    else
        _prepareAndRunSimdAlignment(results, traces, seqH, seqV, scoringScheme, alignConfig, TGapModel(),
                                    SimdAlignVariableLength<Nothing>());
}

}  // namespace seqan
#endif  // #ifndef INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_

