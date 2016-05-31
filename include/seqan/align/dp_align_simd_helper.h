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

// ============================================================================
// Forwards
// ============================================================================

template <unsigned LENGTH>
struct VectorLength_;

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
// Function _checkAndCreateSimdRepresentation()
// ----------------------------------------------------------------------------

template<typename TStringSetH, typename TStringSetV, typename TSimdString>
void inline
_createSimdRepresentation(TStringSetH const & stringsH,
                          TStringSetV const & stringsV,
                          TSimdString & simdH,
                          TSimdString & simdV,
                          SimdAlignEqualLength const & /*tag*/)
{
    resize(simdH, length(stringsH[0]));
    resize(simdV, length(stringsV[0]));
    _createSimdRepImpl(simdH, stringsH);
    _createSimdRepImpl(simdV, stringsV);
}

template<typename TStringSetH,
         typename TStringSetV,
         typename TSimdString,
         typename TPosString,
         typename TSimdVector>
void inline
_createSimdRepresentation(TStringSetH const & stringsH,
                          TStringSetV const & stringsV,
                          TSimdString & simdH,
                          TSimdString & simdV,
                          TSimdString & masksH,
                          TSimdString & masksV,
                          TSimdString & masks,
                          TPosString & endsH,
                          TPosString & endsV,
                          SimdAlignVariableLength<TSimdVector> const & /*tag*/)
{
    using TStringH = typename std::decay<typename Value<TStringSetH>::Type>::type;
    using TStringV = typename std::decay<typename Value<TStringSetV>::Type>::type;
    // check if all sequences have the same length
    unsigned int numAlignments = LENGTH<TSimdVector>::VALUE;

    // otherwise we have to copy the sequences to be able to add a
    // padding character before calling _createSimdRepresentation,
    // because all sequences must have the same length
    auto maxLengthLambda = [](auto& seqLhs, auto& seqRhs) { return length(seqLhs) < length(seqRhs); };
    size_t maxH = length(*std::max_element(begin(stringsH, Standard()), end(stringsH, Standard()), maxLengthLambda));
    size_t maxV = length(*std::max_element(begin(stringsV, Standard()), end(stringsV, Standard()), maxLengthLambda));

    // and we have to prepare the bit masks of the DPScoutState
    resize(masks, maxV, createVector<TSimdVector>(0));
    resize(masksV, maxV, createVector<TSimdVector>(0));
    resize(masksH, maxH, createVector<TSimdVector>(0));

    // copy strings and add padding chars
    StringSet<TStringH> paddedH;
    StringSet<TStringV> paddedV;
    resize(paddedH, numAlignments);
    resize(paddedV, numAlignments);

    for(unsigned i = 0; i < numAlignments; ++i)
    {
        // add padding: the padding character should be part of the amino acid alphabet
        // otherwise a possible score matrix look-up fails, we use 'A' therefor
        paddedH[i] = stringsH[i];
        paddedV[i] = stringsV[i];
        resize(paddedH[i], maxH, 'A');
        resize(paddedV[i], maxV, 'A');

        // mark the original end position of the alignment in the masks (with -1, all bits set)
        assignValue(masksH[length(stringsH[i]) - 1], i, -1);
        assignValue(masksV[length(stringsV[i]) - 1], i, -1);
        endsH.push_back(length(stringsH[i]) - 1);
        endsV.push_back(length(stringsV[i]) - 1);
    }

    // sort the end positions, remove duplicates
    std::sort(endsH.begin(), endsH.end());
    std::unique(endsH.begin(), endsH.end());
    std::sort(endsV.begin(), endsV.end());
    std::unique(endsV.begin(), endsV.end());

    // now create SIMD representation
    resize(simdH, maxH);
    resize(simdV, maxV);
    _createSimdRepImpl(simdH, paddedH);
    _createSimdRepImpl(simdV, paddedV);
}

template<typename TStringSetH, typename TStringSetV, typename TSimdString>
void inline
_checkAndCreateSimdRepresentation(TStringSetH const & stringsH,
                                  TStringSetV const & stringsV,
                                  TSimdString & simdH,
                                  TSimdString & simdV,
                                  TSimdString & masksH,
                                  TSimdString & masksV,
                                  TSimdString & masks,
                                  std::vector<size_t> & endsH,
                                  std::vector<size_t> & endsV)
{
    using TSimdVector = typename Value<TSimdString>::Type;

    auto seqLengthH = length(stringsH[0]);
    auto seqLengthV = length(stringsV[0]);
    auto zipView = makeZipView(stringsH, stringsV);
    bool allSameLength = std::all_of(begin(zipView, Standard()), end(zipView, Standard()),
                                     [seqLengthH, seqLengthV](auto param)
                                     {
                                         return (length(std::get<0>(param)) == seqLengthH) &&
                                                (length(std::get<1>(param)) == seqLengthV);
                                     });

    // if yes, create SIMD representation without doing anything else
    if(!allSameLength)
        _createSimdRepresentation(stringsH, stringsV, simdH, simdV, masksH, masksV, masks, endsH, endsV,
                                  SimdAlignVariableLength<TSimdVector>());
    else
        _createSimdRepresentation(stringsH, stringsV, simdH, simdV, SimdAlignEqualLength());
}

// ----------------------------------------------------------------------------
// Function _createSimdRepImpl()
// ----------------------------------------------------------------------------

#define SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos, chrPos)    data[strPos][chrPos], data[strPos + 1][chrPos]
#define SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos, chrPos)    SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_2(data, strPos + 2, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos, chrPos)    SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_4(data, strPos + 4, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos, chrPos)   SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos, chrPos),  SEQAN_CREATE_SIMD_REP_IMPL_8(data, strPos + 8, chrPos)
#define SEQAN_CREATE_SIMD_REP_IMPL_32(data, strPos, chrPos)   SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos, chrPos), SEQAN_CREATE_SIMD_REP_IMPL_16(data, strPos + 16, chrPos)

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
                            TGapModel const & /*gapModel*/)
{
    String<TResult, Alloc<OverAligned> > stringSimdH;
    String<TResult, Alloc<OverAligned> > stringSimdV;
    String<TResult, Alloc<OverAligned> > masksH;
    String<TResult, Alloc<OverAligned> > masksV;
    String<TResult, Alloc<OverAligned> > masks;

    std::vector<decltype(length(seqH))> endsH;
    std::vector<decltype(length(seqV))> endsV;

    // create the SIMD representation of the alignments
    // in case of a variable length alignment the variables masks, endsH, endsV will be filled
    _checkAndCreateSimdRepresentation(seqH, seqV, stringSimdH, stringSimdV, masksH, masksV, masks, endsH, endsV);

    // if alignments have equal dimensions do nothing
    if(endsH.size() == 0)
    {
        DPScoutState_<SimdAlignEqualLength> dpScoutState;
        results = _setUpAndRunAlignment(traces, dpScoutState, stringSimdH, stringSimdV,
                                        scoringScheme, alignConfig, TGapModel());
    }
    else  // otherwise prepare the special DPScoutState
    {
        DPScoutState_<SimdAlignVariableLength<TResult> > dpScoutState;
        dpScoutState.dimV = length(stringSimdV);
        dpScoutState.isLocalAlignment = IsLocalAlignment_<TAlgo>::VALUE;
        dpScoutState.right = IsFreeEndGap_<TFreeEndGaps, DPLastColumn>::VALUE;
        dpScoutState.bottom = IsFreeEndGap_<TFreeEndGaps, DPLastRow>::VALUE;
        swap(dpScoutState.masksH, masksH);
        swap(dpScoutState.masksV, masksV);
        swap(dpScoutState.masks, masks);
        std::swap(dpScoutState.endsH, endsH);
        std::swap(dpScoutState.endsV, endsV);
        results = _setUpAndRunAlignment(traces, dpScoutState, stringSimdH, stringSimdV,
                                        scoringScheme, alignConfig, TGapModel());
    }
}

}  // namespace seqan
#endif  // #ifndef INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_
