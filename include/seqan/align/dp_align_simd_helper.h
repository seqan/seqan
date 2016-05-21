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

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

template <typename T>
struct IsAlignObject_ : False
{};

template <typename TSeq, typename TSpec>
struct IsAlignObject_<Align<TSeq, TSpec> > : True
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _checkAndCreateSimdRepresentation()
// ----------------------------------------------------------------------------

template<typename TStringSetH, typename TStringSetV, typename TSimdString>
void _checkAndCreateSimdRepresentation(TStringSetH const & stringsH,
                                       TStringSetV const & stringsV,
                                       TSimdString & simdH,
                                       TSimdString & simdV,
                                       TSimdString & masksH,
                                       TSimdString & masksV,
                                       TSimdString & masks,
                                       std::vector<size_t> & endsH,
                                       std::vector<size_t> & endsV)
{
    using TStringH = typename Value<TStringSetH>::Type;
    using TStringV = typename Value<TStringSetV>::Type;
    using TSimdVector = typename Value<TSimdString>::Type;
    // check if all sequences have the same length
    unsigned int numAlignments = LENGTH<TSimdVector>::VALUE;
    bool allEqualH = true, allEqualV = true;
    std::vector<size_t> seqLengthsH(numAlignments), seqLengthsV(numAlignments);
    seqLengthsH[0] = length(stringsH[0]);
    seqLengthsV[0] = length(stringsV[0]);

    for(unsigned i = 1; i < numAlignments; ++i)
    {
        seqLengthsH[i] = length(stringsH[i]);
        seqLengthsV[i] = length(stringsV[i]);
        allEqualH &= (seqLengthsH[i] == seqLengthsH[i-1]);
        allEqualV &= (seqLengthsV[i] == seqLengthsV[i-1]);
    }

    // if yes, create SIMD representation without doing anything else
    if(allEqualH && allEqualV)
    {
        resize(simdH, seqLengthsH[0]);
        resize(simdV, seqLengthsV[0]);
        _createSimdRepresentation(simdH, stringsH, seqLengthsH[0]);
        _createSimdRepresentation(simdV, stringsV, seqLengthsV[0]);
        return;
    }

    // otherwise we have to copy the sequences to be able to add a
    // padding character before calling _createSimdRepresentation,
    // because all sequences must have the same length
    size_t maxH = *std::max_element(seqLengthsH.begin(), seqLengthsH.end());
    size_t maxV = *std::max_element(seqLengthsV.begin(), seqLengthsV.end());

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
        assignValue(paddedH, i, stringsH[i]);
        assignValue(paddedV, i, stringsV[i]);
        resize(paddedH[i], maxH, 'A');
        resize(paddedV[i], maxV, 'A');

        // mark the original end position of the alignment in the masks (with -1, all bits set)
        assignValue(masksH[seqLengthsH[i]-1], i, -1);
        assignValue(masksV[seqLengthsV[i]-1], i, -1);
        endsH.push_back(seqLengthsH[i]-1);
        endsV.push_back(seqLengthsV[i]-1);
    }

    // sort the end positions, remove duplicates
    std::sort(endsH.begin(), endsH.end());
    endsH.erase( std::unique( endsH.begin(), endsH.end() ), endsH.end() );
    std::sort(endsV.begin(), endsV.end());
    endsV.erase( std::unique( endsV.begin(), endsV.end() ), endsV.end() );

    // now create SIMD representation
    resize(simdH, maxH);
    resize(simdV, maxV);
    _createSimdRepresentation(simdH, paddedH, maxH);
    _createSimdRepresentation(simdV, paddedV, maxV);
}

template<typename TContainer, typename TSimdString, typename TSize>
inline SEQAN_FUNC_ENABLE_IF(IsAlignObject_<typename Value<TContainer>::Type>, void)
_checkAndCreateSimdRepresentation(TContainer & align,
                                  TSimdString & simdH,
                                  TSimdString & simdV,
                                  TSimdString & masksH,
                                  TSimdString & masksV,
                                  TSimdString & masks,
                                  std::vector<TSize> & endsH,
                                  std::vector<TSize> & endsV)
{
    using TAlign    = typename Value<TContainer>::Type;
    using TSequence = typename Source<TAlign>::Type;
    using TSimdVector = typename Value<TSimdString>::Type;

    // check if all sequences have the same length
    unsigned numAlignments = LENGTH<TSimdVector>::VALUE;

    StringSet<TSequence, Dependent<Tight> > depSet1;
    StringSet<TSequence, Dependent<Tight> > depSet2;
    reserve(depSet1, numAlignments);
    reserve(depSet2, numAlignments);
    for (auto& obj : align)
    {
        appendValue(depSet1, source(row(obj, 0)));
        appendValue(depSet2, source(row(obj, 1)));
    }

    _checkAndCreateSimdRepresentation(depSet1, depSet2, simdH, simdV, masksH, masksV, masks, endsH, endsV);
}

// ----------------------------------------------------------------------------
// Function _createSimdRepresentation()
// ----------------------------------------------------------------------------

template<typename TSimdVector, typename TSpec, typename TContainer>
inline SEQAN_FUNC_ENABLE_IF(Is<ContainerConcept<typename Value<TContainer>::Type> >, void)
_createSimdRepresentation(String<TSimdVector, TSpec> & simdRepr,
                          TContainer const & strings,
                          size_t stringLength)
{
    // TODO(rrahn): Make code generic!
    switch((int)LENGTH<TSimdVector>::VALUE)
    {
        case 8:
            for(size_t x = 0; x < stringLength; ++x)
                fillVector(simdRepr[x], strings[0][x], strings[1][x], strings[2][x], strings[3][x],
                           strings[4][x], strings[5][x], strings[6][x], strings[7][x]);
                break;
        case 16:
            for(size_t x = 0; x < stringLength; ++x)
                fillVector(simdRepr[x], strings[0][x], strings[1][x], strings[2][x], strings[3][x],
                           strings[4][x], strings[5][x], strings[6][x], strings[7][x],
                           strings[8][x], strings[9][x], strings[10][x], strings[11][x],
                           strings[12][x], strings[13][x], strings[14][x], strings[15][x]);
                break;
    }
}

template<typename TSimdVector, typename TSpec, typename TString>
inline SEQAN_FUNC_DISABLE_IF(Is<ContainerConcept<typename Value<TString>::Type> >, void)
_createSimdRepresentation(String<TSimdVector, TSpec> & simdRepr,
                          TString const & seq,
                          size_t stringLength)
{
    for(size_t x = 0; x < stringLength; ++x)
        fillVector(simdRepr[x], seq[x]);
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

}
#endif  // #ifndef INCLUDE_SEQAN_ALIGN_DP_ALIGN_SIMD_HELPER_H_
