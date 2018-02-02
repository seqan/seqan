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
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================
// generic SIMD interface for AVX512
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX512_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX512_H_

namespace seqan {

// SimdParams_<64, 64>: 512bit = 64 elements * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector64Char,     char,           64)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector64SChar,    signed char,    64)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector64UChar,    unsigned char,  64)

// SimdParams_<64, 32>: 512bit = 32 elements * 2 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Short,    short,          64)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32UShort,   unsigned short, 64)

// SimdParams_<64, 16>: 512bit = 16 elements * 4 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Int,      int,            64)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UInt,     unsigned int,   64)

// SimdParams_<64, 8>: 512bit = 8 elements * 8 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Int64,     int64_t,        64)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UInt64,    uint64_t,       64)

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// AVX512 wrappers (512bit vectors)
// ============================================================================

// --------------------------------------------------------------------------
// _fillVector (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L, typename TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue> const & x,
            std::index_sequence<0> const &,
            SimdParams_<64, L>)
{
    vector = createVector<TSimdVector>(std::get<0>(x));
}

template <typename TSimdVector, int L, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args,
            std::index_sequence<INDICES...> const &,
            SimdParams_<64, L>)
{
    using TSimdValue = typename Value<TSimdVector>::Type;
    vector = TSimdVector{static_cast<TSimdValue>(std::get<INDICES>(args))...};
}

// --------------------------------------------------------------------------
// _clearVector (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector & vector, SimdParams_<64, L>)
{
    vector = TSimdVector{};
}

// --------------------------------------------------------------------------
// _createVector (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue, int L>
inline TSimdVector _createVector(TValue const x, SimdParams_<64, L>)
{
    using TValue_ = typename Value<TSimdVector>::Type;
    return TSimdVector{} + static_cast<TValue_>(x);
}

// --------------------------------------------------------------------------
// _cmpEq (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _cmpEq(TSimdVector & a, TSimdVector & b, SimdParams_<64, L>)
{
    return a == b;
}

// bad auto-vectorization for gcc
#ifndef __AVX512BW__
template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<64, 32>)
{
    auto aLow = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, a), 0);
    auto bLow = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, b), 0);
    auto cmpLow = _mm256_cmpeq_epi16(aLow, bLow);

    auto aHigh = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, a), 1);
    auto bHigh = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, b), 1);
    auto cmpHigh = _mm256_cmpeq_epi16(aHigh, bHigh);

    auto result = _mm512_broadcast_i64x4(cmpLow);
    result = _mm512_inserti64x4(result, cmpHigh, 1);
    return SEQAN_VECTOR_CAST_(TSimdVector, result);
}
#endif

// --------------------------------------------------------------------------
// _cmpGt (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L, typename TValue>
inline TSimdVector _cmpGt(TSimdVector & a, TSimdVector & b, SimdParams_<64, L, TValue>)
{
    return a > b;
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAndNot(TSimdVector & a, TSimdVector & b, SimdParams_<64, L>)
{
    return (~a & b);
}

// --------------------------------------------------------------------------
// _max (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L, typename TValue>
inline TSimdVector _max(TSimdVector & a, TSimdVector & b, SimdParams_<64, L, TValue>)
{
    return a > b ? a : b;
}

// --------------------------------------------------------------------------
// _min (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L, typename TValue>
inline TSimdVector _min(TSimdVector & a, TSimdVector & b, SimdParams_<64, L, TValue>)
{
    return a < b ? a : b;
}

// --------------------------------------------------------------------------
// _blend (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, int L>
inline TSimdVector _blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask, SimdParams_<64, L>)
{
    return mask ? b : a;
}

// bad auto-vectorization for gcc
#ifndef __AVX512BW__
template <typename TSimdVector, typename TSimdVectorMask>
inline TSimdVector _blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask, SimdParams_<64, 32>)
{
    auto aLow = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, a), 0);
    auto bLow = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, b), 0);
    auto maskLow = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, mask), 0);
    auto blendLow = _mm256_blendv_epi8(aLow, bLow, maskLow);

    auto aHigh = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, a), 1);
    auto bHigh = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, b), 1);
    auto maskHigh = _mm512_extracti64x4_epi64(SEQAN_VECTOR_CAST_(const __m512i&, mask), 1);
    auto blendHigh = _mm256_blendv_epi8(aHigh, bHigh, maskHigh);

    auto result = _mm512_broadcast_i64x4(blendLow);
    result = _mm512_inserti64x4(result, blendHigh, 1);
    return SEQAN_VECTOR_CAST_(TSimdVector, result);
}
#endif

// --------------------------------------------------------------------------
// _storeu (512bit)
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector, int L>
inline void _storeu(T * memAddr, TSimdVector & vec, SimdParams_<64, L>)
{
    constexpr auto length = LENGTH<TSimdVector>::VALUE;
    for (unsigned i = 0; i < length; i++)
        memAddr[i] = vec[i];
}

// ----------------------------------------------------------------------------
// Function _load() 512bit
// ----------------------------------------------------------------------------

template <typename TSimdVector, typename T, int L>
inline TSimdVector _load(T const * memAddr, SimdParams_<64, L>)
{
    constexpr auto length = LENGTH<TSimdVector>::VALUE;
    TSimdVector result;
    for (unsigned i = 0; i < length; i++)
        result[i] = memAddr[i];
    return result;
}

// --------------------------------------------------------------------------
// _shiftRightLogical (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<64, L>)
{
    return vector >> imm;
}

// --------------------------------------------------------------------------
// _gather (512bit)
// --------------------------------------------------------------------------

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE, int L>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & /*scale*/,
        SimdParams_<64, L>)
{
    constexpr auto length = LENGTH<TSimdVector>::VALUE;
    TSimdVector result;
    for (unsigned i = 0; i < length; i++)
        result[i] = memAddr[idx[i]];
    return result;
}

// --------------------------------------------------------------------------
// _shuffleVector (512bit)
// --------------------------------------------------------------------------

template <typename TSimdVector1, typename TSimdVector2, int L>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<64, L>, SimdParams_<64, 64>)
{
    constexpr auto length = seqan::LENGTH<TSimdVector1>::VALUE;
    TSimdVector1 result{};
    for(unsigned i = 0u; i < length; ++i)
        result[i] = vector[indices[i]];
    return result;
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX512_H_
