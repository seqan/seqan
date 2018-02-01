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
// Author: David Weese <david.weese@fu-berlin.de>
//         René Rahn <rene.rahn@fu-berlin.de>
//         Stefan Budach <stefan.budach@fu-berlin.de>
// ==========================================================================
// generic SIMD interface for SSE3 / AVX2
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX2_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX2_H_

namespace seqan {

// SimdParams_<32, 32>: 256bit = 32 elements * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Char,     char,           32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32SChar,    signed char,    32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32UChar,    unsigned char,  32)

// SimdParams_<32, 16>: 256bit = 16 elements * 2 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Short,    short,          32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UShort,   unsigned short, 32)

// SimdParams_<32, 8>: 256bit = 8 elements * 4 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Int,       int,            32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UInt,      unsigned int,   32)

// SimdParams_<32, 4>: 256bit = 4 elements * 8 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int64,     int64_t,        32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt64,    uint64_t,       32)

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// AVX/AVX2 wrappers (256bit vectors)
// ============================================================================

// --------------------------------------------------------------------------
// _fillVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename ...TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &, SimdParams_<32, 32>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi8(std::get<0>(x)));
}

template <typename TSimdVector, typename ...TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &, SimdParams_<32, 16>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi16(std::get<0>(x)));
}

template <typename TSimdVector, typename ...TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &, SimdParams_<32, 8>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi32(std::get<0>(x)));
}

template <typename TSimdVector, typename ...TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &, SimdParams_<32, 4>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi64x(std::get<0>(x)));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 32>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_setr_epi8(std::get<INDICES>(args)...));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 16>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_setr_epi16(std::get<INDICES>(args)...));
}
template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 8>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_setr_epi32(std::get<INDICES>(args)...));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 4>)
{
    // reverse argument list 0, 1, 2, 3 -> 3, 2, 1, 0
    // NOTE(marehr): Intel linux fails to reverse argument list and only
    // _mm256_set_epi64x has no reverse equivalent
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set_epi64x(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...));
}

// --------------------------------------------------------------------------
// _clearVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector & vector, SimdParams_<32, L>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_setzero_si256());
}

// --------------------------------------------------------------------------
// _createVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi8(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi16(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi32(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_< 32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi64x(x));
}

// --------------------------------------------------------------------------
// _cmpEq (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                             SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _cmpGt (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                             SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, uint8_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x80 = ~0x7F (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpgt_epi8(
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_set1_epi8(~0x7F)),
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, b), _mm256_set1_epi8(~0x7F))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, uint16_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x8000 = ~0x7FFF (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpgt_epi16(
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_set1_epi16(~0x7FFF)),
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, b), _mm256_set1_epi16(~0x7FFF))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, uint32_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x80000000 = ~0x7FFFFFFF (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpgt_epi32(
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_set1_epi32(~0x7FFFFFFF)),
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, b), _mm256_set1_epi32(~0x7FFFFFFF))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, int64_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, uint64_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x8000000000000000ul = ~0x7FFFFFFFFFFFFFFFul (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpgt_epi64(
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, a) ,_mm256_set1_epi64x(~0x7FFFFFFFFFFFFFFFul)),
                                  _mm256_xor_si256(SEQAN_VECTOR_CAST_(const __m256i&, b), _mm256_set1_epi64x(~0x7FFFFFFFFFFFFFFFul))));
}

// --------------------------------------------------------------------------
// _bitwiseOr (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseOr(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_or_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                           SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAnd (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAnd(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_and_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                            SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAndNot(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_andnot_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));

}
template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

// --------------------------------------------------------------------------
// _divide (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi8(a, _mm256_set1_epi8(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi16(a, _mm256_set1_epi16(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi32(a, _mm256_set1_epi32(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi64(a, _mm256_set1_epi64x(b)));
}

// --------------------------------------------------------------------------
// _add (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _sub (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _mult (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const &/*b*/, SimdParams_<32, 32>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("AVX2 intrinsics for multiplying 8 bit values not implemented!");
    return a;
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_mullo_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                 SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_mullo_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                 SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const &/*b*/, SimdParams_<32, 4>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("AVX2 intrinsics for multiplying 64 bit values not implemented!");
    return a;
}

// --------------------------------------------------------------------------
// _max (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, uint8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epu8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, uint16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epu16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, uint32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epu32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, int64_t>)
{
    #if defined(__AVX512VL__)
        return SEQAN_VECTOR_CAST_(TSimdVector,
                                  _mm256_max_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                   SEQAN_VECTOR_CAST_(const __m256i&, b)));
    #else // defined(__AVX512VL__)
        return blend(b, a, cmpGt(a, b));
    #endif // defined(__AVX512VL__)
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, uint64_t>)
{
    #if defined(__AVX512VL__)
        return SEQAN_VECTOR_CAST_(TSimdVector,
                                  _mm256_max_epu64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                   SEQAN_VECTOR_CAST_(const __m256i&, b)));
    #else // defined(__AVX512VL__)
        return blend(b, a, cmpGt(a, b));
    #endif // defined(__AVX512VL__)
}


// --------------------------------------------------------------------------
// _min (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 32, uint8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epu8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 16, uint16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epu16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 8, uint32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_min_epu32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, int64_t>)
{
    #if defined(__AVX512VL__)
        return SEQAN_VECTOR_CAST_(TSimdVector,
                                  _mm256_min_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                   SEQAN_VECTOR_CAST_(const __m256i&, b)));
    #else // defined(__AVX512VL__)
        return blend(a, b, cmpGt(a, b));
    #endif // defined(__AVX512VL__)
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<32, 4, uint64_t>)
{
    #if defined(__AVX512VL__)
        return SEQAN_VECTOR_CAST_(TSimdVector,
                                  _mm256_min_epu64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                   SEQAN_VECTOR_CAST_(const __m256i&, b)));
    #else // defined(__AVX512VL__)
        return blend(a, b, cmpGt(a, b));
    #endif // defined(__AVX512VL__)
}

// --------------------------------------------------------------------------
// _blend (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, int L>
inline TSimdVector _blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_blendv_epi8(SEQAN_VECTOR_CAST_(const __m256i &, a),
                                                 SEQAN_VECTOR_CAST_(const __m256i &, b),
                                                 SEQAN_VECTOR_CAST_(const __m256i &, mask)));
}

// --------------------------------------------------------------------------
// _storeu (256bit)
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector, int L>
inline void _storeu(T * memAddr, TSimdVector const & vec, SimdParams_<32, L>)
{
    _mm256_storeu_si256((__m256i*)memAddr, SEQAN_VECTOR_CAST_(const __m256i&, vec));
}

// ----------------------------------------------------------------------------
// Function _load() 256bit
// ----------------------------------------------------------------------------

template <typename TSimdVector, typename T, int L>
inline TSimdVector _load(T const * memAddr, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_load_si256((__m256i const *) memAddr));
}

// --------------------------------------------------------------------------
// _shiftRightLogical (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi16(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm) & _mm256_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi16(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi32(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi64(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}

// --------------------------------------------------------------------------
// Extend sign from integer types 256bit
// --------------------------------------------------------------------------

inline __m256i
seqan_mm256_i16sign_extend_epis8(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xff00)
            _mm256_sub_epi16( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x80 (select msb)
                    v,
                    _mm256_set1_epi16(0x80)
                ),
                _mm256_set1_epi16(1)
            ),
            _mm256_set1_epi16(static_cast<uint16_t>(0xff00u))
        )
    );
}

inline __m256i
seqan_mm256_i32sign_extend_epis8(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xffffff00u)
            _mm256_sub_epi32( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x80 (select msb)
                    v,
                    _mm256_set1_epi32(0x80)
                ),
                _mm256_set1_epi32(1)
            ),
            _mm256_set1_epi32(static_cast<uint32_t>(0xffffff00u))
        )
    );
}

inline __m256i
seqan_mm256_i32sign_extend_epis16(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xffff0000u)
            _mm256_sub_epi32( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x8000 (select msb)
                    v,
                    _mm256_set1_epi32(0x8000)
                ),
                _mm256_set1_epi32(1)
            ),
            _mm256_set1_epi32(static_cast<uint32_t>(0xffff0000u))
        )
    );
}

inline __m256i
seqan_mm256_i64sign_extend_epis8(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xffffffffffffff00ul)
            _mm256_sub_epi64( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x80 (select msb)
                    v,
                    _mm256_set1_epi64x(0x80)
                ),
                _mm256_set1_epi64x(1)
            ),
            _mm256_set1_epi64x(static_cast<uint64_t>(0xffffffffffffff00ul))
        )
    );
}

inline __m256i
seqan_mm256_i64sign_extend_epis16(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xffffffffffff0000ul)
            _mm256_sub_epi64( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x8000 (select msb)
                    v,
                    _mm256_set1_epi64x(0x8000)
                ),
                _mm256_set1_epi64x(1)
            ),
            _mm256_set1_epi64x(static_cast<uint64_t>(0xffffffffffff0000ul))
        )
    );
}

inline __m256i
seqan_mm256_i64sign_extend_epis32(__m256i const & v)
{
    return _mm256_or_si256( // extend sign (v | hi-bits)
        v,
        _mm256_and_si256( // select hi-bits (hi-bits = msk & 0xffffffffffff0000ul)
            _mm256_sub_epi64( // msk = msb - 1
                _mm256_andnot_si256( //msb = ~v & 0x80000000 (select msb)
                    v,
                    _mm256_set1_epi64x(0x80000000)
                ),
                _mm256_set1_epi64x(1)
            ),
            _mm256_set1_epi64x(static_cast<uint64_t>(0xffffffff00000000ul))
        )
    );
}

// --------------------------------------------------------------------------
// _gather (256bit)
// --------------------------------------------------------------------------

template <typename TValue, typename TSize, TSize SCALE>
inline __m256i
seqan_mm256_i8gather_epi(TValue const * memAddr,
                         __m256i const & idx,
                         std::integral_constant<TSize, SCALE> const & /*scale*/)
{
    // mem:    ( 0,  3,  6,  9 | 12, 15, 18, 21 | 24, 27, 30, 33 | 36, 39, 42, 45 || 48, 51, 54, 57 | 60, 63, 66, 69 | 72, 75, 78, 81 | 84, 87, 90, 93)
    // idx:    (31, 30, 29, 28 | 27, 26, 25, 24 | 23, 22, 21, 20 | 19, 18, 17, 16 || 15, 14, 13, 12 | 11, 10,  9,  8 |  7,  6,  5,  4 |  3,  2,  1,  0)
    // pack:   (93, 90, 87, 84 | 81, 78, 75, 72 | 69, 66, 63, 60 | 57, 54, 51, 48 || 45, 42, 39, 36 | 33, 30, 27, 24 | 21, 18, 15, 12 |  9,  6,  3,  0)
    return _mm256_packus_epi16(
        // pckLow: (93,  0, 90,  0 | 87,  0, 84,  0 | 81,  0, 78,  0 | 75,  0, 72,  0 || 45,  0, 42,  0 | 39,  0, 36,  0 | 33,  0, 30,  0 | 27,  0, 24,  0)
        _mm256_packus_epi16(
            // mskLL:  (93,  0,  0,  0 | 90,  0,  0,  0 | 87,  0,  0,  0 | 84,  0,  0,  0 || 45,  0,  0,  0 | 42,  0,  0,  0 | 39,  0,  0,  0 | 36,  0,  0,  0)
            _mm256_and_si256(
                // gtrLL:  (93, 31, 30, 29 | 90, 93, 31, 30 | 87, 90, 93, 31 | 84, 87, 90, 93 || 45, 48, 51, 54 | 42, 45, 48, 51 | 39, 42, 45, 48 | 36, 39, 42, 45)
                _mm256_i32gather_epi32(
                    (const int *) memAddr,
                    // lowlow: (31,  0,  0,  0 | 30,  0,  0,  0 | 29,  0,  0,  0 | 28,  0,  0,  0 || 15,  0,  0,  0 | 14,  0,  0,  0 | 13,  0,  0,  0 | 12,  0,  0,  0)
                    _mm256_shuffle_epi8(idx, __m256i {
                        ~0xFF000000FFl | 0x0100000000, ~0xFF000000FFl | 0x0300000002,
                        ~0xFF000000FFl | 0x0100000000, ~0xFF000000FFl | 0x0300000002
                    }),
                    SCALE
                ),
                _mm256_set1_epi32(0xFF)
            ),
            // mskLH:  (81,  0,  0,  0 | 78,  0,  0,  0 | 75,  0,  0,  0 | 72,  0,  0,  0 || 33,  0,  0,  0 | 30,  0,  0,  0 | 27,  0,  0,  0 | 24,  0,  0,  0)
            _mm256_and_si256(
                // gtrLH:  (81, 84, 87, 90 | 78, 81, 84, 87 | 75, 78, 81, 84 | 72, 75, 78, 81 || 33, 36, 39, 42 | 30, 33, 36, 39 | 27, 30, 33, 36 | 24, 27, 30, 33)
                _mm256_i32gather_epi32(
                    (const int *) memAddr,
                    // lowhig: (27,  0,  0,  0 | 26,  0,  0,  0 | 25,  0,  0,  0 | 24,  0,  0,  0 || 11,  0,  0,  0 | 10,  0,  0,  0 |  9,  0,  0,  0 |  8,  0,  0,  0)
                    _mm256_shuffle_epi8(idx, __m256i {
                        ~0xFF000000FFl | 0x0500000004, ~0xFF000000FFl | 0x0700000006,
                        ~0xFF000000FFl | 0x0500000004, ~0xFF000000FFl | 0x0700000006
                    }),
                    SCALE
                ),
                _mm256_set1_epi32(0xFF)
            )
        ),
        // pckHih: (69,  0, 66,  0 | 63,  0, 60,  0 | 57,  0, 54,  0 | 51,  0, 48,  0 || 21,  0, 18,  0 | 15,  0, 12,  0 |  9,  0,  6,  0 |  3,  0,  0,  0)
        _mm256_packus_epi16(
            // mskHL:  (69,  0,  0,  0 | 66,  0,  0,  0 | 63,  0,  0,  0 | 60,  0,  0,  0 || 21,  0,  0,  0 | 18,  0,  0,  0 | 15,  0,  0,  0 | 12,  0,  0,  0)
            _mm256_and_si256(
                // gtrHL:  (69, 72, 75, 78 | 66, 69, 72, 75 | 63, 66, 69, 72 | 60, 63, 66, 69 || 21, 24, 27, 30 | 18, 21, 24, 27 | 15, 18, 21, 24 | 12, 15, 18, 21)
                _mm256_i32gather_epi32(
                    (const int *) memAddr,
                    // higlow: (23,  0,  0,  0 | 22,  0,  0,  0 | 21,  0,  0,  0 | 20,  0,  0,  0 ||  7,  0,  0,  0 |  6,  0,  0,  0 |  5,  0,  0,  0 |  4,  0,  0,  0)
                    _mm256_shuffle_epi8(idx, __m256i {
                        ~0xFF000000FFl | 0x0900000008, ~0xFF000000FFl | 0x0B0000000A,
                        ~0xFF000000FFl | 0x0900000008, ~0xFF000000FFl | 0x0B0000000A
                    }),
                    SCALE
                ),
                _mm256_set1_epi32(0xFF)
            ),
            // mskHH:  (57,  0,  0,  0 | 54,  0,  0,  0 | 51,  0,  0,  0 | 48,  0,  0,  0 ||  9,  0,  0,  0 |  6,  0,  0,  0 |  3,  0,  0,  0 |  0,  0,  0,  0)
            _mm256_and_si256(
                // gtrHH:  (57, 60, 63, 66 | 54, 57, 60, 63 | 51, 54, 57, 60 | 48, 51, 54, 57 ||  9, 12, 15, 18 |  6,  9, 12, 15 |  3,  6,  9, 12 |  0,  3,  6,  9)
                _mm256_i32gather_epi32(
                    (const int *) memAddr,
                    // highig: (19,  0,  0,  0 | 18,  0,  0,  0 | 17,  0,  0,  0 | 16,  0,  0,  0 ||  3,  0,  0,  0 |  2,  0,  0,  0 |  1,  0,  0,  0 |  0,  0,  0,  0)
                    _mm256_shuffle_epi8(idx, __m256i {
                        ~0xFF000000FFl | 0x0D0000000C, ~0xFF000000FFl | 0x0F0000000E,
                        ~0xFF000000FFl | 0x0D0000000C, ~0xFF000000FFl | 0x0F0000000E
                    }),
                    SCALE
                ),
                _mm256_set1_epi32(0xFF)
            )
        )
    );
}

template <typename TValue, typename TSize, TSize SCALE>
inline __m256i
seqan_mm256_i16gather_epi(TValue const * memAddr,
                          __m256i const & idx,
                          std::integral_constant<TSize, SCALE> const & /*scale*/)
{
    using TUnsignedValue = typename MakeUnsigned<TValue>::Type;

    // The cast makes sure that the max value of TValue = (u)int64_t and
    // (u)int32_t will be max value of int16_t (i.e. `~0` in int16_t), because
    // the resulting __m256i can only hold int16_t values.
    //
    // NOTE(marehr): the masking is only needed for TValue = (u)int8_t and
    // (u)int16_t. It could be omitted if _mm256_packus_epi32 would be exchanged
    // by _mm256_packs_epi32, because for (u)int32_t and (u)int64_t the masking
    // operations are basically the identity function.
    constexpr int const mask = static_cast<uint16_t>(MaxValue<TUnsignedValue>::VALUE);

    // 1. Unpack low idx values and interleave with 0 and gather from memAddr.
    // 2. Unpack high idx values and interleave with 0, than gather from memAddr.
    // 3. Merge 2 8x32 vectors into 1x16 vector by signed saturation. This operation reverts the interleave by the unpack operations above.
    //
    // The following is an example for SimdVector<uint16_t, 16> idx and uint16_t
    // const * memAddr:
    // mem:    ( 0,  0,  3,  0 |  6,  0,  9,  0 | 12,  0, 15,  0 | 18,  0, 21,  0 || 24,  0, 27,  0 | 30,  0, 33,  0 | 36,  0, 39,  0 | 42,  0, 45,  0)
    // idx:    (15,  0, 14,  0 | 13,  0, 12,  0 | 11,  0, 10,  0 |  9,  0,  8,  0 ||  7,  0,  6,  0 |  5,  0,  4,  0 |  3,  0,  2,  0 |  1,  0,  0,  0)
    // pack:   (45,  0, 42,  0 | 39,  0, 36,  0 | 33,  0, 30,  0 | 27,  0, 24,  0 || 21,  0, 18,  0 | 15,  0, 12,  0 |  9,  0,  6,  0 |  3,  0,  0,  0)
    return _mm256_packus_epi32(
        // mskLow: (45,  0,  0,  0 | 42,  0,  0,  0 | 39,  0,  0,  0 | 36,  0,  0,  0 || 21,  0,  0,  0 | 18,  0,  0,  0 | 15,  0,  0,  0 | 12,  0,  0,  0)
        _mm256_and_si256(
            // gtrLow: (45,  0, 15,  0 | 42,  0, 45,  0 | 39,  0, 42,  0 | 36,  0, 39,  0 || 21,  0, 24,  0 | 18,  0, 21,  0 | 15,  0, 18,  0 | 12,  0, 15,  0)
            _mm256_i32gather_epi32(
                (const int *) memAddr,
                // low:    (15,  0,  0,  0 | 14,  0,  0,  0 | 13,  0,  0,  0 | 12,  0,  0,  0 ||  7,  0,  0,  0 |  6,  0,  0,  0 |  5,  0,  0,  0 |  4,  0,  0,  0)
                _mm256_unpacklo_epi16(
                    idx, _mm256_set1_epi16(0)
                ),
                SCALE
            ),
            _mm256_set1_epi32(mask)
        ),
        // mskHih: (33,  0,  0,  0 | 30,  0,  0,  0 | 27,  0,  0,  0 | 24,  0,  0,  0 ||  9,  0,  0,  0 |  6,  0,  0,  0 |  3,  0,  0,  0 |  0,  0,  0,  0)
        _mm256_and_si256(
            // gtrHih: (33,  0, 36,  0 | 30,  0, 33,  0 | 27,  0, 30,  0 | 24,  0, 27,  0 ||  9,  0, 12,  0 |  6,  0,  9,  0 |  3,  0,  6,  0 |  0,  0,  3,  0)
            _mm256_i32gather_epi32(
                (const int *) memAddr,
                // high:   (11,  0,  0,  0 | 10,  0,  0,  0 |  9,  0,  0,  0 |  8,  0,  0,  0 ||  3,  0,  0,  0 |  2,  0,  0,  0 |  1,  0,  0,  0 |  0,  0,  0,  0)
                _mm256_unpackhi_epi16(
                    idx, _mm256_set1_epi16(0)
                ),
                SCALE
            ),
            _mm256_set1_epi32(mask)
        )
    );
}

template <typename TValue, typename TSize, TSize SCALE>
inline __m256i
seqan_mm256_i32gather_epi(TValue const * memAddr,
                          __m256i const & idx,
                          std::integral_constant<TSize, SCALE> const & /*scale*/)
{
    using TUnsignedValue = typename MakeUnsigned<TValue>::Type;
    constexpr auto const mask = static_cast<uint32_t>(MaxValue<TUnsignedValue>::VALUE);

    return _mm256_and_si256(
        _mm256_i32gather_epi32((const int *) memAddr, idx, SCALE),
        _mm256_set1_epi32(mask)
    );
}

template <typename TValue, typename TSize, TSize SCALE>
inline __m256i
seqan_mm256_i64gather_epi(TValue const * memAddr,
                          __m256i const & idx,
                          std::integral_constant<TSize, SCALE> const & /*scale*/)
{
    using TUnsignedValue = typename MakeUnsigned<TValue>::Type;
    constexpr auto const mask = static_cast<uint64_t>(MaxValue<TUnsignedValue>::VALUE);

    return _mm256_and_si256(
        _mm256_i64gather_epi64((const long long *) memAddr, idx, SCALE),
        _mm256_set1_epi64x(mask)
    );
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
        seqan_mm256_i8gather_epi(
            memAddr,
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            scale
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int8_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 16>)
{
    // Note that memAddr is a signed integer type, thus a cast would extend the
    // sign. E.g., -3 = 253 in 8 bit, but would be 65533 in 16 bit.
    // Use _gather(uint8_t) and extend the sign to [u]int16_t.
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i16sign_extend_epis8(
            seqan_mm256_i16gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i16gather_epi(
            memAddr,
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            scale
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int8_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 8>)
{
    // Note that memAddr is a signed integer type, thus a cast would extend the
    // sign.
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i32sign_extend_epis8(
            seqan_mm256_i32gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int16_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 8>)
{
    // Note that memAddr is a signed integer type, thus a cast would extend the
    // sign.
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i32sign_extend_epis16(
            seqan_mm256_i32gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i32gather_epi(
            memAddr,
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            scale
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int8_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i64sign_extend_epis8(
            seqan_mm256_i64gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int16_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i64sign_extend_epis16(
            seqan_mm256_i64gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(int32_t const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i64sign_extend_epis32(
            seqan_mm256_i64gather_epi(
                memAddr,
                SEQAN_VECTOR_CAST_(__m256i const &, idx),
                scale
            )
        )
    );
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & scale,
        SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector,
        seqan_mm256_i64gather_epi(
            memAddr,
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            scale
        )
    );
}

// --------------------------------------------------------------------------
// _shuffleVector (256bit)
// --------------------------------------------------------------------------

inline __m256i
seqan_m256_shuffle_epi8(__m256i const & vector, __m256i const & indices)
{
    return _mm256_xor_si256(
        // shuffle bytes from the lower bytes of vector
        _mm256_shuffle_epi8(
            // repeat twice the low bytes of vector in a new __m256i vector i.e.
            //   vh[127:0] = v[127:0]
            //   vh[255:128] = v[127:0]
            _mm256_broadcastsi128_si256(
                 _mm256_extracti128_si256(vector, 0)
            ),
            // ((indices[i] << 3) & 0b1000 0000) ^ indices[i]:
            //   Adds the 5th bit of indices[i] as most significant bit. If the
            //   5th bit is set, that means that indices[i] >= 16.
            //   r = _mm256_shuffle_epi8(vl, indices) will set r[i] = 0 if the
            //   most significant bit of indices[i] is 1. Since this bit is the
            //   5th bit, r[i] = 0 if indices[i] >= 16 and r[i] = vl[indices[i]]
            //   if indices[i] < 16.
            _mm256_xor_si256(
                _mm256_and_si256(
                    _mm256_slli_epi16(indices, 3),
                    _mm256_set1_epi8(-127) // 0b1000 0000
                ),
                indices
            )
        ),
        // shuffle bytes from the higher bytes of vector
        _mm256_shuffle_epi8(
            // repeat twice the higher bytes of vector in a new __m256i vector
            // i.e.
            //   vh[127:0] = v[255:128]
            //   vh[255:128] = v[255:128]
            _mm256_broadcastsi128_si256(
                 _mm256_extracti128_si256(vector, 1)
            ),
            // indices[i] - 16:
            //   r = _mm256_shuffle_epi8(vh, indices)
            //   will return r[i] = 0 if the most significant bit of the byte
            //   indices[i] is 1. Thus, indices[i] - 16 will select all high
            //   bytes in vh, i.e. r[i] = vh[indices[i] - 16], if indices[i] >=
            //   16 and r[i] = 0 if indices[i] < 16.
            _mm256_sub_epi8(
                indices,
                _mm256_set1_epi8(16)
            )
        )
    );
}

inline __m256i
seqan_m256_shuffle_epi16(const __m256i a, const __m256i b)
{
    // multiply by 2
    __m256i idx = _mm256_slli_epi16(
        _mm256_permute4x64_epi64(b, 0b01010000),
        1
    );
    // _print(_mm256_add_epi8(idx, _mm256_set1_epi8(1)));
    // _print(        _mm256_unpacklo_epi8(
    //             idx,
    //             _mm256_add_epi8(idx, _mm256_set1_epi8(1))
    //         ));
    return seqan_m256_shuffle_epi8(
        a,
        // interleave idx[15:0]   = 2*indices[15],   ..., 2*indices[0]
        // with       idx[15:0]+1 = 2*indices[15]+1, ..., 2*indices[0]+1
        // => 2*indices[15]+1, 2*indices[15], ..., 2*indices[0]+1, 2*indices[0]
        _mm256_unpacklo_epi8(
            idx,
            _mm256_add_epi8(idx, _mm256_set1_epi8(1))
        )
    );
}

inline __m256i
seqan_m256_shuffle_epi32(const __m256i a, const __m256i b)
{
    // multiply by 4
    __m256i idx = _mm256_slli_epi16(
        _mm256_permutevar8x32_epi32(b, __m256i {0x0, 0x0, 0x1, 0x0}),
        2
    );
    return seqan_m256_shuffle_epi8(
        a,
        // interleave 4*indices[7]+1, 4*indices[7]+0; ..., 4*indices[0]+1, 4*indices[0]+0
        // with       4*indices[7]+3, 4*indices[7]+2; ..., 4*indices[0]+3, 4*indices[0]+2
        // => 4*indices[7]+3, 4*indices[7]+2; 4*indices[7]+1, 4*indices[7]+0;
        //    ...
        //    4*indices[0]+3, 4*indices[0]+2; 4*indices[0]+1, 4*indices[0]+0
        _mm256_unpacklo_epi16(
            // interleave idx[7:0]+0 = 4*indices[7]+0; ...; 4*indices[0]+0
            // with       idx[7:0]+1 = 4*indices[7]+1; ...; 4*indices[0]+1
            // => 4*indices[7]+1; 4*indices[7]+0; ...; 4*indices[0]+1; 4*indices[0]+0
            _mm256_unpacklo_epi8(
                idx,
                _mm256_add_epi8(idx, _mm256_set1_epi8(1))
            ),
            // interleave idx[7:0]+2 = 4*indices[7]+2; ...; 4*indices[0]+2
            // with       idx[7:0]+3 = 4*indices[7]+3; ...; 4*indices[0]+3
            // => 4*indices[7]+3; 4*indices[7]+2; ...; 4*indices[0]+3; 4*indices[0]+2
            _mm256_unpacklo_epi8(
                _mm256_add_epi8(idx, _mm256_set1_epi8(2)),
                _mm256_add_epi8(idx, _mm256_set1_epi8(3))
            )
    ));
}

#define seqan_mm256_set_m128i(v0, v1) _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)

inline __m256i
seqan_m256_shuffle_epi64(const __m256i a, const __m256i b)
{
    __m128i lowidx = _mm256_extracti128_si256(
        // multiply by 8
        _mm256_slli_epi16(b, 3),
        0
    );

    __m256i idx = seqan_mm256_set_m128i(
        _mm_srli_si128(lowidx, 2),
        lowidx
    );

    return seqan_m256_shuffle_epi8(
        a,
        _mm256_unpacklo_epi32(
            // interleave 8*indices[3]+1, 8*indices[3]+0; ..., 8*indices[0]+1, 8*indices[0]+0
            // with       8*indices[3]+3, 8*indices[3]+2; ..., 8*indices[0]+3, 8*indices[0]+2
            // => 8*indices[3]+3, 8*indices[3]+2; 8*indices[3]+1, 8*indices[3]+0;
            //    ...
            //    8*indices[0]+3, 8*indices[0]+2; 8*indices[0]+1, 8*indices[0]+0
            _mm256_unpacklo_epi16(
                // interleave idx[3:0]+0 = 8*indices[3]+0; ...; 8*indices[0]+0
                // with       idx[3:0]+1 = 8*indices[3]+1; ...; 8*indices[0]+1
                // => 8*indices[3]+1; 8*indices[3]+0; ...; 8*indices[0]+1; 8*indices[0]+0
               _mm256_unpacklo_epi8(
                   idx,
                   _mm256_add_epi8(idx, _mm256_set1_epi8(1))
               ),
               // interleave idx[3:0]+2 = 8*indices[3]+2; ...; 8*indices[0]+2
               // with       idx[3:0]+3 = 8*indices[3]+3; ...; 8*indices[0]+3
               // => 8*indices[3]+3; 8*indices[3]+2; ...; 8*indices[0]+3; 8*indices[0]+2
               _mm256_unpacklo_epi8(
                   _mm256_add_epi8(idx, _mm256_set1_epi8(2)),
                   _mm256_add_epi8(idx, _mm256_set1_epi8(3))
               )
           ),
           // interleave 8*indices[3]+5, 8*indices[3]+4; ..., 8*indices[0]+5, 8*indices[0]+4
           // with       8*indices[3]+7, 8*indices[3]+6; ..., 8*indices[0]+7, 8*indices[0]+6
           // => 8*indices[3]+7, 8*indices[3]+6; 8*indices[3]+5, 8*indices[3]+4;
           //    ...
           //    8*indices[0]+7, 8*indices[0]+6; 8*indices[0]+5, 8*indices[0]+4
            _mm256_unpacklo_epi16(
                // interleave idx[3:0]+4 = 8*indices[3]+4; ...; 8*indices[0]+4
                // with       idx[3:0]+5 = 8*indices[3]+5; ...; 8*indices[0]+5
                // => 8*indices[3]+5; 8*indices[3]+4; ...; 8*indices[0]+5; 8*indices[0]+4
                _mm256_unpacklo_epi8(
                    _mm256_add_epi8(idx, _mm256_set1_epi8(4)),
                    _mm256_add_epi8(idx, _mm256_set1_epi8(5))
                ),
                // interleave idx[3:0]+6 = 8*indices[3]+6; ...; 8*indices[0]+6
                // with       idx[3:0]+7 = 8*indices[3]+7; ...; 8*indices[0]+7
                // => 8*indices[3]+7; 8*indices[3]+6; ...; 8*indices[0]+7; 8*indices[0]+6
                _mm256_unpacklo_epi8(
                    _mm256_add_epi8(idx, _mm256_set1_epi8(6)),
                    _mm256_add_epi8(idx, _mm256_set1_epi8(7))
                )
            )
        )
    );
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<32, 16>, SimdParams_<16, 16>)
{
    // copy 2nd 64bit word to 3rd, compute 2*idx
    __m256i idx = _mm256_slli_epi16(_mm256_permute4x64_epi64(_mm256_castsi128_si256(SEQAN_VECTOR_CAST_(const __m128i &, indices)), 0x50), 1);

    // interleave with 2*idx+1 and call shuffle
    return SEQAN_VECTOR_CAST_(TSimdVector1,
        _mm256_shuffle_epi8(
            SEQAN_VECTOR_CAST_(const __m256i &, vector),
            _mm256_unpacklo_epi8(
                idx,
                _mm256_add_epi8(
                    idx, _mm256_set1_epi8(1)
                )
            )
        )
    );
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<32, 32>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi8(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<32, 16>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi16(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<32, 8>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi32(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<32, 4>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi64(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

// --------------------------------------------------------------------------
// _transposeMatrix (256bit)
// --------------------------------------------------------------------------

// emulate missing _mm256_unpacklo_epi128/_mm256_unpackhi_epi128 instructions
inline __m256i _mm256_unpacklo_epi128(__m256i const & a, __m256i const & b)
{
    return _mm256_permute2x128_si256(a, b, 0x20);
//    return _mm256_inserti128_si256(a, _mm256_extracti128_si256(b, 0), 1);
}

inline __m256i _mm256_unpackhi_epi128(__m256i const & a, __m256i const & b)
{
    return _mm256_permute2x128_si256(a, b, 0x31);
//    return _mm256_inserti128_si256(b, _mm256_extracti128_si256(a, 1), 0);
}

template <typename TSimdVector>
inline void
_transposeMatrix(TSimdVector matrix[], SimdMatrixParams_<32, 32, 8>)
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    static const unsigned char bitRev[] = { 0, 8, 4,12, 2,10, 6,14, 1, 9, 5,13, 3,11, 7,15,
                                           16,24,20,28,18,26,22,30,17,25,21,29,19,27,23,31};

    // transpose a 32x32 byte matrix
    __m256i tmp1[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i]    = _mm256_unpacklo_epi8(
            SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i]),
            SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i+1])
        );
        tmp1[i+16] = _mm256_unpackhi_epi8(
            SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i]),
            SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i+1])
        );
    }
    __m256i  tmp2[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i]    = _mm256_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+16] = _mm256_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i]    = _mm256_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+16] = _mm256_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i]    = _mm256_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+16] = _mm256_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        matrix[bitRev[i]]    = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_unpacklo_epi128(tmp2[2*i],tmp2[2*i+1]));
        matrix[bitRev[i+16]] = SEQAN_VECTOR_CAST_(TSimdVector, _mm256_unpackhi_epi128(tmp2[2*i],tmp2[2*i+1]));
    }
}

// --------------------------------------------------------------------------
// Function _testAllZeros (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, int)
inline _testAllZeros(TSimdVector const & vector, TSimdVector const & mask, SimdParams_<32>)
{
    return _mm256_testz_si256(SEQAN_VECTOR_CAST_(const __m256i &, vector),
                              SEQAN_VECTOR_CAST_(const __m256i &, mask));
}

// --------------------------------------------------------------------------
// Function _testAllOnes (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline int _testAllOnes(TSimdVector const & vector, SimdParams_<32>)
{
    __m256i vec = SEQAN_VECTOR_CAST_(const __m256i &, vector);
    return _mm256_testc_si256(vec, _mm256_cmpeq_epi32(vec, vec));
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX2_H_
