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

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_SSE4_2_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_SSE4_2_H_

namespace seqan {

// SimdParams_<8, 8>: 64bit = 8 elements * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Char,      char,           8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8SChar,     signed char,    8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UChar,     unsigned char,  8)

// SimdParams_<8, 4>: 64bit = 4 elements * 2 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Short,     short,          8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UShort,    unsigned short, 8)

// SimdParams_<8, 2>: 64bit = 2 elements * 4 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int,       int,            8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt,      unsigned int,   8)

// SimdParams_<16, 16>: 128bit = 16 elements * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Char,     char,           16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16SChar,    signed char,    16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UChar,    unsigned char,  16)

// SimdParams_<16, 8>: 128bit = 8 elements * 2 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Short,     short,          16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UShort,    unsigned short, 16)

// SimdParams_<16, 4>: 128bit = 4 elements * 4 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int,       int,            16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt,      unsigned int,   16)

// SimdParams_<16, 2>: 128bit = 2 elements * 8 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int64,     int64_t,        16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt64,    uint64_t,       16)

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// _fillVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename... TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &,
            SimdParams_<16, 16> const &)
{
  vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi8(std::get<0>(x)));
}

template <typename TSimdVector, typename... TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &,
            SimdParams_<16, 8> const &)
{
  vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi16(std::get<0>(x)));
}

template <typename TSimdVector, typename... TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &,
            SimdParams_<16, 4> const &)
{
  vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi32(std::get<0>(x)));
}

template <typename TSimdVector, typename... TValue>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & x,
            std::index_sequence<0> const &,
            SimdParams_<16, 2> const &)
{
  vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi64x(std::get<0>(x)));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args,
            std::index_sequence<INDICES...> const &,
            SimdParams_<16, 16> const &)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_setr_epi8(std::get<INDICES>(args)...));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args,
            std::index_sequence<INDICES...> const &,
            SimdParams_<16, 8> const &)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_setr_epi16(std::get<INDICES>(args)...));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args,
            std::index_sequence<INDICES...> const &,
            SimdParams_<16, 4> const &)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_setr_epi32(std::get<INDICES>(args)...));
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void
_fillVector(TSimdVector & vector,
            std::tuple<TValue...> const & args,
            std::index_sequence<INDICES...> const &,
            SimdParams_<16, 2> const &)
{
    // reverse argument list 0, 1 -> 1, 0
    // NOTE(marehr): Intel linux fails to reverse argument list and only
    // _mm_set_epi64x has no reverse equivalent
    // NOTE(rrahn): For g++-4.9 the set_epi function is a macro, which does not work with parameter pack expansion.
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_set_epi64x(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...));
}

// --------------------------------------------------------------------------
// _clearVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector & vector, SimdParams_<16, L>)
{
    vector = SEQAN_VECTOR_CAST_(TSimdVector, _mm_setzero_si128());
}

// --------------------------------------------------------------------------
// _createVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi8(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi16(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi32(x));
}

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue const x, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_set1_epi64x(x));
}

// --------------------------------------------------------------------------
// cmpEq (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                             SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _cmpGt (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                             SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, uint8_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x80 = ~0x7F (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi8(
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, a), _mm_set1_epi8(~0x7F)),
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, b), _mm_set1_epi8(~0x7F))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, uint16_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x8000 = ~0x7FFF (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi16(
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, a), _mm_set1_epi16(~0x7FFF)),
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, b), _mm_set1_epi16(~0x7FFF))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, uint32_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x80000000 = ~0x7FFFFFFF (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi32(
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, a), _mm_set1_epi32(~0x7FFFFFFF)),
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, b), _mm_set1_epi32(~0x7FFFFFFF))));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, int64_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, uint64_t>)
{
    // There is no unsigned cmpgt, we reduce it to the signed case.
    // Note that 0x8000000000000000ul = ~0x7FFFFFFFFFFFFFFFul (prevent overflow messages).
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpgt_epi64(
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, a) ,_mm_set1_epi64x(~0x7FFFFFFFFFFFFFFFul)),
                                  _mm_xor_si128(SEQAN_VECTOR_CAST_(const __m128i&, b), _mm_set1_epi64x(~0x7FFFFFFFFFFFFFFFul))));
}

// --------------------------------------------------------------------------
// _bitwiseOr (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseOr(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_or_si128(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAnd (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAnd(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_and_si128(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAndNot(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_andnot_si128(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                               SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseNot (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                             _mm_setzero_si128()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              _mm_setzero_si128()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              _mm_setzero_si128()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector const & a, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              _mm_setzero_si128()));
}

// --------------------------------------------------------------------------
// _divide (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_div_epi8(a, _mm_set1_epi8(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_div_epi16(a, _mm_set1_epi16(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_div_epi32(a, _mm_set1_epi32(b)));
}

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector const & a, int b, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_div_epi64(a, _mm_set1_epi64x(b)));
}

// --------------------------------------------------------------------------
// _add (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_add_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_add_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_add_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_add_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _sub (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_sub_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_sub_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_sub_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_sub_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

// --------------------------------------------------------------------------
// _mult (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const &/*b*/, SimdParams_<16, 16>)
{
    SEQAN_ASSERT_FAIL("SSE intrinsics for multiplying 8 bit values not implemented!");
    return a;
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_mullo_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_mullo_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector const & a, TSimdVector const &/*b*/, SimdParams_<16, 2>)
{
    SEQAN_ASSERT_FAIL("SSE intrinsics for multiplying 64 bit values not implemented!");
    return a;
}

// --------------------------------------------------------------------------
// _max (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, uint8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epu8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, uint16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epu16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, uint32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epu32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, int64_t>)
{
#if defined(__AVX512VL__)
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
#else // defined(__AVX512VL__)
    return blend(b, a, cmpGt(a, b));
#endif // defined(__AVX512VL__)
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, uint64_t>)
{
#if defined(__AVX512VL__)
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_max_epu64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
#else // defined(__AVX512VL__)
    return blend(b, a, cmpGt(a, b));
#endif // defined(__AVX512VL__)
}


// --------------------------------------------------------------------------
// _min (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, int8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 16, uint8_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epu8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                           SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, int16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epi16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 8, uint16_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epu16(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, int32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epi32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 4, uint32_t>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epu32(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, int64_t>)
{
#if defined(__AVX512VL__)
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epi64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
#else // defined(__AVX512VL__)
    return blend(a, b, cmpGt(a, b));
#endif // defined(__AVX512VL__)
}

template <typename TSimdVector>
inline TSimdVector _min(TSimdVector const & a, TSimdVector const & b, SimdParams_<16, 2, uint64_t>)
{
#if defined(__AVX512VL__)
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_min_epu64(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                            SEQAN_VECTOR_CAST_(const __m128i&, b)));
#else // defined(__AVX512VL__)
    return blend(a, b, cmpGt(a, b));
#endif // defined(__AVX512VL__)
}

// --------------------------------------------------------------------------
// _blend (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, int L>
inline TSimdVector _blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask, SimdParams_<16, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm_blendv_epi8(SEQAN_VECTOR_CAST_(const __m128i&, a),
                                              SEQAN_VECTOR_CAST_(const __m128i&, b),
                                              SEQAN_VECTOR_CAST_(const __m128i&, mask)));
}

// --------------------------------------------------------------------------
// _storeu (128bit)
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector, int L>
inline void _storeu(T * memAddr, TSimdVector const & vec, SimdParams_<16, L>)
{
    _mm_storeu_si128((__m128i*)memAddr, reinterpret_cast<const __m128i &>(vec));
}

// ----------------------------------------------------------------------------
// Function _load() 128bit
// ----------------------------------------------------------------------------

template <typename TSimdVector, typename T, int L>
inline TSimdVector _load(T const * memAddr, SimdParams_<16, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_load_si128((__m128i const *) memAddr));
}

// --------------------------------------------------------------------------
// _shiftRightLogical (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_srli_epi16(SEQAN_VECTOR_CAST_(const __m128i &, vector), imm) & _mm_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<16, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_srli_epi16(SEQAN_VECTOR_CAST_(const __m128i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<16, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_srli_epi32(SEQAN_VECTOR_CAST_(const __m128i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const & vector, const int imm, SimdParams_<16, 2>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm_srli_epi64(SEQAN_VECTOR_CAST_(const __m128i &, vector), imm));
}

// --------------------------------------------------------------------------
// _gather (128bit)
// --------------------------------------------------------------------------

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE, typename TSimdParams>
inline TSimdVector
_gather(TValue const * memAddr,
        TSimdVector const & idx,
        std::integral_constant<TSize, SCALE> const & /*scale*/,
        TSimdParams)
{
    TSimdVector ret;
    for (auto i = 0u; i < LENGTH<TSimdVector>::VALUE; ++i)
    {
        ret[i] = memAddr[idx[i]];
    }
    return ret;
}

// --------------------------------------------------------------------------
// _shuffleVector (128bit)
// --------------------------------------------------------------------------

inline __m128i
seqan_mm_shuffle_epi16(const __m128i a, const __m128i b)
{
    // multiply by 2
    __m128i idx = _mm_slli_epi16(b, 1);
    return _mm_shuffle_epi8(
        a,
        // interleave idx[7:0]   = 2*indices[7],   ..., 2*indices[0]
        // with       idx[7:0]+1 = 2*indices[7]+1, ..., 2*indices[0]+1
        // => 2*indices[7]+1, 2*indices[7], ..., 2*indices[0]+1, 2*indices[0]
        _mm_unpacklo_epi8(
            idx,
            _mm_add_epi8(idx, _mm_set1_epi8(1))
        )
    );
}

inline __m128i
seqan_mm_shuffle_epi32(const __m128i a, const __m128i b)
{
    // multiply by 4
    __m128i idx = _mm_slli_epi16(b, 2);
    return _mm_shuffle_epi8(
        a,
        // interleave 4*indices[3]+1, 4*indices[3]+0; ..., 4*indices[0]+1, 4*indices[0]+0
        // with       4*indices[3]+3, 4*indices[3]+2; ..., 4*indices[0]+3, 4*indices[0]+2
        // => 4*indices[3]+3, 4*indices[3]+2; 4*indices[3]+1, 4*indices[3]+0;
        //    ...
        //    4*indices[0]+3, 4*indices[0]+2; 4*indices[0]+1, 4*indices[0]+0
        _mm_unpacklo_epi16(
            // interleave idx[3:0]+0 = 4*indices[3]+0; ...; 4*indices[0]+0
            // with       idx[3:0]+1 = 4*indices[3]+1; ...; 4*indices[0]+1
            // => 4*indices[3]+1; 4*indices[3]+0; ...; 4*indices[0]+1; 4*indices[0]+0
            _mm_unpacklo_epi8(
                idx,
                _mm_add_epi8(idx, _mm_set1_epi8(1))
            ),
            // interleave idx[3:0]+2 = 4*indices[3]+2; ...; 4*indices[0]+2
            // with       idx[3:0]+3 = 4*indices[3]+3; ...; 4*indices[0]+3
            // => 4*indices[3]+3; 4*indices[3]+2; ...; 4*indices[0]+3; 4*indices[0]+2
            _mm_unpacklo_epi8(
                _mm_add_epi8(idx, _mm_set1_epi8(2)),
                _mm_add_epi8(idx, _mm_set1_epi8(3))
            )
    ));
}

inline __m128i
seqan_mm_shuffle_epi64(const __m128i a, const __m128i b)
{
    // multiply by 8
    __m128i idx = _mm_slli_epi16(b, 3);
    return _mm_shuffle_epi8(
        a,
        _mm_unpacklo_epi32(
            // interleave 8*indices[1]+1, 8*indices[1]+0; ..., 8*indices[0]+1, 8*indices[0]+0
            // with       8*indices[1]+3, 8*indices[1]+2; ..., 8*indices[0]+3, 8*indices[0]+2
            // => 8*indices[1]+3, 8*indices[1]+2; 8*indices[1]+1, 8*indices[1]+0;
            //    ...
            //    8*indices[0]+3, 8*indices[0]+2; 8*indices[0]+1, 8*indices[0]+0
            _mm_unpacklo_epi16(
                // interleave idx[1:0]+0 = 8*indices[1]+0; ...; 8*indices[0]+0
                // with       idx[1:0]+1 = 8*indices[1]+1; ...; 8*indices[0]+1
                // => 8*indices[1]+1; 8*indices[1]+0; ...; 8*indices[0]+1; 8*indices[0]+0
                _mm_unpacklo_epi8(
                    idx,
                    _mm_add_epi8(idx, _mm_set1_epi8(1))
                ),
                // interleave idx[1:0]+2 = 8*indices[1]+2; ...; 8*indices[0]+2
                // with       idx[1:0]+3 = 8*indices[1]+3; ...; 8*indices[0]+3
                // => 8*indices[1]+3; 8*indices[1]+2; ...; 8*indices[0]+3; 8*indices[0]+2
                _mm_unpacklo_epi8(
                    _mm_add_epi8(idx, _mm_set1_epi8(2)),
                    _mm_add_epi8(idx, _mm_set1_epi8(3))
                )
            ),
            // interleave 8*indices[1]+5, 8*indices[1]+4; ..., 8*indices[0]+5, 8*indices[0]+4
            // with       8*indices[1]+7, 8*indices[1]+6; ..., 8*indices[0]+7, 8*indices[0]+6
            // => 8*indices[1]+7, 8*indices[1]+6; 8*indices[1]+5, 8*indices[1]+4;
            //    ...
            //    8*indices[0]+7, 8*indices[0]+6; 8*indices[0]+5, 8*indices[0]+4
            _mm_unpacklo_epi16(
                // interleave idx[1:0]+4 = 8*indices[1]+4; ...; 8*indices[0]+4
                // with       idx[1:0]+5 = 8*indices[1]+5; ...; 8*indices[0]+5
                // => 8*indices[1]+5; 8*indices[1]+4; ...; 8*indices[0]+5; 8*indices[0]+4
                _mm_unpacklo_epi8(
                    _mm_add_epi8(idx, _mm_set1_epi8(4)),
                    _mm_add_epi8(idx, _mm_set1_epi8(5))
                ),
                // interleave idx[1:0]+6 = 8*indices[1]+6; ...; 8*indices[0]+6
                // with       idx[1:0]+7 = 8*indices[1]+7; ...; 8*indices[0]+7
                // => 8*indices[1]+7; 8*indices[1]+6; ...; 8*indices[0]+7; 8*indices[0]+6
                _mm_unpacklo_epi8(
                    _mm_add_epi8(idx, _mm_set1_epi8(6)),
                    _mm_add_epi8(idx, _mm_set1_epi8(7))
                )
            )
        )
    );
}

template <typename TSimdVector1, typename TSimdVector2>
[[deprecated("Here be dragons")]]
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<16, 8>, SimdParams_<8, 8>)
{
#if SEQAN_IS_32_BIT
    __m128i idx = _mm_slli_epi16(
        _mm_unpacklo_epi32(
            _mm_cvtsi32_si128(reinterpret_cast<const uint32_t &>(indices)),
            _mm_cvtsi32_si128(reinterpret_cast<const uint64_t &>(indices) >> 32)
        ),
        1
    );
#else
    __m128i idx = _mm_slli_epi16(_mm_cvtsi64_si128(reinterpret_cast<const uint64_t &>(indices)), 1);
#endif  // SEQAN_IS_32_BIT
    return SEQAN_VECTOR_CAST_(TSimdVector1,
        _mm_shuffle_epi8(
            SEQAN_VECTOR_CAST_(const __m128i &, vector),
            _mm_unpacklo_epi8(idx, _mm_add_epi8(idx, _mm_set1_epi8(1)))
        ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<16, 16>, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector1,
        _mm_shuffle_epi8(
            SEQAN_VECTOR_CAST_(const __m128i &, vector),
            SEQAN_VECTOR_CAST_(const __m128i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<16, 8>, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector1,
        seqan_mm_shuffle_epi16(
            SEQAN_VECTOR_CAST_(const __m128i &, vector),
            SEQAN_VECTOR_CAST_(const __m128i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<16, 4>, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector1,
        seqan_mm_shuffle_epi32(
            SEQAN_VECTOR_CAST_(const __m128i &, vector),
            SEQAN_VECTOR_CAST_(const __m128i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices, SimdParams_<16, 2>, SimdParams_<16, 16>)
{
    return SEQAN_VECTOR_CAST_(
        TSimdVector1,
        seqan_mm_shuffle_epi64(
            SEQAN_VECTOR_CAST_(const __m128i &, vector),
            SEQAN_VECTOR_CAST_(const __m128i &, indices)
    ));
}

// --------------------------------------------------------------------------
// _transposeMatrix (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline void
_transposeMatrix(TSimdVector matrix[], SimdMatrixParams_<8, 8, 8>)
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    static const unsigned char bitRev[] = {0,4,2,6,1,5,3,7};

    // transpose a 8x8 byte matrix
    __m64 tmp1[8];
    for (int i = 0; i < 4; ++i)
    {
        tmp1[i]   = _mm_unpacklo_pi8(SEQAN_VECTOR_CAST_(const __m64 &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m64 &, matrix[2*i+1]));
        tmp1[i+4] = _mm_unpackhi_pi8(SEQAN_VECTOR_CAST_(const __m64 &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m64 &, matrix[2*i+1]));
    }
    __m64 tmp2[8];
    for (int i = 0; i < 4; ++i)
    {
        tmp2[i]   = _mm_unpacklo_pi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+4] = _mm_unpackhi_pi16(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 4; ++i)
    {
        matrix[bitRev[i]]   = SEQAN_VECTOR_CAST_(TSimdVector, _mm_unpacklo_pi32(tmp2[2*i], tmp2[2*i+1]));
        matrix[bitRev[i+4]] = SEQAN_VECTOR_CAST_(TSimdVector, _mm_unpackhi_pi32(tmp2[2*i], tmp2[2*i+1]));
    }
}

template <typename TSimdVector>
inline void
_transposeMatrix(TSimdVector matrix[], SimdMatrixParams_<16, 16, 8>)
{
    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    static const unsigned char bitRev[] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};

    // transpose a 16x16 byte matrix
    //
    // matrix =
    // A0 A1 A2 ... Ae Af
    // B0 B1 B2 ... Be Bf
    // ...
    // P0 P1 P2 ... Pe Pf
    __m128i tmp1[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi8(SEQAN_VECTOR_CAST_(const __m128i &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m128i &, matrix[2*i+1]));
        tmp1[i+8] = _mm_unpackhi_epi8(SEQAN_VECTOR_CAST_(const __m128i &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m128i &, matrix[2*i+1]));
    }
    // tmp1[0]  = A0 B0 A1 B1 ... A7 B7
    // tmp1[1]  = C0 D0 C1 D1 ... C7 D7
    // ...
    // tmp1[7]  = O0 P0 O1 P1 ... O7 P7
    // tmp1[8]  = A8 B8 A9 B9 ... Af Bf
    // ...
    // tmp1[15] = O8 P8 O9 P9 ... Of Pf
    __m128i tmp2[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp2[i]   = _mm_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+8] = _mm_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    // tmp2[0]  = A0 B0 C0 D0 ... A3 B3 C3 D3
    // tmp2[1]  = E0 F0 G0 H0 ... E3 F3 G3 H3
    // ...
    // tmp2[3]  = M0 N0 O0 P0 ... M3 N3 O3 P3
    // tmp2[4]  = A8 B8 C8 D8 ... Ab Bb Cb Db
    // ...
    // tmp2[7]  = M8 N8 O8 P8 ... Mb Nb Ob Pb
    // tmp2[8]  = A4 B4 C4 D4 ... A7 B7 C7 D7
    // ..
    // tmp2[12] = Ac Bc Cc Dc ... Af Bf Cf Df
    // ...
    // tmp2[15] = Mc Nc Oc Pc ... Mf Nf Of Pf
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+8] = _mm_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    // tmp1[0]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    // ...
    // tmp1[4]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    for (int i = 0; i < 8; ++i)
    {
        matrix[bitRev[i]]   = SEQAN_VECTOR_CAST_(TSimdVector, _mm_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]));
        matrix[bitRev[i+8]] = SEQAN_VECTOR_CAST_(TSimdVector, _mm_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]));
    }
}

// --------------------------------------------------------------------------
// Function _testAllZeros (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, int)
inline _testAllZeros(TSimdVector const & vector, TSimdVector const & mask, SimdParams_<16>)
{
    return _mm_testz_si128(SEQAN_VECTOR_CAST_(const __m128i &, vector),
                           SEQAN_VECTOR_CAST_(const __m128i &, mask));
}

// --------------------------------------------------------------------------
// Function _testAllOnes (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline
SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, int)
_testAllOnes(TSimdVector const & vector, SimdParams_<16>)
{
    return _mm_test_all_ones(SEQAN_VECTOR_CAST_(const __m128i &, vector));
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_SSE4_2_H_
