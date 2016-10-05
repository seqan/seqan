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
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Float,     float,          32)

// SimdParams_<32, 4>: 256bit = 4 elements * 8 * 8bit
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int64,     int64_t,        32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt64,    uint64_t,       32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Double,    double,         32)

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
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & x, std::index_sequence<0> const &, SimdParams_<32, 32>) { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_epi8(std::get<0>(x)); }
template <typename TSimdVector, typename ...TValue>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & x, std::index_sequence<0> const &, SimdParams_<32, 16>) { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_epi16(std::get<0>(x)); }
template <typename TSimdVector, typename ...TValue>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & x, std::index_sequence<0> const &, SimdParams_<32, 8>)  { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_epi32(std::get<0>(x)); }
template <typename TSimdVector, typename ...TValue>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & x, std::index_sequence<0> const &, SimdParams_<32, 4>)  { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_epi64x(std::get<0>(x)); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, std::tuple<float> const & x,  std::index_sequence<0> const &, SimdParams_<32, 8>)  { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_ps(std::get<0>(x)); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, std::tuple<double> const & x, std::index_sequence<0> const &, SimdParams_<32, 4>)  { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set1_pd(std::get<0>(x)); }

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 32>)
{
    // reverse argument list 0, 1, 2, 3 -> 3, 2, 1, 0
    SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set_epi8(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...);
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 16>)
{
    // reverse argument list 0, 1, 2, 3 -> 3, 2, 1, 0
    SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set_epi16(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...);
}
template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 8>)
{
    // reverse argument list 0, 1, 2, 3 -> 3, 2, 1, 0
    SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set_epi32(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...);
}

template <typename TSimdVector, typename ...TValue, size_t ...INDICES>
inline void _fillVector(TSimdVector &vector, std::tuple<TValue...> const & args, std::index_sequence<INDICES...> const &, SimdParams_<32, 4>)
{
    // reverse argument list 0, 1, 2, 3 -> 3, 2, 1, 0
    SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_set_epi64x(std::get<sizeof...(INDICES) - 1 - INDICES>(args)...);
}

// --------------------------------------------------------------------------
// _clearVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, L>) { SEQAN_VECTOR_CAST_LVALUE_(__m256i&, vector) = _mm256_setzero_si256(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, 8>) { SEQAN_VECTOR_CAST_LVALUE_(__m256&, vector) = _mm256_setzero_ps(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, 4>) { SEQAN_VECTOR_CAST_LVALUE_(__m256d&, vector) = _mm256_setzero_pd(); }

// --------------------------------------------------------------------------
// _createVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<32, 32>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi8(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<32, 16>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi16(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<32, 8>)  { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi32(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_< 32, 4>)  { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_epi64x(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(float x,  SimdParams_<32, 8>)  { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_ps(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(double x, SimdParams_<32, 4>)  { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_set1_pd(x)); }

// --------------------------------------------------------------------------
// _cmpEq (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                             SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _cmpGt (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                             SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_cmpgt_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseOr (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseOr(TSimdVector &a, TSimdVector &b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_or_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                           SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAnd (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_and_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                            SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_and_ps(SEQAN_VECTOR_CAST_(const __m256&, a),
                                                         SEQAN_VECTOR_CAST_(const __m256&, b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_and_pd(SEQAN_VECTOR_CAST_(const __m256d&, a),
                                                         SEQAN_VECTOR_CAST_(const __m256d&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, int L>
inline TSimdVector _bitwiseAndNot(TSimdVector &a, TSimdVector &b, SimdParams_<32, L>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_andnot_si256(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _bitwiseNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));

}
template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_cmpeq_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a), _mm256_setzero_si256()));
}

// --------------------------------------------------------------------------
// _divide (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<32, 32>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi8(a, _mm256_set1_epi8(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<32, 16>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi16(a, _mm256_set1_epi16(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<32, 8>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi32(a, _mm256_set1_epi32(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<32, 4>) { return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_div_epi64(a, _mm256_set1_epi64x(b))); }

// --------------------------------------------------------------------------
// _add (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_add_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _sub (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_sub_epi64(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

// --------------------------------------------------------------------------
// _mult (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &/*b*/, SimdParams_<32, 32>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("AVX2 intrinsics for multiplying 8 bit values not implemented!");
    return a;
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_mullo_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                 SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_mullo_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                                 SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &/*b*/, SimdParams_<32, 4>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("AVX2 intrinsics for multiplying 64 bit values not implemented!");
    return a;
}

// --------------------------------------------------------------------------
// _max (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &b, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi8(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                              SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &b, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi16(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &b, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector,
                              _mm256_max_epi32(SEQAN_VECTOR_CAST_(const __m256i&, a),
                                               SEQAN_VECTOR_CAST_(const __m256i&, b)));
}

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &/*b*/, SimdParams_<32, 4>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("AVX2 intrinsics for max on 64 bit values not implemented!");
    return a;
}

// --------------------------------------------------------------------------
// _blend (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, int L>
inline TSimdVector _blend(TSimdVector const &a, TSimdVector const &b, TSimdVectorMask const &mask, SimdParams_<32, L>)
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
inline void _storeu(T * memAddr, TSimdVector &vec, SimdParams_<32, L>)
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
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi16(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm) & _mm256_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 16>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi16(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 8>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi32(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 4>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector, _mm256_srli_epi64(SEQAN_VECTOR_CAST_(const __m256i &, vector), imm));
}

// --------------------------------------------------------------------------
// _gather (256bit)
// --------------------------------------------------------------------------

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector _gather(TValue const * /*memAddr*/,
                           TSimdVector const & /*idx*/,
                           std::integral_constant<TSize, SCALE> const & /*scale*/,
                           SimdParams_<32, 32>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("SSE intrinsics for gather 8 bit values not implemented!");
    return createVector<TSimdVector>(0);
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector _gather(TValue const * memAddr,
                           TSimdVector const & idx,
                           std::integral_constant<TSize, SCALE> const & /*scale*/,
                           SimdParams_<32, 16>)
{
    // 1. Unpack low idx values and interleave with 0 and gather from memAddr.
    // 2. Unpack high idx values and interleave with 0, than gather from memAddr.
    // 3. Merge 2 8x32 vectors into 1x16 vector by signed saturation. This operation reverts the interleave by the unpack operations above.
    return SEQAN_VECTOR_CAST_(
      TSimdVector,
      _mm256_packs_epi32(
        _mm256_i32gather_epi32(
          static_cast<int32_t const *>(memAddr),
          _mm256_unpacklo_epi16(
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            _mm256_set1_epi16(0)
          ),
          SCALE
        ),
        _mm256_i32gather_epi32(
          static_cast<int32_t const *>(memAddr),
          _mm256_unpackhi_epi16(
            SEQAN_VECTOR_CAST_(__m256i const &, idx),
            _mm256_set1_epi16(0)
          ),
          SCALE
        )
      )
    );
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector _gather(TValue const * /*memAddr*/,
                           TSimdVector const & /*idx*/,
                           std::integral_constant<TSize, SCALE> const & /*scale*/,
                           SimdParams_<32, 8>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("SSE intrinsics for gather 32 bit values not implemented!");
    return createVector<TSimdVector>(0);
}

template <typename TValue, typename TSimdVector, typename TSize, TSize SCALE>
inline TSimdVector _gather(TValue const * /*memAddr*/,
                           TSimdVector const & /*idx*/,
                           std::integral_constant<TSize, SCALE> const & /*scale*/,
                           SimdParams_<32, 4>)
{
    SEQAN_SKIP_TEST;
    SEQAN_ASSERT_FAIL("SSE intrinsics for gather 64 bit values not implemented!");
    return createVector<TSimdVector>(0);
}

// --------------------------------------------------------------------------
// _shuffleVector (256bit)
// --------------------------------------------------------------------------

inline __m256i
seqan_m256_shuffle_epi8(__m256i const &vector, __m256i const &indices)
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
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 16>, SimdParams_<16, 16>)
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
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 32>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi8(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 16>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi16(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 8>, SimdParams_<32, 32>)
{
    return SEQAN_VECTOR_CAST_(TSimdVector1, seqan_m256_shuffle_epi32(
        SEQAN_VECTOR_CAST_(const __m256i &, vector),
        SEQAN_VECTOR_CAST_(const __m256i &, indices)
    ));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 4>, SimdParams_<32, 32>)
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
inline __m256i _mm256_unpacklo_epi128(__m256i const &a, __m256i const &b)
{
    return _mm256_permute2x128_si256(a, b, 0x20);
//    return _mm256_inserti128_si256(a, _mm256_extracti128_si256(b, 0), 1);
}

inline __m256i _mm256_unpackhi_epi128(__m256i const &a, __m256i const &b)
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
        tmp1[i]    = _mm256_unpacklo_epi8(SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i+1]));
        tmp1[i+16] = _mm256_unpackhi_epi8(SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i]), SEQAN_VECTOR_CAST_(const __m256i &, matrix[2*i+1]));
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
inline _testAllZeros(TSimdVector const &vector, TSimdVector const &mask, SimdParams_<32>)
{
    return _mm256_testz_si256(SEQAN_VECTOR_CAST_(const __m256i &, vector),
                              SEQAN_VECTOR_CAST_(const __m256i &, mask));
}

// --------------------------------------------------------------------------
// Function _testAllOnes (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline int _testAllOnes(TSimdVector const &vector, SimdParams_<32>)
{
    __m256i vec = SEQAN_VECTOR_CAST_(const __m256i &, vector);
    return _mm256_testc_si256(vec, _mm256_cmpeq_epi32(vec, vec));
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_AVX2_H_
