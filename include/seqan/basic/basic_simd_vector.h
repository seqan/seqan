// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// ==========================================================================
// generic SIMD interface for SSE3 / AVX2
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_

#include <x86intrin.h>

#ifndef __SSE3__
 #warning "At least SSE3 instruction set must be available"
#endif

namespace seqan {

// ============================================================================
// Useful Macros
// ============================================================================

//this kind of access might not be supported by all compilers
#define SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)     \
template <typename TPosition>                            \
inline typename InnerValue<TSimdVector>::Type            \
value(TSimdVector &vector, TPosition pos)                \
{                                                        \
    return vector[pos];                                  \
}

#define SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)         \
template <typename TPosition, typename TValue2>                    \
inline void                                                        \
assignValue(TSimdVector &vector, TPosition pos, TValue2 value)     \
{                                                                  \
    vector[pos] = value;                                           \
}

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// define a concept and its models
// they allow us to define generic vector functions
SEQAN_CONCEPT(SimdVectorConcept, (T)) {};

// size of vector in byte
#if defined(__AVX2__)
 #define SEQAN_SIZEOF_MAX_VECTOR 32
#elif defined(__SSE3__)
 #define SEQAN_SIZEOF_MAX_VECTOR 16
#endif

// a metafunction returning the biggest supported SIMD vector
template <typename TValue, int LENGTH = SEQAN_SIZEOF_MAX_VECTOR / sizeof(TValue)>
struct SimdVector;

// internal struct to specialize for vector parameters
// TYPE = InnerValue<TVector>::Type
// SIZEOF = sizeof(TVector)
template <typename TYPE, int SIZEOF>
struct SimdParams_ {};

// internal struct to specialize for matrix parameters
template <int ROWS, int COLS, int BITS_PER_VALUE>
struct SimdMatrixParams_ {};

// struct to get the type that is in the vector
template <typename T, const int I = 0>
struct InnerValue;

#define SEQAN_DEFINE_SIMD_VECTOR_(TSimdVector, TValue, SIZEOF_VECTOR)                                           \
        typedef TValue TSimdVector __attribute__ ((__vector_size__(SIZEOF_VECTOR)));                            \
        template <> struct SimdVector<TValue, SIZEOF_VECTOR / sizeof(TValue)> {  typedef TSimdVector Type; };   \
        template <> struct InnerValue<TSimdVector>           { typedef TValue Type; };                          \
        template <> struct InnerValue<TSimdVector const>:  public InnerValue<TSimdVector> {};                   \
        template <> struct Value<TSimdVector>           { typedef TSimdVector Type; };                          \
        template <> struct Value<TSimdVector const>:  public Value<TSimdVector> {};                             \
        template <> struct LENGTH<TSimdVector>          { enum { VALUE = SIZEOF_VECTOR / sizeof(TValue) }; };   \
        template <> struct LENGTH<TSimdVector const>: public LENGTH<TSimdVector> {};                            \
        SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                            \
        SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector const)                                                      \
        SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                                      \
        template <>                                                                                             \
        SEQAN_CONCEPT_IMPL((TSimdVector),       (SimdVectorConcept));                                           \
        template <>                                                                                             \
        SEQAN_CONCEPT_IMPL((TSimdVector const), (SimdVectorConcept));

#ifdef __AVX2__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Char,     char,           32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32SChar,    signed char,    32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32UChar,    unsigned char,  32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Short,    short,          32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UShort,   unsigned short, 32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Int,       int,            32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UInt,      unsigned int,   32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int64,     __int64,        32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt64,    __uint64,       32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Float,    float,          32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Double,   double,         32)
#endif

#ifdef __SSE3__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Char,     char,           16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16SChar,    signed char,    16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UChar,    unsigned char,  16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Short,     short,          16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UShort,    unsigned short, 16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int,       int,            16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt,      unsigned int,   16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int64,     __int64,        16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt64,    __uint64,       16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Float,    float,          16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Double,   double,         16)
#endif


// vector type for DP alignments
#if defined(__AVX2__)
 typedef typename SimdVector<int16_t, 16>::Type TSimdAlign;
#elif defined(__SSE3__)
 typedef typename SimdVector<int16_t, 8>::Type TSimdAlign;
#endif

// ============================================================================
// Functions
// ============================================================================

// ============================================================================
// AVX2 wrappers (256bit vectors)
// ============================================================================

#ifdef __AVX2__

// --------------------------------------------------------------------------
// _fillVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int8_t, 32>) { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi8(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int16_t, 32>) { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi16(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int32_t, 32>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi32(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int64_t, 32>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi64x(x); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, float x,  SimdParams_<float, 32>)  { reinterpret_cast<__m256&>(vector) = _mm256_set1_ps(x); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, double x, SimdParams_<double, 32>)  { reinterpret_cast<__m256d&>(vector) = _mm256_set1_pd(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4, TValue x5,
                        TValue x6, TValue x7, TValue x8, TValue x9, TValue x10, TValue x11, TValue x12,
                        TValue x13, TValue x14, TValue x15, TValue x16, SimdParams_<int16_t, 32>)
{
    reinterpret_cast<__m256i&>(vector) = _mm256_set_epi16(x16,x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1);
}
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, SimdParams_<int32_t, 32>)
{
    reinterpret_cast<__m256i&>(vector) = _mm256_set_epi32(x8,x7,x6,x5,x4,x3,x2,x1);
}
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, SimdParams_<int16_t, 32>)
{
    //no-op, but we get unused parameter compiler warnings without this line
    reinterpret_cast<__m256i&>(vector) = _mm256_set_epi32(x8,x7,x6,x5,x4,x3,x2,x1);
}

// --------------------------------------------------------------------------
// _clearVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline void _clearVector(TSimdVector &vector, SimdParams_<TValue, 32>) { reinterpret_cast<__m256i&>(vector) = _mm256_setzero_si256(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<float, 32>) { reinterpret_cast<__m256&>(vector) = _mm256_setzero_ps(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<double, 32>) { reinterpret_cast<__m256d&>(vector) = _mm256_setzero_pd(); }

// --------------------------------------------------------------------------
// _createVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int8_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_set1_epi8(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int16_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_set1_epi16(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int32_t, 32>)  { return reinterpret_cast<TSimdVector>(_mm256_set1_epi32(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int64_t, 32>)  { return reinterpret_cast<TSimdVector>(_mm256_set1_epi64x(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(float x,  SimdParams_<float, 32>)  { return reinterpret_cast<TSimdVector>(_mm256_set1_ps(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(double x, SimdParams_<double, 32>)  { return reinterpret_cast<TSimdVector>(_mm256_set1_pd(x)); }

// --------------------------------------------------------------------------
// _cmpEq (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi8(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi32(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi64(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _cmpGt (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpgt_epi8(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpgt_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpgt_epi32(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpgt_epi64(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseOr (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseOr(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_or_si256(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseAnd (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_and_si256(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<double, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_and_pd(
        reinterpret_cast<const __m256d&>(a),
        reinterpret_cast<const __m256d&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<float, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_and_ps(
        reinterpret_cast<const __m256&>(a),
        reinterpret_cast<const __m256&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseAndNot(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_andnot_si256(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseNot (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi8(
        reinterpret_cast<const __m256i&>(a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi16(
        reinterpret_cast<const __m256i&>(a), _mm256_setzero_si256()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi32(
        reinterpret_cast<const __m256i&>(a), _mm256_setzero_si256()));

}
template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_cmpeq_epi64(
        reinterpret_cast<const __m256i&>(a), _mm256_setzero_si256()));
}

// --------------------------------------------------------------------------
// _divide (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int8_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_div_epi8(a, _mm256_set1_epi8(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int16_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_div_epi16(a, _mm256_set1_epi16(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int32_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_div_epi32(a, _mm256_set1_epi32(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int64_t, 32>) { return reinterpret_cast<TSimdVector>(_mm256_div_epi64(a, _mm256_set1_epi64x(b))); }

// --------------------------------------------------------------------------
// _add (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_epi8(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_epi32(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_epi64(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<float, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_ps(
        reinterpret_cast<const __m256&>(a),
        reinterpret_cast<const __m256&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<double, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_add_pd(
        reinterpret_cast<const __m256d&>(a),
        reinterpret_cast<const __m256d&>(b)));
}

// --------------------------------------------------------------------------
// _sub (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_sub_epi8(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_sub_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_sub_epi32(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_sub_epi64(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _mult (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_mullo_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _max (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_max_epi16(
        reinterpret_cast<const __m256i&>(a),
        reinterpret_cast<const __m256i&>(b)));
}

// --------------------------------------------------------------------------
// _blend (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, typename TValue>
inline TSimdVector _blend(TSimdVector const &a, TSimdVector const &b, TSimdVectorMask const &mask, SimdParams_<TValue, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_blendv_epi8(
        reinterpret_cast<const __m256i &>(a),
        reinterpret_cast<const __m256i &>(b),
        reinterpret_cast<const __m256i &>(mask)));
}

// --------------------------------------------------------------------------
// _storeu (256bit)
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector, typename TValue>
inline void _storeu(T * memAddr, TSimdVector &vec, SimdParams_<TValue, 32>)
{
    _mm256_storeu_si256((__m256i*)memAddr, reinterpret_cast<const __m256i&>(vec));
}

// --------------------------------------------------------------------------
// _shuffleVector (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<int8_t, 32>, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector1>(_mm256_shuffle_epi8(
        reinterpret_cast<const __m256i&>(vector),
        reinterpret_cast<const __m256i&>(indices)));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<int16_t, 32>, SimdParams_<int8_t, 16>)
{
    // copy 2nd 64bit word to 3rd, compute 2*idx
    __m256i idx = _mm256_slli_epi16(_mm256_permute4x64_epi64(_mm256_castsi128_si256(reinterpret_cast<const __m128i &>(indices)), 0x50), 1);
    // interleave with 2*idx+1 and call shuffle
    return reinterpret_cast<TSimdVector1>(_mm256_shuffle_epi8(
        reinterpret_cast<const __m256i &>(vector),
        _mm256_unpacklo_epi8(idx, _mm256_add_epi8(idx, _mm256_set1_epi8(1)))));
}

// --------------------------------------------------------------------------
// _shiftRightLogical (256bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int8_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi16(reinterpret_cast<const __m256i &>(vector), imm) & _mm256_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int16_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi16(reinterpret_cast<const __m256i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int32_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi32(reinterpret_cast<const __m256i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int64_t, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi64(reinterpret_cast<const __m256i &>(vector), imm));
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
        tmp1[i]    = _mm256_unpacklo_epi8(reinterpret_cast<const __m256i &>(matrix[2*i]), reinterpret_cast<const __m256i &>(matrix[2*i+1]));
        tmp1[i+16] = _mm256_unpackhi_epi8(reinterpret_cast<const __m256i &>(matrix[2*i]), reinterpret_cast<const __m256i &>(matrix[2*i+1]));
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
        matrix[bitRev[i]]    = reinterpret_cast<TSimdVector>(_mm256_unpacklo_epi128(tmp2[2*i],tmp2[2*i+1]));
        matrix[bitRev[i+16]] = reinterpret_cast<TSimdVector>(_mm256_unpackhi_epi128(tmp2[2*i],tmp2[2*i+1]));
    }
}

#endif  // #ifdef __AVX2__

// ============================================================================
// SSE3 wrappers (128bit vectors)
// ============================================================================

#ifdef __SSE3__

// --------------------------------------------------------------------------
// _fillVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int8_t, 16>) { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi8(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int16_t, 16>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi16(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int32_t, 16>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi32(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<int64_t, 16>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi64x(x); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, float x,  SimdParams_<float, 16>)  { reinterpret_cast<__m128&>(vector) = _mm_set1_ps(x); }
template <typename TSimdVector>
inline void _fillVector(TSimdVector &vector, double x, SimdParams_<double, 16>)  { reinterpret_cast<__m128d&>(vector) = _mm_set1_pd(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, TValue x9, TValue x10,
                        TValue x11, TValue x12, TValue x13, TValue x14, TValue x15, TValue x16, SimdParams_<int8_t, 16>)
{
    reinterpret_cast<__m128i&>(vector) = _mm_set_epi8(x16,x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1);
}
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, SimdParams_<int16_t, 16>)
{
    reinterpret_cast<__m128i&>(vector) = _mm_set_epi16(x8,x7,x6,x5,x4,x3,x2,x1);
}
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, SimdParams_<int8_t, 16>)
{
    //no-op, but we get unused parameter compiler warnings without this line
    reinterpret_cast<__m128i&>(vector) = _mm_set_epi16(x8,x7,x6,x5,x4,x3,x2,x1);
}
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                        TValue x5, TValue x6, TValue x7, TValue x8, TValue x9, TValue x10,
                        TValue x11, TValue x12, TValue x13, TValue x14, TValue x15, TValue x16, SimdParams_<int16_t, 16>)
{
    //no-op, but we get unused parameter compiler warnings without this line
    reinterpret_cast<__m128i&>(vector) = _mm_set_epi8(x16,x15,x14,x13,x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1);
}

// --------------------------------------------------------------------------
// _clearVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline void _clearVector(TSimdVector &vector, SimdParams_<TValue, 16>) { reinterpret_cast<__m128i&>(vector) = _mm_setzero_si128(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<float, 16>) { reinterpret_cast<__m128&>(vector) = _mm_setzero_ps(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<double, 16>) { reinterpret_cast<__m128d&>(vector) = _mm_setzero_pd(); }

// --------------------------------------------------------------------------
// _createVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int8_t, 16>) { return reinterpret_cast<TSimdVector>(_mm_set1_epi8(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int16_t, 16>)  { return reinterpret_cast<TSimdVector>(_mm_set1_epi16(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int32_t, 16>)  { return reinterpret_cast<TSimdVector>(_mm_set1_epi32(x)); }
template <typename TSimdVector, typename TValue>
inline TSimdVector _createVector(TValue x, SimdParams_<int64_t, 16>)  { return reinterpret_cast<TSimdVector>(_mm_set1_epi64x(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(float x,  SimdParams_<float, 16>)  { return reinterpret_cast<TSimdVector>(_mm_set1_ps(x)); }
template <typename TSimdVector>
inline TSimdVector _createVector(double x, SimdParams_<double, 16>)  { return reinterpret_cast<TSimdVector>(_mm_set1_pd(x)); }

// --------------------------------------------------------------------------
// cmpEq (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi8(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpEq(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi32(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _cmpGt (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpgt_epi8(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpgt_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _cmpGt(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpgt_epi32(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseOr (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseOr(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_or_si128(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseAnd (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_and_si128(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<double, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_and_pd(
        reinterpret_cast<const __m128d&>(a),
        reinterpret_cast<const __m128d&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseAnd(TSimdVector &a, TSimdVector &b, SimdParams_<float, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_and_ps(
        reinterpret_cast<const __m128&>(a),
        reinterpret_cast<const __m128&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseAndNot (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline TSimdVector _bitwiseAndNot(TSimdVector &a, TSimdVector &b, SimdParams_<TValue, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_andnot_si128(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _bitwiseNot (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi8(
        reinterpret_cast<const __m128i&>(a), _mm_setzero_si128()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi16(
    reinterpret_cast<const __m128i&>(a), _mm_setzero_si128()));
}

template <typename TSimdVector>
inline TSimdVector _bitwiseNot(TSimdVector &a, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_cmpeq_epi32(
        reinterpret_cast<const __m128i&>(a), _mm_setzero_si128()));
}

// --------------------------------------------------------------------------
// _divide (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int8_t, 16>) { return reinterpret_cast<TSimdVector>(_mm_div_epi8(a, _mm_set1_epi8(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int16_t, 16>){ return reinterpret_cast<TSimdVector>(_mm_div_epi16(a, _mm_set1_epi16(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int32_t, 16>) { return reinterpret_cast<TSimdVector>(_mm_div_epi32(a, _mm_set1_epi32(b))); }

template <typename TSimdVector>
inline TSimdVector _divide(TSimdVector &a, int b, SimdParams_<int64_t, 16>) { return reinterpret_cast<TSimdVector>(_mm_div_epi64(a, _mm_set1_epi64x(b))); }

// --------------------------------------------------------------------------
// _add (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_epi8(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_epi32(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_epi64(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<float, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_ps(
        reinterpret_cast<const __m128&>(a),
        reinterpret_cast<const __m128&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _add(TSimdVector &a, TSimdVector &b, SimdParams_<double, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_add_pd(
        reinterpret_cast<const __m128d&>(a),
        reinterpret_cast<const __m128d&>(b)));
}

// --------------------------------------------------------------------------
// _sub (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_sub_epi8(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_sub_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_sub_epi32(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

template <typename TSimdVector>
inline TSimdVector _sub(TSimdVector &a, TSimdVector &b, SimdParams_<int64_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_sub_epi64(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _mult (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _mult(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_mullo_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _max (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _max(TSimdVector &a, TSimdVector &b, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_max_epi16(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b)));
}

// --------------------------------------------------------------------------
// _blend (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask, typename TValue>
inline TSimdVector _blend(TSimdVector const &a, TSimdVector const &b, TSimdVectorMask const &mask, SimdParams_<TValue, 16>)
{
#ifdef __SSE4_1__
    return reinterpret_cast<TSimdVector>(_mm_blendv_epi8(
        reinterpret_cast<const __m128i&>(a),
        reinterpret_cast<const __m128i&>(b),
        reinterpret_cast<const __m128i&>(mask)));
#else
    return reinterpret_cast<TSimdVector>(_mm_or_si128(_mm_and_si128(
                                    reinterpret_cast<const __m128i&>(mask),
                                    reinterpret_cast<const __m128i&>(b)),
                                                      _mm_andnot_si128(
                                    reinterpret_cast<const __m128i&>(mask),
                                    reinterpret_cast<const __m128i&>(a))));
#endif
}

// --------------------------------------------------------------------------
// _storeu (128bit)
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector, typename TValue>
inline void _storeu(T * memAddr, TSimdVector &vec, SimdParams_<TValue, 16>)
{
    _mm_storeu_si128((__m128i*)memAddr, reinterpret_cast<const __m128i &>(vec));
}

// --------------------------------------------------------------------------
// _shuffleVector (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<int8_t, 16>, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector1>(_mm_shuffle_epi8(vector,indices));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<int16_t, 16>, SimdParams_<int8_t, 8>)
{
    __m128i idx = _mm_slli_epi16(_mm_cvtsi64_si128(reinterpret_cast<const __uint64 &>(indices)), 1);
    return reinterpret_cast<TSimdVector1>(_mm_shuffle_epi8(vector, _mm_unpacklo_epi8(idx, _mm_add_epi8(idx, _mm_set1_epi8(1)))));
}

// --------------------------------------------------------------------------
// _shiftRightLogical (128bit)
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int8_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi16(reinterpret_cast<const __m128i &>(vector), imm) & _mm_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int16_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi16(reinterpret_cast<const __m128i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int32_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi32(reinterpret_cast<const __m128i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<int64_t, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi64(reinterpret_cast<const __m128i &>(vector), imm));
}

// --------------------------------------------------------------------------
// _transposeMatrix (128bit)
// --------------------------------------------------------------------------

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
        tmp1[i]   = _mm_unpacklo_epi8(reinterpret_cast<const __m128i &>(matrix[2*i]), reinterpret_cast<const __m128i &>(matrix[2*i+1]));
        tmp1[i+8] = _mm_unpackhi_epi8(reinterpret_cast<const __m128i &>(matrix[2*i]), reinterpret_cast<const __m128i &>(matrix[2*i+1]));
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
        matrix[bitRev[i]]   = reinterpret_cast<TSimdVector>(_mm_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]));
        matrix[bitRev[i+8]] = reinterpret_cast<TSimdVector>(_mm_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]));
    }
}

#endif  // #ifdef __SSE3__

// ============================================================================
//
// INTERFACE FUNCTIONS
// - these should be used in the actual code, they will call one of the wrapper
//   functions defined above based on the vector type
//
// ============================================================================

template <int ROWS, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
transpose(TSimdVector matrix[ROWS])
{
    typedef typename InnerValue<TSimdVector>::Type TValue;
    _transposeMatrix(matrix, SimdMatrixParams_<ROWS, LENGTH<TSimdVector>::VALUE, BitsPerValue<TValue>::VALUE>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
clearVector(TSimdVector &vector)
{
    _clearVector(vector, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
createVector(TValue x)
{
    return _createVector<TSimdVector>(x, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector &vector, TValue x)
{
    _fillVector(vector, x, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                                TValue x5, TValue x6, TValue x7, TValue x8)
{
    _fillVector(vector, x1, x2, x3, x4, x5, x6, x7, x8, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector &vector, TValue x1, TValue x2, TValue x3, TValue x4,
                                       TValue x5, TValue x6, TValue x7, TValue x8,
                                       TValue x9, TValue x10, TValue x11, TValue x12,
                                       TValue x13, TValue x14, TValue x15, TValue x16)
{
    _fillVector(vector, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16,
                SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}


template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
cmpEq (TSimdVector const &a, TSimdVector const &b)
{
    return _cmpEq(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator == (TSimdVector const &a, TSimdVector const &b)
{
    return _cmpEq(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
cmpGt (TSimdVector const &a, TSimdVector const &b)
{
    return _cmpGt(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator > (TSimdVector const &a, TSimdVector const &b)
{
    return _cmpGt(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
max(TSimdVector const &a, TSimdVector const &b)
{
    return _max(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator | (TSimdVector const &a, TSimdVector const &b)
{
    return _bitwiseOr(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator |= (TSimdVector &a, TSimdVector const &b)
{
    a = a | b;
    return a;
}


template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator & (TSimdVector const &a, TSimdVector const &b)
{
    return _bitwiseAnd(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type,
                                         sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator &= (TSimdVector &a, TSimdVector const &b)
{
    a = a & b;
    return a;
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator ~ (TSimdVector const &a)
{
    return _bitwiseNot(a, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator + (TSimdVector const &a, TSimdVector const &b)
{
    return _add(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}


template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator - (TSimdVector const &a, TSimdVector const &b)
{
    return _sub(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator * (TSimdVector const &a, TSimdVector const &b)
{
    return _mult(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
andNot(TSimdVector const &a, TSimdVector const &b)
{
    return _bitwiseAndNot(a, b, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector1, typename TSimdVector2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector1> >, TSimdVector1)
shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices)
{
    return _shuffleVector(
                vector,
                indices,
                SimdParams_<typename InnerValue<TSimdVector1>::Type, sizeof(TSimdVector1)>(),
                SimdParams_<typename InnerValue<TSimdVector2>::Type, sizeof(TSimdVector2)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
shiftRightLogical(TSimdVector const &vector, const int imm)
{
    return _shiftRightLogical(vector, imm, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector, typename TSimdVectorMask>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
blend(TSimdVector const &a, TSimdVector const &b, TSimdVectorMask const & mask)
{
    return _blend(a, b, mask, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename T, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
storeu(T * memAddr, TSimdVector const &vec)
{
    _storeu(memAddr, vec, SimdParams_<typename InnerValue<TSimdVector>::Type, sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, std::ostream &)
print(std::ostream &stream, TSimdVector const &vector)
{
    stream << '<';
    for (int i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
        stream << '\t' << vector[i];
    stream << "\t>\n";
    return stream;
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
