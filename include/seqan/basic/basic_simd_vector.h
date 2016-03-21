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
// ==========================================================================
// generic SIMD interface for SSE4/AVX
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_

#if defined(__SSE4_1__) || defined(__AVX__)
 #include <immintrin.h>
#else
// SSE4.1 or greater required
#ifdef _MSC_VER
 #pragma message("SSE4.1 instruction set not enabled")
#else
 #warning "SSE4.1 instruction set not enabled"
#endif  // _MSC_VER
#endif


namespace seqan {

// ============================================================================
// Useful Macros
// ============================================================================

#define SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector)                                                 \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
getValue(TSimdVector &vector, TPosition pos)                                                            \
{                                                                                                       \
/*                                                                                                      \
    typedef typename Value<TSimdVector>::Type TValue;                                                   \
    TValue val = (reinterpret_cast<TValue*>(&vector))[pos];                                    \
    return val;                                                                                         \
*/                                                                                                      \
    return vector[pos];                                                                                 \
}


#define SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                    \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
value(TSimdVector &vector, TPosition pos)                                                               \
{                                                                                                       \
    return getValue(vector, pos);                                                                       \
}

#define SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                              \
template <typename TPosition, typename TValue2>                                                         \
inline void                                                                                             \
assignValue(TSimdVector &vector, TPosition pos, TValue2 value)                                          \
{                                                                                                       \
/*                                                                                                      \
    typedef typename Value<TSimdVector>::Type TValue;                                                   \
    (reinterpret_cast<TValue*>(&vector))[pos] = value;                                                  \
*/                                                                                                      \
    vector[pos] = value;                                                                                \
}


// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// define a concept and its models
// they allow us to define generic vector functions
SEQAN_CONCEPT(SimdVectorConcept, (T)) {};

#if defined(__AVX2__)
#define SEQAN_SIZEOF_MAX_VECTOR 32
#elif defined(__SSE3__)
#define SEQAN_SIZEOF_MAX_VECTOR 16
#else
#define SEQAN_SIZEOF_MAX_VECTOR 8
#endif

// a metafunction returning the biggest supported SIMD vector
template <typename TValue, int LENGTH = SEQAN_SIZEOF_MAX_VECTOR / sizeof(TValue)>
struct SimdVector;

// internal struct to specialize for vector parameters (SIZEOF=sizeof(TVector), LENGTH=LENGTH<TVector>::VALUE)
template <int SIZEOF, int LENGTH = 0>
struct SimdParams_ {};

// internal struct to specialize for matrix parameters
template <int ROWS, int COLS, int BITS_PER_VALUE>
struct SimdMatrixParams_
{
};


#define SEQAN_DEFINE_SIMD_VECTOR_(TSimdVector, TValue, SIZEOF_VECTOR)                                       \
    typedef TValue TSimdVector __attribute__ ((__vector_size__ (SIZEOF_VECTOR)));                           \
    template <> struct SimdVector<TValue, SIZEOF_VECTOR / sizeof(TValue)> { typedef TSimdVector Type; };    \
    template <> struct Value<TSimdVector>           { typedef TValue Type; };                               \
    template <> struct LENGTH<TSimdVector>          { enum { VALUE = SIZEOF_VECTOR / sizeof(TValue) }; };   \
    template <> struct Value<TSimdVector const>:  public Value<TSimdVector> {};                             \
    template <> struct LENGTH<TSimdVector const>: public LENGTH<TSimdVector> {};                            \
    SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector const)                                                   \
    SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                            \
    SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector const)                                                      \
    SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                                      \
    template <>                                                                                             \
    SEQAN_CONCEPT_IMPL((TSimdVector),       (SimdVectorConcept));                                           \
    template <>                                                                                             \
    SEQAN_CONCEPT_IMPL((TSimdVector const), (SimdVectorConcept));

#ifdef __AVX__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32Char,     char,           32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32SChar,    signed char,    32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector32UChar,    unsigned char,  32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Short,    short,          32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UShort,   unsigned short, 32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Int,       int,            32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UInt,      unsigned int,   32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int64,     int64_t,        32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt64,    uint64_t,       32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Float,     float,          32)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Double,    double,         32)
#endif

#ifdef __SSE3__
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Char,      char,           8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8SChar,     signed char,    8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UChar,     unsigned char,  8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Short,     short,          8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UShort,    unsigned short, 8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int,       int,            8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt,      unsigned int,   8)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Float,     float,          8)

SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16Char,     char,           16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16SChar,    signed char,    16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector16UChar,    unsigned char,  16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8Short,     short,          16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector8UShort,    unsigned short, 16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Int,       int,            16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4UInt,      unsigned int,   16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Int64,     int64_t,        16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2UInt64,    uint64_t,       16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector4Float,     float,          16)
SEQAN_DEFINE_SIMD_VECTOR_(SimdVector2Double,    double,         16)
#endif

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// AVX/AVX2 wrappers
// --------------------------------------------------------------------------

#ifdef __AVX__

template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<32, 32>) { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi8(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<32, 16>) { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi16(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<32, 8>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi32(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<32, 4>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_epi64x(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, float x,  SimdParams_<32, 8>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_ps(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, double x, SimdParams_<32, 4>)  { reinterpret_cast<__m256i&>(vector) = _mm256_set1_pd(x); }

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, L>) { reinterpret_cast<__m256i&>(vector) = _mm256_setzero_si256(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, 8>) { reinterpret_cast<__m256&>(vector) = _mm256_setzero_ps(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<32, 4>) { reinterpret_cast<__m256d&>(vector) = _mm256_setzero_pd(); }

#ifdef __AVX2__

template <typename TSimdVector, int L>
inline TSimdVector _blend(TSimdVector const &a, TSimdVector const &b, TSimdVector const &mask, SimdParams_<32, L>)
{
    return reinterpret_cast<TSimdVector>(_mm256_blendv_epi8(
        reinterpret_cast<const __m256i &>(a),
        reinterpret_cast<const __m256i &>(b),
        reinterpret_cast<const __m256i &>(mask)));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 32>, SimdParams_<32, 32>)
{
    return reinterpret_cast<TSimdVector1>(_mm256_shuffle_epi8(
        reinterpret_cast<const __m256i &>(vector),
        reinterpret_cast<const __m256i &>(indices)));
}
template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 16>, SimdParams_<16, 16>)
{
    // copy 2nd 64bit word to 3rd, compute 2*idx
    __m256i idx = _mm256_slli_epi16(_mm256_permute4x64_epi64(_mm256_castsi128_si256(reinterpret_cast<const __m128i &>(indices)), 0x50), 1);
    // interleave with 2*idx+1 and call shuffle
    return reinterpret_cast<TSimdVector1>(_mm256_shuffle_epi8(
        reinterpret_cast<const __m256i &>(vector),
         _mm256_unpacklo_epi8(idx, _mm256_add_epi8(idx, _mm256_set1_epi8(1)))));
}

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 32>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi16(reinterpret_cast<const __m256i &>(vector), imm) & _mm256_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi16(reinterpret_cast<const __m256i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 8>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi32(reinterpret_cast<const __m256i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<32, 4>)
{
    return reinterpret_cast<TSimdVector>(_mm256_srli_epi64(reinterpret_cast<const __m256i &>(vector), imm));
}

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

#else   // #ifdef __AVX2__
template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<32, 32>, SimdParams_<32, 32>)
{
    return reinterpret_cast<TSimdVector1>(_mm256_permute2f128_si256(
        _mm256_castsi128_si256 (_mm_shuffle_epi8(
                                    _mm256_castsi256_si128(reinterpret_cast<const __m256i &>(vector)),
                                    _mm256_castsi256_si128(reinterpret_cast<const __m256i &>(indices)))),
        _mm256_castsi128_si256 (_mm_shuffle_epi8(
                                    _mm256_castsi256_si128(reinterpret_cast<const __m256i &>(vector)),
                                    _mm256_extractf128_si256(reinterpret_cast<const __m256i &>(indices), 1))),
        0x20));
}

inline SimdVector32Char   shiftRightLogical(SimdVector32Char   const &vector, const int imm)
{
    return reinterpret_cast<SimdVector32Char>(_mm256_permute2f128_si256(
        _mm256_castsi128_si256 (_mm_srli_epi16(
                                    _mm256_castsi256_si128(reinterpret_cast<const __m256i &>(vector)),
                                    imm)),
        _mm256_castsi128_si256 (_mm_srli_epi16(
                                    _mm256_extractf128_si256(reinterpret_cast<const __m256i &>(vector), 1),
                                    imm)),
        0x20) & _mm256_set1_epi8(0xff >> imm));
}

#endif  // #ifdef __AVX2__


template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline _testAllZeros(TSimdVector const &vector, TSimdVector const &mask, SimdParams_<32>)
{
#ifdef __AVX2__
    return _mm256_testz_si256(vector, mask);
#else   // #ifdef __AVX2__
    return
        _mm_testz_si128(_mm256_castsi256_si128(vector), _mm256_castsi256_si128(mask)) &
        _mm_testz_si128(_mm256_extractf128_si256(vector, 1), _mm256_extractf128_si256(mask, 1));
#endif  // #ifdef __AVX2__
}

template <typename TSimdVector>
inline int _testAllOnes(TSimdVector const &vector, SimdParams_<32>)
{
    __m256i vec = reinterpret_cast<const __m256i &>(vector);
#ifdef __AVX2__
    return _mm256_testc_si256(vec, _mm256_cmpeq_epi32(vec, vec));
#else   // #ifdef __AVX2__
    return
        _mm_test_all_ones(_mm256_castsi256_si128(vec)) &
        _mm_test_all_ones(_mm256_extractf128_si256(vec, 1));
#endif  // #ifdef __AVX2__
}

#endif  // #ifdef __AVX__


// --------------------------------------------------------------------------
// SSE3 wrappers
// --------------------------------------------------------------------------

#ifdef __SSE3__

template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<16, 16>) { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi8(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<16, 8>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi16(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<16, 4>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi32(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, TValue x, SimdParams_<16, 2>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_epi64x(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, float x,  SimdParams_<16, 4>)   { reinterpret_cast<__m128i&>(vector) = _mm_set1_ps(x); }
template <typename TSimdVector, typename TValue>
inline void _fillVector(TSimdVector &vector, double x, SimdParams_<16, 2>)  { reinterpret_cast<__m128i&>(vector) = _mm_set1_pd(x); }

template <typename TSimdVector, int L>
inline void _clearVector(TSimdVector &vector, SimdParams_<16, L>) { reinterpret_cast<__m128i&>(vector) = _mm_setzero_si128(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<16, 4>)  { reinterpret_cast<__m128&>(vector) = _mm_setzero_ps(); }
template <typename TSimdVector>
inline void _clearVector(TSimdVector &vector, SimdParams_<16, 2>)  { reinterpret_cast<__m128d&>(vector) = _mm_setzero_pd(); }


template <typename TSimdVector, int L>
inline TSimdVector _blend(TSimdVector const &a, TSimdVector const &b, TSimdVector const &mask, SimdParams_<16, L>)
{
    return reinterpret_cast<TSimdVector>(_mm_blendv_epi8(
        reinterpret_cast<const __m128i &>(a),
        reinterpret_cast<const __m128i &>(b),
        reinterpret_cast<const __m128i &>(mask)));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<16, 16>, SimdParams_<16, 16>)
{
    return reinterpret_cast<TSimdVector1>(_mm_shuffle_epi8(
        reinterpret_cast<const __m128i &>(vector),
        reinterpret_cast<const __m128i &>(indices)));
}

template <typename TSimdVector1, typename TSimdVector2>
inline TSimdVector1
_shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices, SimdParams_<16, 8>, SimdParams_<8, 8>)
{
#if defined(SEQAN_IS_32_BIT)
    __m128i idx = _mm_slli_epi16(_mm_unpacklo_epi32(_mm_cvtsi32_si128(reinterpret_cast<const uint32_t &>(indices)),
                                                    _mm_cvtsi32_si128(reinterpret_cast<const uint64_t &>(indices) >> 32)), 1);
#else
    __m128i idx = _mm_slli_epi16(_mm_cvtsi64_si128(reinterpret_cast<const uint64_t &>(indices)), 1);
#endif  // defined(SEQAN_IS_32_BIT)
    return reinterpret_cast<TSimdVector1>(_mm_shuffle_epi8(
        reinterpret_cast<const __m128i &>(vector),
        _mm_unpacklo_epi8(idx, _mm_add_epi8(idx, _mm_set1_epi8(1)))));
}

template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<16, 16>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi16(reinterpret_cast<const __m128i &>(vector), imm) & _mm_set1_epi8(0xff >> imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<16, 8>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi16(reinterpret_cast<const __m128i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<16, 4>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi32(reinterpret_cast<const __m128i &>(vector), imm));
}
template <typename TSimdVector>
inline TSimdVector _shiftRightLogical(TSimdVector const &vector, const int imm, SimdParams_<16, 2>)
{
    return reinterpret_cast<TSimdVector>(_mm_srli_epi64(reinterpret_cast<const __m128i &>(vector), imm));
}



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
        tmp1[i]   = _mm_unpacklo_pi8(reinterpret_cast<const __m64 &>(matrix[2*i]), reinterpret_cast<const __m64 &>(matrix[2*i+1]));
        tmp1[i+4] = _mm_unpackhi_pi8(reinterpret_cast<const __m64 &>(matrix[2*i]), reinterpret_cast<const __m64 &>(matrix[2*i+1]));
    }
    __m64 tmp2[8];
    for (int i = 0; i < 4; ++i)
    {
        tmp2[i]   = _mm_unpacklo_pi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+4] = _mm_unpackhi_pi16(tmp1[2*i], tmp1[2*i+1]);
    }
    for (int i = 0; i < 4; ++i)
    {
        matrix[bitRev[i]]   = reinterpret_cast<TSimdVector>(_mm_unpacklo_pi32(tmp2[2*i], tmp2[2*i+1]));
        matrix[bitRev[i+4]] = reinterpret_cast<TSimdVector>(_mm_unpackhi_pi32(tmp2[2*i], tmp2[2*i+1]));
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

#ifdef __SSE4_1__
template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline _testAllZeros(TSimdVector const &vector, TSimdVector const &mask, SimdParams_<16>)
{
    return _mm_testz_si128(vector, mask);
}

template <typename TSimdVector>
inline int _testAllOnes(TSimdVector const &vector, SimdParams_<16>)
{
    return _mm_test_all_ones(reinterpret_cast<const __m128i &>(vector));
}


#endif  // #ifdef __SSE3__
//#endif  // #ifdef __AVX__
#endif


template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllZeros(TSimdVector const &vector, TSimdVector const &mask)
{
    return _testAllZeros(vector, mask, SimdParams_<sizeof(TSimdVector)>());
}

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    int)
inline testAllZeros(TSimdVector const &vector)
{
    return _testAllZeros(vector, vector, SimdParams_<sizeof(TSimdVector)>());
}

template <typename TSimdVector>
inline int _testAllOnes(TSimdVector const &vector, True)
{
    return _testAllOnes(vector, SimdParams_<sizeof(vector)>());
}

template <int ROWS, typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    void)
inline transpose(TSimdVector matrix[ROWS])
{
    typedef typename Value<TSimdVector>::Type TValue;
    _transposeMatrix(matrix, SimdMatrixParams_<
                                    ROWS,
                                    LENGTH<TSimdVector>::VALUE,
                                    BitsPerValue<TValue>::VALUE>());
}


template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    void)
inline clearVector(TSimdVector &vector)
{
    _clearVector(vector, SimdParams_<sizeof(TSimdVector), LENGTH<TSimdVector>::VALUE>());
}


template <typename TSimdVector, typename TValue>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    void)
inline fillVector(TSimdVector &vector, TValue x)
{
    _fillVector(vector, x, SimdParams_<sizeof(TSimdVector), LENGTH<TSimdVector>::VALUE>());
}

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    TSimdVector)
inline blend(TSimdVector const &a, TSimdVector const &b, TSimdVector const &mask)
{
    return _blend(a, b, mask, SimdParams_<sizeof(TSimdVector), LENGTH<TSimdVector>::VALUE>());
}

template <typename TSimdVector1, typename TSimdVector2>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector1> >,
    TSimdVector1)
inline shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices)
{
    return _shuffleVector(
                vector,
                indices,
                SimdParams_<sizeof(TSimdVector1), LENGTH<TSimdVector1>::VALUE>(),
                SimdParams_<sizeof(TSimdVector2), LENGTH<TSimdVector2>::VALUE>());
}

template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    TSimdVector)
inline shiftRightLogical(TSimdVector const &vector, const int imm)
{
    return _shiftRightLogical(vector, imm, SimdParams_<sizeof(TSimdVector), LENGTH<TSimdVector>::VALUE>());
}


template <typename TSimdVector>
SEQAN_FUNC_ENABLE_IF(
    Is<SimdVectorConcept<TSimdVector> >,
    std::ostream &)
inline print(std::ostream &stream, TSimdVector const &vector)
{
    stream << '<';
    for (int i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
        stream << '\t' << (unsigned)vector[i];
    stream << "\t>";
    return stream;
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_BASIC_SIMD_VECTOR_H_
