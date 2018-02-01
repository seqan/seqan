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

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_H_

#include <utility>
#include <tuple>

#if defined(PLATFORM_WINDOWS_VS)
  /* Microsoft C/C++-compatible compiler */
  #include <intrin.h>
#elif defined(PLATFORM_GCC) && (defined(__x86_64__) || defined(__i386__))
  /* GCC-compatible compiler, targeting x86/x86-64 */
  #include <x86intrin.h>
#elif defined(SEQAN_SIMD_ENABLED)
  #pragma message "You are trying to build with -DSEQAN_SIMD_ENABLED, which might be " \
  "auto-defined if AVX or SSE was enabled (e.g. -march=native, -msse4, ...), " \
  "but we only support x86/x86-64 architectures for SIMD vectorization! " \
  "You might want to use UME::SIMD (https://github.com/edanor/umesimd) combined " \
  "with -DSEQAN_UMESIMD_ENABLED for a different SIMD backend."
#endif

namespace seqan {

#ifdef COMPILER_LINTEL
#include <type_traits>
#define SEQAN_VECTOR_CAST_(T, v) static_cast<typename std::decay<T>::type>(v)
#define SEQAN_VECTOR_CAST_LVALUE_(T, v) static_cast<T>(v)
#else
#define SEQAN_VECTOR_CAST_(T, v) reinterpret_cast<T>(v)
#define SEQAN_VECTOR_CAST_LVALUE_(T, v) reinterpret_cast<T>(v)
#endif

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Useful Macros
// ============================================================================

#define SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector)                                                 \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
getValue(TSimdVector & vector, TPosition const pos)                                                     \
{                                                                                                       \
    return vector[pos];                                                                                 \
}

#define SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                    \
template <typename TPosition>                                                                           \
inline typename Value<TSimdVector>::Type                                                                \
value(TSimdVector & vector, TPosition const pos)                                                        \
{                                                                                                       \
    return getValue(vector, pos);                                                                       \
}

#define SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                              \
template <typename TPosition, typename TValue2>                                                         \
inline void                                                                                             \
assignValue(TSimdVector & vector, TPosition const pos, TValue2 const value)                             \
{                                                                                                       \
    vector[pos] = value;                                                                                \
}

// Only include following code if simd instructions are enabled.
#ifdef SEQAN_SIMD_ENABLED

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// a metafunction returning the biggest supported SIMD vector
template <typename TValue, int LENGTH = SEQAN_SIZEOF_MAX_VECTOR / sizeof(TValue)>
struct SimdVector;

// internal struct to specialize for vector parameters
// VEC_SIZE    = Vector size := sizeof(vec)
// LENGTH      = number of elements := VEC_SIZE / sizeof(InnerValue<TVec>::Type)
// SCALAR_TYPE = the scalar type of the vector (maybe optional, if the type
//               doesn't matter for the operation)
template <int VEC_SIZE, int LENGTH = 0, typename SCALAR_TYPE = void>
struct SimdParams_
{};

// internal traits meta-function to capture correct the mask type.
template <typename TSimdVector, typename TSimdParams>
struct SimdVectorTraits
{
    using MaskType = TSimdVector;
};

// internal struct to specialize for matrix parameters
template <int ROWS, int COLS, int BITS_PER_VALUE>
struct SimdMatrixParams_
{};

#define SEQAN_DEFINE_SIMD_VECTOR_(TSimdVector, TValue, SIZEOF_VECTOR)                                           \
        typedef TValue TSimdVector __attribute__ ((__vector_size__(SIZEOF_VECTOR)));                            \
        template <> struct SimdVector<TValue, SIZEOF_VECTOR / sizeof(TValue)> {  typedef TSimdVector Type; };   \
        template <> struct Value<TSimdVector>           { typedef TValue Type; };                               \
        template <> struct Value<TSimdVector const>:  public Value<TSimdVector> {};                             \
        template <> struct LENGTH<TSimdVector>          { enum { VALUE = SIZEOF_VECTOR / sizeof(TValue) }; };   \
        template <> struct LENGTH<TSimdVector const>: public LENGTH<TSimdVector> {};                            \
        SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector)                                                         \
        SEQAN_DEFINE_SIMD_VECTOR_GETVALUE_(TSimdVector const)                                                   \
        SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector)                                                            \
        SEQAN_DEFINE_SIMD_VECTOR_VALUE_(TSimdVector const)                                                      \
        SEQAN_DEFINE_SIMD_VECTOR_ASSIGNVALUE_(TSimdVector)                                                      \
        template <>                                                                                             \
        SEQAN_CONCEPT_IMPL((TSimdVector),       (SimdVectorConcept));                                           \
        template <>                                                                                             \
        SEQAN_CONCEPT_IMPL((TSimdVector const), (SimdVectorConcept));
#endif  // SEQAN_SIMD_ENABLED

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_IMPL_H_
