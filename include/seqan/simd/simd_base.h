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
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================
// generic SIMD interface for SSE3 / AVX2
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_H_

// Define global macro to check if simd instructions are enabled.
#define SEQAN_SIMD_ENABLED 1

namespace seqan
{

// a metafunction returning the biggest supported SIMD vector
template <typename TValue, int LENGTH>
struct SimdVector;

template <typename TValue, int LENGTH>
struct Value<SimdVector<TValue, LENGTH> >
{
    typedef TValue Type;
};

template <typename TValue, int LENGTH>
struct Value<SimdVector<TValue, LENGTH> const>:
    public Value<SimdVector<TValue, LENGTH> >
{
};

template <typename TValue, int LENGTH_>
struct LENGTH<SimdVector<TValue, LENGTH_> >
{
    enum { VALUE = LENGTH_ };
};

template <typename TValue, int LENGTH_>
struct LENGTH<SimdVector<TValue, LENGTH_> const>:
    public LENGTH<SimdVector<TValue, LENGTH_> >
{
};

// define a concept and its models
// they allow us to define generic vector functions
SEQAN_CONCEPT(SimdVectorConcept, (T)) {};

template <typename TSimdVector, typename TPosition>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename Value<TSimdVector>::Type)
getValue(TSimdVector &vector, TPosition pos);

template <typename TSimdVector, typename TPosition>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename Value<TSimdVector>::Type)
value(TSimdVector &vector, TPosition pos);

template <typename TSimdVector, typename TPosition, typename TValue2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
assignValue(TSimdVector &vector, TPosition pos, TValue2 value);

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_H_
