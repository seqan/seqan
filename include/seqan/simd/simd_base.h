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
SEQAN_CONCEPT(SimdVectorConcept, (TSimdVector)) {
    typedef typename Reference<TSimdVector>::Type TReference;

    TSimdVector a;

    SEQAN_CONCEPT_USAGE(SimdVectorConcept)
    {
        static_assert(IsSameType<decltype(a[0]), TReference>::VALUE, "Type of a[] should be the same as the reference type of a.");
    }
};

template <typename TSimdVector, typename TIsSimdVec>
struct SimdSwizzleVectorImpl;

template <typename TSimdVector>
struct SimdSwizzleVector : SimdSwizzleVectorImpl<TSimdVector, typename Is<SimdVectorConcept<TSimdVector> >::Type >
{
};

/**
 * ```
 * getValue(a, pos);
 *
 * // same as
 *
 * a[pos];
 * ```
 */
template <typename TSimdVector, typename TPosition>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename Value<TSimdVector>::Type)
getValue(TSimdVector const &vector, TPosition const pos);

/**
 * ```
 * value(a, pos);
 *
 * // same as
 *
 * a[pos];
 * ```
 */
template <typename TSimdVector, typename TPosition>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename Value<TSimdVector>::Type)
value(TSimdVector const &vector, TPosition const pos);

/**
 * ```
 * assignValue(a, pos, value);
 *
 * // same as
 *
 * a[pos] = value;
 * ```
 */
template <typename TSimdVector, typename TPosition, typename TValue2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
assignValue(TSimdVector &vector, TPosition const pos, TValue2 const value);

template <int ROWS, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
transpose(TSimdVector matrix[ROWS]);

/**
 * ```
 * clearVector(a);
 *
 * // same as
 *
 * for(auto i = 0u; i < LENGTH; ++i)
 *     c[i] = 0;
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
clearVector(TSimdVector &vector);

/**
 * ```
 * auto c = createVector<SimdVector4Int>(a);
 *
 * // same as
 *
 * for(auto i = 0u; i < LENGTH; ++i)
 *     c[i] = a;
 * ```
 */
template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
createVector(TValue x);

/**
 * ```
 * fillVector(a, 1, 3, 23, 1337);
 *
 * // same as
 *
 * a[0] = 1;
 * a[1] = 3;
 * a[2] = 13;
 * a[3] = 1337;
 * ```
 */
template <typename TSimdVector, typename ...TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector &vector, TValue const... args);

/**
 * ```
 * c = cmpEq(a, b);
 *
 * // same as
 *
 * c = a == b;
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
cmpEq (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator == (TSimdVector const &a, TSimdVector const &b);

/**
 * ```
 * c = cmpGt(a, b);
 *
 * // same as
 *
 * c = a > b;
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
cmpGt (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator > (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
max(TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator | (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator |= (TSimdVector &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator & (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator &= (TSimdVector &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator ~ (TSimdVector const &a);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator + (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator - (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator * (TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator/ (TSimdVector const &a, TSimdVector const &b);

/**
 * ```
 * c = andNot(a, b);
 *
 * // same as
 *
 * for(auto i = 0u; i < LENGTH; ++i)
 *     c[i] = (~a[i]) & b[i];
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
andNot(TSimdVector const &a, TSimdVector const &b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
shiftRightLogical(TSimdVector const &vector, const int imm);

template <typename TSimdVector, typename TSimdVectorMask>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
blend(TSimdVector const &a, TSimdVector const &b, TSimdVectorMask const & mask);

/**
 * Unaligned store, i.e. memAddr does not need to be aligned (e.g. SEE4.2 16byte
 * aligned, AVX2 32byte aligned).
 */
template <typename T, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
storeu(T * memAddr, TSimdVector const &vec);

/**
 * Aligned load, i.e. memAddr MUST be aligned (e.g. SEE4.2 16byte
 * aligned, AVX2 32byte aligned).
 */
template <typename TSimdVector, typename T>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
load(T const * memAddr);

template <typename TValue, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
gather(TValue const * memAddr, TSimdVector const & idx);

template <typename TSimdVector1, typename TSimdVector2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector1> >, TSimdVector1)
shuffleVector(TSimdVector1 const &vector, TSimdVector2 const &indices);

// NOTE(rmaerker): Make this function available, also if SIMD is not enabled.
template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_DISABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
createVector(TValue x)
{
    return x;
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_H_
