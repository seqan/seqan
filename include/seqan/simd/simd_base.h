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
struct Value<SimdVector<TValue, LENGTH> const> :
    public Value<SimdVector<TValue, LENGTH> >
{};

template <typename TValue, int LENGTH_>
struct LENGTH<SimdVector<TValue, LENGTH_> >
{
    enum { VALUE = LENGTH_ };
};

template <typename TValue, int LENGTH_>
struct LENGTH<SimdVector<TValue, LENGTH_> const> :
    public LENGTH<SimdVector<TValue, LENGTH_> >
{};

// define a concept and its models
// they allow us to define generic vector functions
SEQAN_CONCEPT(SimdMaskVectorConcept, (TSimdMaskVector))
{
    typedef typename Reference<TSimdMaskVector>::Type TReference;

    TSimdMaskVector a;

    SEQAN_CONCEPT_USAGE(SimdMaskVectorConcept)
    {
        static_assert(IsSameType<decltype(a[0]), TReference>::VALUE, "Type of a[] should be the same as the reference type of a.");
    }
};

SEQAN_CONCEPT_REFINE(SimdVectorConcept, (TSimdVector), (SimdMaskVectorConcept))
{
    SEQAN_CONCEPT_USAGE(SimdVectorConcept)
    {}
};

template <typename TSimdVector, typename TIsSimdVec>
struct SimdMaskVectorImpl {
    using Type = Nothing;
};

/**
 * SimdMaskVector is the return type of all logical operations of simd vectors
 * like comparisons.
 *
 * ```
 * using TSimdVector = SimdVector<uint32_t, 4>::Type;
 * using TSimdMaskVector = SimdMaskVector<TSimdVector>::Type;
 *
 * TSimdVector vec1 {2, 4, 8, 16}, vec2 {16, 8, 4, 2};
 * TSimdMaskVector cmp = vec1 > vec2; // cmp = {false, false, true, true}
 * ```
 */
template <typename TSimdVector>
struct SimdMaskVector : SimdMaskVectorImpl<TSimdVector, typename Is<SimdVectorConcept<TSimdVector> >::Type >
{
};

template <typename TSimdVector, typename TIsSimdVec>
struct SimdSwizzleVectorImpl;

/**
 * SimdSwizzleVector is needed for shuffleVector() as index type.
 *
 * ```
 * using TSimdVector = SimdVector<uint32_t, 4>::Type;
 * using TSimdSwizzleVector = SimdSwizzleVector<TSimdVector>::Type;
 *
 * TSimdVector vec {2, 4, 8, 16}, res;
 * TSimdSwizzleVector swizzle {3, 2, 0, 2};
 *
 * res = shuffleVector(vec, swizzle); // res = {16, 8, 2, 8}
 * ```
 */
template <typename TSimdVector>
struct SimdSwizzleVector : SimdSwizzleVectorImpl<TSimdVector, typename Is<SimdVectorConcept<TSimdVector> >::Type >
{};

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
getValue(TSimdVector const & vector, TPosition const pos);

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
value(TSimdVector const & vector, TPosition const pos);

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
assignValue(TSimdVector & vector, TPosition const pos, TValue2 const value);

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
clearVector(TSimdVector & vector);

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
createVector(TValue const x);

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
fillVector(TSimdVector & vector, TValue const... args);

/**
 * ```
 * auto c = cmpEq(a, b);
 *
 * // same as
 *
 * auto c = a == b;
 * ```
 *
 * NOTE:
 * The type of c might change from unsigned to signed if auto is used
 *
 * ```
 * using TSimdVector = SimdVector<uint32_t, 4>::Type;
 * TSimdVector a, b;
 *
 * auto c = a == b; // type of c might change to SimdVector<int32_t, 4>::Type
 * TSimdVector d = a == b; // has the same type
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpEq (TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator==(TSimdVector const & a, TSimdVector const & b);

/**
 * ```
 * auto c = cmpGt(a, b);
 *
 * // same as
 *
 * auto c = a > b;
 * ```
 *
 * NOTE:
 * The type of c might change from unsigned to signed if auto is used
 *
 * ```
 * using TSimdVector = SimdVector<uint32_t, 4>::Type;
 * using TSimdMaskVector = SimdMaskVector<TSimdVector>::Type;
 * TSimdVector a, b;
 *
 * auto c = a > b; // type of c might change to SimdVector<int32_t, 4>::Type
 * TSimdMaskVector d = a > b; // has the same type
 * ```
 */
template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpGt (TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator>(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
max(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
min(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator|(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator|=(TSimdVector & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator&(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator&=(TSimdVector & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator~(TSimdVector const & a);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator+(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator-(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator*(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator/(TSimdVector const & a, TSimdVector const & b);

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
andNot(TSimdVector const & a, TSimdVector const & b);

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
shiftRightLogical(TSimdVector const & vector, const int imm);

template <typename TSimdVector, typename TSimdVectorMask>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask);

/**
 * Unaligned store, i.e. memAddr does not need to be aligned (e.g. SEE4.2 16byte
 * aligned, AVX2 32byte aligned).
 */
template <typename T, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
storeu(T * memAddr, TSimdVector const & vec);

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
shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices);

// NOTE(rmaerker): Make this function available, also if SIMD is not enabled.
template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<NumberConcept<TSimdVector>>, TSimdVector)
createVector(TValue const x)
{
    return x;
}

// --------------------------------------------------------------------------
// Function print()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdMaskVectorConcept<TSimdVector> >, std::ostream &)
print(std::ostream & stream, TSimdVector const & vector)
{
    stream << '<';
    for (int i = 0; i < LENGTH<TSimdVector>::VALUE; ++i)
        stream << '\t' << vector[i];
    stream << "\t>\n";
    return stream;
}

}  // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_H_
