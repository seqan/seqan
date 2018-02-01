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
// SIMD implementation of umesimd
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_UMESIMD_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_UMESIMD_IMPL_H_

#include "umesimd/UMESimd.h"

namespace seqan
{

template <typename TSimdVector>
struct SimdMaskVectorImpl<TSimdVector, True>
{
    using Type = typename UME::SIMD::SIMDTraits<TSimdVector>::MASK_T;
};

template <typename TSimdVector>
struct SimdSwizzleVectorImpl<TSimdVector, True>
{
    using Type = typename UME::SIMD::SIMDTraits<TSimdVector>::SWIZZLE_T;
};

template <typename TValue, int LENGTH>
struct SimdVector
{
    typedef UME::SIMD::SIMDVec<TValue, LENGTH> Type;
};

// // 64 bit
// using SimdVector8Char   = UME::SIMD::SIMDVec<char, 8>;
using SimdVector8SChar  = UME::SIMD::SIMDVec<signed char, 8>;
using SimdVector8UChar  = UME::SIMD::SIMDVec<unsigned char, 8>;
using SimdVector4Short  = UME::SIMD::SIMDVec<short, 4>;
using SimdVector4UShort = UME::SIMD::SIMDVec<unsigned short, 4>;
using SimdVector2Int    = UME::SIMD::SIMDVec<int, 2>;
using SimdVector2UInt   = UME::SIMD::SIMDVec<unsigned int, 2>;

// 128 bit
// using SimdVector16Char  = UME::SIMD::SIMDVec<char, 16>;
using SimdVector16SChar = UME::SIMD::SIMDVec<signed char, 16>;
using SimdVector16UChar = UME::SIMD::SIMDVec<unsigned char, 16>;
using SimdVector8Short  = UME::SIMD::SIMDVec<short, 8>;
using SimdVector8UShort = UME::SIMD::SIMDVec<unsigned short, 8>;
using SimdVector4Int    = UME::SIMD::SIMDVec<int, 4>;
using SimdVector4UInt   = UME::SIMD::SIMDVec<unsigned int, 4>;
using SimdVector2Int64  = UME::SIMD::SIMDVec<int64_t, 2>;
using SimdVector2UInt64 = UME::SIMD::SIMDVec<uint64_t, 2>;

// 256 bit
// using SimdVector32Char   = UME::SIMD::SIMDVec<char, 32>;
using SimdVector32SChar  = UME::SIMD::SIMDVec<signed char, 32>;
using SimdVector32UChar  = UME::SIMD::SIMDVec<unsigned char, 32>;
using SimdVector16Short  = UME::SIMD::SIMDVec<short, 16>;
using SimdVector16UShort = UME::SIMD::SIMDVec<unsigned short, 16>;
using SimdVector8Int     = UME::SIMD::SIMDVec<int, 8>;
using SimdVector8UInt    = UME::SIMD::SIMDVec<unsigned int, 8>;
using SimdVector4Int64   = UME::SIMD::SIMDVec<int64_t, 4>;
using SimdVector4UInt64  = UME::SIMD::SIMDVec<uint64_t, 4>;

// 512 bit
// using SimdVector64Char   = UME::SIMD::SIMDVec<char, 64>;
using SimdVector64SChar  = UME::SIMD::SIMDVec<signed char, 64>;
using SimdVector64UChar  = UME::SIMD::SIMDVec<unsigned char, 64>;
using SimdVector32Short  = UME::SIMD::SIMDVec<short, 32>;
using SimdVector32UShort = UME::SIMD::SIMDVec<unsigned short, 32>;
using SimdVector16Int    = UME::SIMD::SIMDVec<int, 16>;
using SimdVector16UInt   = UME::SIMD::SIMDVec<unsigned int, 16>;
using SimdVector8Int64   = UME::SIMD::SIMDVec<int64_t, 8>;
using SimdVector8UInt64  = UME::SIMD::SIMDVec<uint64_t, 8>;

// ============================================================================
// SIMDMaskVector
// ============================================================================

template <uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVecMask<LENGTH>),       (SimdMaskVectorConcept));

template <uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVecMask<LENGTH> const), (SimdMaskVectorConcept));

template <uint32_t LENGTH>
struct Value<UME::SIMD::SIMDVecMask<LENGTH> >
{
    typedef bool Type;
};

template <uint32_t LENGTH_>
struct LENGTH<UME::SIMD::SIMDVecMask<LENGTH_> >
{
    enum { VALUE = LENGTH_ };
};

template <uint32_t LENGTH, typename TPosition>
inline typename Value<UME::SIMD::SIMDVecMask<LENGTH> >::Type
getValue(UME::SIMD::SIMDVecMask<LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <uint32_t LENGTH, typename TPosition>
inline typename Value<UME::SIMD::SIMDVecMask<LENGTH> >::Type
value(UME::SIMD::SIMDVecMask<LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <uint32_t LENGTH, typename TPosition, typename TValue2>
inline void
assignValue(UME::SIMD::SIMDVecMask<LENGTH> &vector, TPosition const pos, TValue2 const value)
{
    vector.insert(pos, value);
}

// ============================================================================
// SIMDSwizzle
// ============================================================================

template <uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDSwizzle<LENGTH>),       (SimdVectorConcept));

template <uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDSwizzle<LENGTH> const), (SimdVectorConcept));

template <uint32_t LENGTH>
struct Value<UME::SIMD::SIMDSwizzle<LENGTH> >
{
    typedef uint32_t Type;
};

template <uint32_t LENGTH_>
struct LENGTH<UME::SIMD::SIMDSwizzle<LENGTH_> >
{
    enum { VALUE = LENGTH_ };
};

template <uint32_t LENGTH, typename TPosition>
inline typename Value<UME::SIMD::SIMDSwizzle<LENGTH> >::Type
getValue(UME::SIMD::SIMDSwizzle<LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <uint32_t LENGTH, typename TPosition>
inline typename Value<UME::SIMD::SIMDSwizzle<LENGTH> >::Type
value(UME::SIMD::SIMDSwizzle<LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <uint32_t LENGTH, typename TPosition, typename TValue2>
inline void
assignValue(UME::SIMD::SIMDSwizzle<LENGTH> &vector, TPosition const pos, TValue2 const value)
{
    vector.insert(pos, value);
}

// ============================================================================
// SIMDVec_u
// ============================================================================

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_u<TValue, LENGTH>),       (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_u<TValue, LENGTH> const), (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
struct Value<UME::SIMD::SIMDVec_u<TValue, LENGTH> >
{
    typedef TValue Type;
};

template <typename TValue, uint32_t LENGTH_>
struct LENGTH<UME::SIMD::SIMDVec_u<TValue, LENGTH_> > {
    enum { VALUE = LENGTH_ };
};

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
getValue(UME::SIMD::SIMDVec_u<TValue, LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
value(UME::SIMD::SIMDVec_u<TValue, LENGTH> const & vector, TPosition const pos)
{

    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition, typename TValue2>
inline void
assignValue(UME::SIMD::SIMDVec_u<TValue, LENGTH> &vector, TPosition const pos, TValue2 const value)
{
    vector[pos] = value;
}

// ============================================================================
// SIMDVec_i
// ============================================================================

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_i<TValue, LENGTH>),       (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_i<TValue, LENGTH> const), (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
struct Value<UME::SIMD::SIMDVec_i<TValue, LENGTH> >
{
    typedef TValue Type;
};

template <typename TValue, uint32_t LENGTH_>
struct LENGTH<UME::SIMD::SIMDVec_i<TValue, LENGTH_> > {
    enum { VALUE = LENGTH_ };
};

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
getValue(UME::SIMD::SIMDVec_i<TValue, LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
value(UME::SIMD::SIMDVec_i<TValue, LENGTH> const & vector, TPosition const pos)
{

    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition, typename TValue2>
inline void
assignValue(UME::SIMD::SIMDVec_i<TValue, LENGTH> &vector, TPosition const pos, TValue2 const value)
{
    vector[pos] = value;
}

// ============================================================================
// SIMDVec_f
// ============================================================================

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_f<TValue, LENGTH>),       (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
SEQAN_CONCEPT_IMPL((typename UME::SIMD::SIMDVec_f<TValue, LENGTH> const), (SimdVectorConcept));

template <typename TValue, uint32_t LENGTH>
struct Value<UME::SIMD::SIMDVec_f<TValue, LENGTH> >
{
    typedef TValue Type;
};

template <typename TValue, uint32_t LENGTH_>
struct LENGTH<UME::SIMD::SIMDVec_f<TValue, LENGTH_> > {
    enum { VALUE = LENGTH_ };
};

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
getValue(UME::SIMD::SIMDVec_f<TValue, LENGTH> const & vector, TPosition const pos)
{
    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition>
inline TValue
value(UME::SIMD::SIMDVec_f<TValue, LENGTH> const & vector, TPosition const pos)
{

    return vector[pos];
}

template <typename TValue, uint32_t LENGTH, typename TPosition, typename TValue2>
inline void
assignValue(UME::SIMD::SIMDVec_f<TValue, LENGTH> &vector, TPosition const pos, TValue2 const value)
{
    vector[pos] = value;
}

} // namespace seqan

namespace UME
{
namespace SIMD
{
    template <typename TStream,
              typename TVector, typename TScalar>
    inline TStream & operator<<(TStream & stream,
               IntermediateIndex<TVector, TScalar> const & pInterIndex)
    {
        stream << static_cast<TScalar>(pInterIndex);
        return stream;
    }
}
}

namespace seqan
{

// --------------------------------------------------------------------------
// Function clearVector()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
clearVector(TSimdVector & vector)
{
    vector = 0;
}

// --------------------------------------------------------------------------
// Function createVector()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(And<Is<SimdMaskVectorConcept<TSimdVector>>,
                                Not<Is<SimdVectorConcept<TSimdVector>>>>, TSimdVector)
createVector(TValue const x)
{
    return TSimdVector(static_cast<bool>(x));
}

// --------------------------------------------------------------------------
// Function createVector()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
createVector(TValue const x)
{
    return TSimdVector(x);
}

// --------------------------------------------------------------------------
// Function fillVector()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename ...TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector & vector, TValue const... args)
{
    vector = TSimdVector(args...);
}

// --------------------------------------------------------------------------
// Function cmpEq()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpEq (TSimdVector const & a, TSimdVector const & b)
{
    return a.cmpeq(b);
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator==(TSimdVector const & a, TSimdVector const & b)
{
    return a.cmpeq(b);
}

// --------------------------------------------------------------------------
// Function operatorGt()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpGt (TSimdVector const & a, TSimdVector const & b)
{
    return a.cmpgt(b);
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator>(TSimdVector const & a, TSimdVector const & b)
{
    return a.cmpgt(b);
}

// --------------------------------------------------------------------------
// Function max()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
max(TSimdVector const & a, TSimdVector const & b)
{
    return a.max(b);
}

// --------------------------------------------------------------------------
// Function min()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
min(TSimdVector const & a, TSimdVector const & b)
{
    return a.min(b);
}

// --------------------------------------------------------------------------
// Function operator|()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator|(TSimdVector const & a, TSimdVector const & b)
{
    return a.bor(b);
}

// --------------------------------------------------------------------------
// Function operator|=()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator|=(TSimdVector & a, TSimdVector const & b)
{
    return a.bora(b);
}

// --------------------------------------------------------------------------
// Function operator&()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator&(TSimdVector const & a, TSimdVector const & b)
{
    return a.band(b);
}

// --------------------------------------------------------------------------
// Function operator&=()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator&=(TSimdVector & a, TSimdVector const & b)
{
    return a.banda(b);
}

// --------------------------------------------------------------------------
// Function operator~()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator~(TSimdVector const & a)
{
    return a.bnot();
}

// --------------------------------------------------------------------------
// Function operator+()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator+(TSimdVector const & a, TSimdVector const & b)
{
    return a.add(b);
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator-(TSimdVector const & a, TSimdVector const & b)
{
    return a.sub(b);
}

// --------------------------------------------------------------------------
// Function operator*()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator*(TSimdVector const & a, TSimdVector const & b)
{
    return a.mul(b);
}

// --------------------------------------------------------------------------
// Function operator/()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator/(TSimdVector const & a, TSimdVector const & b)
{
    return a.div(b);
}

// --------------------------------------------------------------------------
// Function andNot
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
andNot(TSimdVector const & a, TSimdVector const & b)
{
    return a.bandnot(b);
}


// --------------------------------------------------------------------------
// Function shiftRightLogical()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
shiftRightLogical(TSimdVector const & vector, const int imm)
{
    return vector.rsh(imm);
}

// --------------------------------------------------------------------------
// Function blend()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask)
{
    return a.blend(mask, b);
}

// --------------------------------------------------------------------------
// Function storeu()
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
storeu(T * memAddr, TSimdVector const & vec)
{
    vec.store(memAddr);
}

// --------------------------------------------------------------------------
// Function load()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename T>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
load(T const * memAddr)
{
    return TSimdVector(memAddr);
}

// --------------------------------------------------------------------------
// Function gather()
// --------------------------------------------------------------------------

template <typename TValue, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(IsSameType<TValue, typename Value<TSimdVector>::Type>, TSimdVector)
_gather(TValue const * memAddr, TSimdVector const & idx)
{
    using TIndexVector = typename UME::SIMD::SIMDTraits<TSimdVector>::UINT_VEC_T;

    TSimdVector a;
    a.gather(memAddr, static_cast<TIndexVector>(idx));
    return a;
}

template <typename TValue, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Not<IsSameType<TValue, typename Value<TSimdVector>::Type> >, TSimdVector)
_gather(TValue const * memAddr, TSimdVector const & idx)
{
    using TIndexVector = typename UME::SIMD::SIMDTraits<TSimdVector>::UINT_VEC_T;

    TSimdVector a;
    for (auto i = 0u; i < TIndexVector::length(); ++i)
    {
        a[i] = memAddr[idx[i]];
    }
    return a;
}

template <typename TValue, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
gather(TValue const * memAddr, TSimdVector const & idx)
{
    return _gather(memAddr, idx);
}

// --------------------------------------------------------------------------
// Function shuffleVector()
// --------------------------------------------------------------------------

template <typename TSimdVector1, typename TSimdVector2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector1> >, TSimdVector1)
shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices)
{
    return vector.swizzle(indices);
}

}

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_UMESIMD_IMPL_H_
