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

#ifndef SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_INTERFACE_H_

namespace seqan {

template <typename TSimdVector>
struct SimdMaskVectorImpl<TSimdVector, True>
{
    using Type = typename SimdVectorTraits<TSimdVector, SimdParams_<sizeof(TSimdVector), LENGTH<TSimdVector>::VALUE>>::MaskType;
};

template <typename TSimdVector>
struct SimdSwizzleVectorImpl<TSimdVector, True>
{
    typedef typename SimdVector<uint8_t, sizeof(TSimdVector)>::Type Type;
};

// ============================================================================
//
// INTERFACE FUNCTIONS
// - these should be used in the actual code, they will call one of the wrapper
//   functions defined above based on the vector type
//
// ============================================================================

// --------------------------------------------------------------------------
// Function transpose()
// --------------------------------------------------------------------------

template <int ROWS, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
transpose(TSimdVector matrix[ROWS])
{
    typedef typename Value<TSimdVector>::Type TValue;
    _transposeMatrix(matrix, SimdMatrixParams_<ROWS, LENGTH<TSimdVector>::VALUE, BitsPerValue<TValue>::VALUE>());
}

// --------------------------------------------------------------------------
// Function clearVector()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
clearVector(TSimdVector & vector)
{
    typedef typename Value<TSimdVector>::Type TValue;
    _clearVector(vector, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function createVector()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
createVector(TValue const x)
{
    typedef typename Value<TSimdVector>::Type TIVal;
    return _createVector<TSimdVector>(x, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TIVal)>());
}

// --------------------------------------------------------------------------
// Function fillVector()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename ...TValue>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
fillVector(TSimdVector & vector, TValue const... args)
{
    // On clang (<= 4.0)
    // std::make_tuple(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) reaches the
    // template recursion limit of 256 (e.g. -ftemplate-depth=256 is default)
    //
    // See same issue asked on http://stackoverflow.com/q/23374953
    // See also discussion to increase -ftemplate-depth to 1024 by default in
    // clang https://llvm.org/bugs/show_bug.cgi?id=18417
    typedef typename Value<TSimdVector>::Type TIVal;
    _fillVector(vector, std::make_tuple(args...),
                std::make_index_sequence<sizeof...(args)>{},
                SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TIVal)>());
}

// --------------------------------------------------------------------------
// Function cmpEq()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpEq (TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _cmpEq(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator==(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _cmpEq(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operatorGt()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
cmpGt (TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _cmpGt(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue), TValue>());
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, typename SimdMaskVector<TSimdVector>::Type)
operator>(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _cmpGt(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function max()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
max(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _max(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue), TValue>());
}

// --------------------------------------------------------------------------
// Function min()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
min(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _min(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue), TValue>());
}

// --------------------------------------------------------------------------
// Function operator|()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator|(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _bitwiseOr(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator|=()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator|=(TSimdVector & a, TSimdVector const & b)
{
    a = a | b;
    return a;
}

// --------------------------------------------------------------------------
// Function operator&()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator&(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _bitwiseAnd(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator&=()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector &)
operator&=(TSimdVector & a, TSimdVector const & b)
{
    a = a & b;
    return a;
}

// --------------------------------------------------------------------------
// Function operator~()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator~(TSimdVector const & a)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _bitwiseNot(a, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator+()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator+(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _add(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator-(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _sub(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator*()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator*(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _mult(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function operator/()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
operator/(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _div(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function andNot
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
andNot(TSimdVector const & a, TSimdVector const & b)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _bitwiseAndNot(a, b, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function shiftRightLogical()
// --------------------------------------------------------------------------

template <typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
shiftRightLogical(TSimdVector const & vector, const int imm)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _shiftRightLogical(vector, imm, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function blend()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename TSimdVectorMask>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
blend(TSimdVector const & a, TSimdVector const & b, TSimdVectorMask const & mask)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _blend(a, b, mask, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function storeu()
// --------------------------------------------------------------------------

template <typename T, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, void)
storeu(T * memAddr, TSimdVector const & vec)
{
    typedef typename Value<TSimdVector>::Type TValue;
    _storeu(memAddr, vec, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function load()
// --------------------------------------------------------------------------

template <typename TSimdVector, typename T>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
load(T const * memAddr)
{
    typedef typename Value<TSimdVector>::Type TValue;
    return _load<TSimdVector>(memAddr, SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TValue)>());
}

// --------------------------------------------------------------------------
// Function gather()
// --------------------------------------------------------------------------

template <typename TValue, typename TSimdVector>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector> >, TSimdVector)
gather(TValue const * memAddr, TSimdVector const & idx)
{
    typedef typename Value<TSimdVector>::Type TInnerValue;
    return _gather(memAddr, idx, std::integral_constant<size_t, sizeof(TValue)>(), SimdParams_<sizeof(TSimdVector), sizeof(TSimdVector) / sizeof(TInnerValue)>());
}

// --------------------------------------------------------------------------
// Function shuffleVector()
// --------------------------------------------------------------------------

template <typename TSimdVector1, typename TSimdVector2>
inline SEQAN_FUNC_ENABLE_IF(Is<SimdVectorConcept<TSimdVector1> >, TSimdVector1)
shuffleVector(TSimdVector1 const & vector, TSimdVector2 const & indices)
{
    typedef typename Value<TSimdVector1>::Type TValue1;
    typedef typename Value<TSimdVector2>::Type TValue2;
    return _shuffleVector(
                vector,
                indices,
                SimdParams_<sizeof(TSimdVector1), sizeof(TSimdVector1) / sizeof(TValue1)>(),
                SimdParams_<sizeof(TSimdVector2), sizeof(TSimdVector2) / sizeof(TValue2)>());
}

} // namespace seqan

#endif // SEQAN_INCLUDE_SEQAN_SIMD_SIMD_BASE_SEQAN_INTERFACE_H_
