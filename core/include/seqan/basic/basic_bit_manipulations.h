// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Implements generic interfaces for basic bit operations. We need the
// adaption to allow larger objects to be manipulated by bitwise operations
// to avoid an overhead when returning the result. Note, when using C++11 or
// higher, one could use rvalue references and move constructors to get rid
// of these additional interfaces.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BASIC_BASIC_BIT_MANIPULATIONS_H_
#define CORE_INCLUDE_SEQAN_BASIC_BASIC_BIT_MANIPULATIONS_H_

namespace seqan
{

// ----------------------------------------------------------------------------
// Function bitwiseAnd()
// ----------------------------------------------------------------------------

/*!
 * @fn bitwiseAnd
 * @headerfile <seqan/basic.h>
 * @brief Performs logical and for each pair of corresponding bits on the binary representation of two values of
 *        equal length.
 *
 * @signature void bitwiseAnd(res, lVal, rVal)
 *
 * @param[out] res  The result of the operation.
 * @param[in]  lVal The value on the left hand side of the operation.
 * @param[in]  rVal The value on the right hand side of the operation.
 *
 * @see bitwiseNot
 * @see bitwiseOr
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValueL, typename TValueR>
inline void
bitwiseAnd(TResult & res, TValueL const & valL, TValueR const & valR)
{
    res = valL & valR;
}

// ----------------------------------------------------------------------------
// Function bitwiseAndNot()
// ----------------------------------------------------------------------------

/*!
 * @fn bitwiseAndNot
 * @headerfile <seqan/basic.h>
 * @brief Performs logical and not for each pair of corresponding bits on the binary representation of two values of
 *        equal length.
 *
 * @signature void bitwiseAndNot(res, lVal, rVal)
 *
 * @param[out] res  The result of the operation.
 * @param[in]  lVal The value on the left hand side of the operation.
 * @param[in]  rVal The value on the right hand side of the operation.
 *
 * @see bitwiseNot
 * @see bitwiseOr
 * @see bitwiseAnd
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValueL, typename TValueR>
inline void
bitwiseAndNot(TResult & res, TValueL const & valL, TValueR const & valR)
{
    res = valL & ~valR;
}

// ----------------------------------------------------------------------------
// Function bitwiseOr()
// ----------------------------------------------------------------------------

/*!
 * @fn bitwiseOr
 * @headerfile <seqan/basic.h>
 * @brief Performs logical or for each pair of corresponding bits on the binary representation of two values of
 *        equal length.
 *
 * @signature void bitwiseOr(res, lVal, rVal)
 *
 * @param[out] res  The result of the operation.
 * @param[in]  lVal The value on the left hand side of the operation.
 * @param[in]  rVal The value on the right hand side of the operation.
 *
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValueL, typename TValueR>
inline void
bitwiseOr(TResult & res, TValueL const & valL, TValueR const & valR)
{
    res = valL | valR;
}

// ----------------------------------------------------------------------------
// Function bitwiseNot()
// ----------------------------------------------------------------------------

/*!
 * @fn bitwiseNot
 * @headerfile <seqan/basic.h>
 * @brief Performs logical not for each pair of corresponding bits on the binary representation of two values of
 *        equal length.
 *
 * @signature void bitwiseNot(res, lVal, rVal)
 *
 * @param[out] res  The result of the operation.
 * @param[in]  lVal The value on the left hand side of the operation.
 * @param[in]  rVal The value on the right hand side of the operation.
 *
 * @see bitwiseOr
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValue>
inline void
bitwiseNot(TResult & res, TValue const & val)
{
    res = ~val;
}

// ----------------------------------------------------------------------------
// Function testAllZeros()
// ----------------------------------------------------------------------------

/*!
 * @fn testAllZeros
 * @headerfile <seqan/basic.h>
 * @brief Tests whether all bits are set to 0.
 *
 * @signature bool testAllZeros(val)
 *
 * @param[in]  val The value to test the bits for.
 *
 * @return bool  True if all bits are set to 0, false otherwise.
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllOnes
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TValue>
inline bool
testAllZeros(TValue const & val)
{
    return val == static_cast<TValue>(0);
}

// ----------------------------------------------------------------------------
// Function setAllZeros()
// ----------------------------------------------------------------------------

/*!
 * @fn setAllZeros
 * @headerfile <seqan/basic.h>
 * @brief Sets all bits to 0.
 *
 * @signature void setAllZeros(val)
 *
 * @param[in, out]  val The value for which all bits are set to 0.
 *
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllZeros
 * @see testAllOnes
 * @see setAllOnes
 */

template <typename TValue>
inline void
setAllZeros(TValue & val)
{
    val = 0;
}

// ----------------------------------------------------------------------------
// Function testAllOnes()
// ----------------------------------------------------------------------------

/*!
 * @fn testAllOnes
 * @headerfile <seqan/basic.h>
 * @brief Tests whether all bits are set to 1.
 *
 * @signature bool testAllOnes(val)
 *
 * @param[in]  val The value to test the bits for.
 *
 * @return bool  True if all bits are set to 1, false otherwise.
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TValue>
inline bool
testAllOnes(TValue const & val)
{
    return val == ~static_cast<TValue>(0);
}

// ----------------------------------------------------------------------------
// Function setAllOnes()
// ----------------------------------------------------------------------------

/*!
 * @fn setAllOnes
 * @headerfile <seqan/basic.h>
 * @brief Sets all bits to 1.
 *
 * @signature void setAllOnes(val)
 *
 * @param[in, out]  val The value for which all bits are set to 1.
 *
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see bitScanReverse
 * @see testAllZeros
 * @see testAllOnes
 * @see setAllZeros
 */

template <typename TValue>
inline void
setAllOnes(TValue & val)
{
    val = ~static_cast<TValue>(0);
}

// ----------------------------------------------------------------------------
// Function bitScanReverse()
// ----------------------------------------------------------------------------

/*!
 * @fn bitScanReverse
 * @headerfile <seqan/basic.h>
 * @brief Returns the index of the last set 1 in the binary representation of the given value.
 *
 * @signature void bitScanReverse(idx, val)
 *
 * @param[out] idx The index of the last set 1.
 * @param[in]  val The value to scan. Has to be non-zero.
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanForward
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValue>
inline void
bitScanReverse(TResult & res, TValue const & in)
{
   SEQAN_ASSERT_NEQ(in, static_cast<TValue>(0));
   __asm__ ("bsr %1,%0" : "=r"(res) : "r"(in));   // Form: %0-> output; first param, %1 -> input; second param
}

// ----------------------------------------------------------------------------
// Function bitScanForward()
// ----------------------------------------------------------------------------

/*!
 * @fn bitScanForward
 * @headerfile <seqan/basic.h>
 * @brief Returns the index of the first set 1 in the binary representation of the given value.
 *
 * @signature void bitScanForward(idx, val)
 *
 * @param[out] idx The index of the first set 1.
 * @param[in]  val The value to scan. Has to be non-zero.
 *
 * @see bitwiseOr
 * @see bitwiseNot
 * @see bitwiseAnd
 * @see bitwiseAndNot
 * @see bitScanReverse
 * @see testAllOnes
 * @see testAllZeros
 * @see setAllZeros
 * @see setAllOnes
 */

template <typename TResult, typename TValue>
inline void
bitScanForward(TResult & res, TValue const & in)
{
   SEQAN_ASSERT_NEQ(in, static_cast<TValue>(0));
   __asm__ ("bsf %1,%0" : "=r"(res) : "r"(in));   // Form: %0-> output; first param, %1 -> input; second param
}

}


#endif // CORE_INCLUDE_SEQAN_BASIC_BASIC_BIT_MANIPULATIONS_H_
