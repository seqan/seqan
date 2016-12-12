// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Various useful bit-twiddling routines, mostly taken from the website
// https://graphics.stanford.edu/~seander/bithacks.html
// ==========================================================================

#ifndef SEQAN_MISC_BIT_TWIDDLING_H_
#define SEQAN_MISC_BIT_TWIDDLING_H_

#ifdef STDLIB_VS

// Make intrinsics visible.  It appears that this is not necessary with VS 10
// any more, for VS 9, it must be included.
#include <intrin.h>

#ifdef __SSE4_2__
#include <nmmintrin.h>
#endif

#endif  // #ifdef STDLIB_VS

// TODO(holtgrew): Test this!

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TWord>
struct SimdVectorConcept;

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

// DeBruijn sequence for 64 bit bitScanReverse.
static const int DeBruijnMultiplyLookupBSR64[64] =
{
    0, 47,  1, 56, 48, 27,  2, 60,
   57, 49, 41, 37, 28, 16,  3, 61,
   54, 58, 35, 52, 50, 42, 21, 44,
   38, 32, 29, 23, 17, 11,  4, 62,
   46, 55, 26, 59, 40, 36, 15, 53,
   34, 51, 20, 43, 31, 22, 10, 45,
   25, 39, 14, 33, 19, 30,  9, 24,
   13, 18,  8, 12,  7,  6,  5, 63
};

// DeBruijn sequence for 64 bit bitScanForward.
static const int DeBruijnMultiplyLookupBSF64[64] =
{
    0,  1, 48,  2, 57, 49, 28,  3,
   61, 58, 50, 42, 38, 29, 17,  4,
   62, 55, 59, 36, 53, 51, 43, 22,
   45, 39, 33, 30, 24, 18, 12,  5,
   63, 47, 56, 27, 60, 41, 37, 16,
   54, 35, 52, 21, 44, 32, 23, 11,
   46, 26, 40, 15, 34, 20, 31, 10,
   25, 14, 19,  9, 13,  8,  7,  6
};

// DeBruijn sequence for 32 bit bitScanForward and bitScanReverse.
static const int DeBruijnMultiplyLookup[32] =
{
  0,   1, 28,  2, 29, 14, 24, 3,
  30, 22, 20, 15, 25, 17,  4, 8,
  31, 27, 13, 23, 21, 19, 16, 7,
  26, 12, 18,  6, 11,  5, 10, 9
};

// ----------------------------------------------------------------------------
// Tag WordSize_
// ----------------------------------------------------------------------------
// This parametrized tag is used for selecting a _popCountImpl() implementation.

template <unsigned int NUM_BITS>
struct WordSize_ {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setBitTo()
// ----------------------------------------------------------------------------

/*!
 * @fn setBitTo
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Set the bit with the given index to the given value.
 *
 * @signature void setBitTo(word, index, value);
 *
 * @param[in,out] word  The machine word (number) to set bits of (@link IntegerConcept @endlink).
 * @param[in]     index The index of the bit in the word to set (@link IntegerConcept @endlink).
 * @param[in]     value The value to set to, <tt>bool</tt>.
 */

template <typename TWord, typename TPos>
inline void
setBitTo(TWord & word, TPos index, bool value)
{
    // See http://graphics.stanford.edu/~seander/bithacks.html#ConditionalSetOrClearBitsWithoutBranching
    word = (word & ~(1u << index)) | (-value & (1u << index));
}

// ----------------------------------------------------------------------------
// Function setBit()
// ----------------------------------------------------------------------------

/*!
 * @fn setBit
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Set the bit with the given index to 1.
 *
 * @signature void setBit(word, index);
 *
 * @param[in,out] word  The word to set the bit of (@link IntegerConcept @endlink).
 * @param[in]     index The index of the bit to set (@link IntegerConcept @endlink).
 */

template <typename TWord, typename TPos>
inline void
setBit(TWord & word, TPos index)
{
    word |= (1u << index);
}

// ----------------------------------------------------------------------------
// Function clearBit()
// ----------------------------------------------------------------------------

/*!
 * @fn clearBit
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Set the bit with the given index to 0.
 *
 * @signature void clearBit(word, index);
 *
 * @param[in,out] word  The machine word to set the bit of (@link IntegerConcept @endlink).
 * @param[in]     index The index of the bit to set to 0 (@link IntegerConcept @endlink).
 */

template <typename TWord, typename TPos>
inline void
clearBit(TWord & word, TPos index)
{
    word &= ~(1u << index);
}

// ----------------------------------------------------------------------------
// Function clearAllBits()
// ----------------------------------------------------------------------------

/*!
 * @fn clearAllBits
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Set all bits to 0.
 *
 * @signature void clearAllBits(word);
 *
 * @param[in,out] word The word to clear all bits of (@link IntegerConcept @endlink).
 */

template <typename TWord>
inline void
clearBits(TWord & word)
{
    word = 0;
}

// ----------------------------------------------------------------------------
// Function isBitSet()
// ----------------------------------------------------------------------------

/*!
 * @fn isBitSet
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Returns whether the bit with the given index is set to 1.
 *
 * @signature bool isBitSet(word, index);
 *
 * @param[in] word  The word to check (@link IntegerConcept @endlink).
 * @param[in] index The index of the bit to check (@link IntegerConcept @endlink).
 *
 * @return bool Whether the bit with the given index is set in <tt>word</tt>.
 */

template <typename TWord, typename TIndex>
inline bool
isBitSet(TWord word, TIndex index)
{
    typedef typename MakeUnsigned<TWord>::Type TUnsignedWord;
    return (word & (TUnsignedWord(1) << index)) != static_cast<TWord>(0);
}

// ----------------------------------------------------------------------------
// Function hiBits()
// ----------------------------------------------------------------------------

template <typename TWord, typename TPos>
inline TWord
hiBits(TWord word, TPos index)
{
    return word & ~((TWord(1) << (BitsPerValue<TWord>::VALUE - index)) - TWord(1));
}

// ----------------------------------------------------------------------------
// Function popCount()
// ----------------------------------------------------------------------------
// The compiler-dependent implementations of _popCountImpl() follow.

/*!
 * @fn popCount
 * @headerfile <seqan/misc/bit_twiddling.h>
 * @brief Returns number of set bits in an integer.
 *
 * @signature unsigned popCount(words);
 *
 * @param[in] word The word to count the number of set bits of (@link IntegerConcept @endlink).
 *
 * @return unsigned The number of set bits in <tt>word</tt>.
 */

template <typename TWord>
inline unsigned
popCount(TWord word)
{
    return _popCountImpl(word, WordSize_<BitsPerValue<TWord>::VALUE>());
}

// Implementing this platform-independent is tricky.  There are two points to platform independentness. First, the
// choice of compiler and second the used CPU.  Currently, we do not perform any checks for the CPU and assume that
// the Intel intrinsic POPCNT is available.  The function is implemented to work on the supported compilers GCC/MINGW,
// CLANG (which has the same interface as GCC here) and Visual C++.
//
// GCC, MINGW and CLANG provide the intrinsics __builtin_popcount, __builtin_popcountl, and __builtin_popcountll for
// the types unsigned, unsigned long, and unsigned long long.  Starting with version 2008, Visual C++ provides the
// intrinsics __popcnt16, __popcnt, and __popcnt64 for 16, 32, and 64 bit words.

// MSVC >= 2008, has intrinsic
#if defined(COMPILER_MSVC) || defined(COMPILER_WINTEL)

// ----------------------------------------------------------------------------
// Function _popCountImpl()
// ----------------------------------------------------------------------------
// MSVC implementations.

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<64> const & /*tag*/)
{
#if defined(_WIN64)

#if defined(__SSE4_2__)
    // 64-bit Windows, SSE4.2 bit intrinsic available
    return _mm_popcnt_u64(static_cast<uint64_t>(word));
#elif defined(COMPILER_WINTEL)
    return _popcnt64(static_cast<uint64_t>(word));
#else
    // 64-bit Windows, 64 bit intrinsic available
    return __popcnt64(static_cast<uint64_t>(word));
#endif

#else // #if defined(_WIN64)

    // 32-bit Windows, 64 bit intrinsic not available
    return  _popCountImpl(static_cast<const uint32_t>(word), WordSize_<32>())
          + _popCountImpl(static_cast<const uint32_t>(word >> 32), WordSize_<32>());

#endif // #if defined(_WIN64)
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<32> const & /*tag*/)
{
#if defined(__SSE4_2__)
    // SSE4.2 bit intrinsic available
    return _mm_popcnt_u32(static_cast<uint32_t>(word));
#elif defined(COMPILER_WINTEL)
    return _popcnt32(static_cast<uint32_t>(word));
#else
    return __popcnt(static_cast<uint32_t>(word));
#endif
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<16> const & /*tag*/)
{
    return __popcnt16(static_cast<uint16_t>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<8> const & /*tag*/)
{
    return _popCountImpl(static_cast<const uint16_t>(word), WordSize_<16>());
}

// GCC or CLANG
#elif !(defined(COMPILER_MSVC) || defined(COMPILER_WINTEL))

// ----------------------------------------------------------------------------
// Function _popCountImpl()
// ----------------------------------------------------------------------------
// GCC or CLANG implementations.
// SSE4.2 popcnt is emitted when compiling with -mpopcnt or -march=corei7

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<64> const & /*tag*/)
{
    return __builtin_popcountll(static_cast<uint64_t>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<32> const & /*tag*/)
{
    return __builtin_popcount(static_cast<uint32_t>(word));
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<16> const & /*tag*/)
{
    return _popCountImpl(static_cast<uint32_t>(word), WordSize_<32>());
}

template <typename TWord>
inline unsigned
_popCountImpl(TWord word, WordSize_<8> const & /*tag*/)
{
    return _popCountImpl(static_cast<uint32_t>(word), WordSize_<32>());
}

#endif // #if !(defined(COMPILER_MSVC) || defined(COMPILER_WINTEL))

// ----------------------------------------------------------------------------
// Function printBits()
// ----------------------------------------------------------------------------

//template <typename TValue>
//inline void printBits(TValue word)
//{
//    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
//    TValue one = 1;
//    for (TValue i = 0; i < bitsPerValue; ++i)
//        std::cout << ((word >> i) & one);
//    std::cout << std::endl;
//}

//template <typename TValue, typename TSize>
//inline std::ostream & printBits(std::ostream & stream, TValue word, TSize blockSize)
//{
//    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
//    bool temp;
//    for (int i = bitsPerValue - 1; i >= 0; --i)
//    {
//        temp = (word >> i) & 1;
//        stream << temp;
//        if ((bitsPerValue - i) % blockSize == 0)
//            stream << " ";
//    }
//    return stream;
//}

// ----------------------------------------------------------------------------
// Function testAllZeros()
// ----------------------------------------------------------------------------

/*!
 * @fn testAllZeros
 * @headerfile <seqan/misc.h>
 * @brief Tests whether all bits of the given value are set to <b>0</b>.
 *
 * @signature bool testAllZeros(val)
 *
 * @param[in]  val The value to check the bits for. Must be of type @link IntegerConcept @endlink.
 *
 * @return bool  <tt>true</tt> if all bits are set to <b>0</b>, <tt>false</tt> otherwise.
 *
 * @see testAllOnes
 */

template <typename TWord>
inline bool testAllZeros(TWord const & val)
{
    return val == 0;
}

// ----------------------------------------------------------------------------
// Function testAllOnes()
// ----------------------------------------------------------------------------

/*!
 * @fn testAllOnes
 * @headerfile <seqan/misc.h>
 * @brief Tests whether all bits of the given value are set to <b>1</b>.
 *
 * @signature bool testAllOnes(val)
 *
 * @param[in]  val The value to check the bits for. Must be of type @link IntegerConcept @endlink.
 *
 * @return bool <tt>true</tt> if all bits are set to <b>1</b>, <tt>false</tt> otherwise.
 *
 * @see testAllZeros
 */

template <typename TWord>
inline bool _testAllOnes(TWord const & val, False)
{
    return val == ~static_cast<TWord>(0);
}

template <typename TWord>
inline bool testAllOnes(TWord const & val)
{
    return _testAllOnes(val, typename Is<SimdVectorConcept<TWord> >::Type());
}

// ----------------------------------------------------------------------------
// Function _bitScanReverseGeneric()                     [Platform independent]
// ----------------------------------------------------------------------------

// bitScanReverse for 32 bit integers using DeBruijn sequence by Eric Cole, January 8, 2006.
template <typename TWord>
inline TWord
_bitScanReverseGeneric(TWord word, WordSize_<32>)
{
    return DeBruijnMultiplyLookup[static_cast<uint32_t>(word * 0x077CB531U) >> 27];
}

// bitScanReverse for 64 bit integers using DeBruijn sequence by Kim Walisch and Mark Dickinson.
template <typename TWord>
inline TWord
_bitScanReverseGeneric(TWord word, WordSize_<64>)
{
    word |= word >> 1; word |= word >> 2; word |= word >> 4; word |= word >> 8; word |= word >> 16; word |= word >> 32;
    return DeBruijnMultiplyLookupBSR64[static_cast<uint64_t>(word * 0x03f79d71b4cb0a89ULL) >> 58];
}

// ----------------------------------------------------------------------------
// Function _bitScanReverse()                            [Platform independent]
// ----------------------------------------------------------------------------

template <typename TWord, unsigned int NUM_BITS>
inline TWord
_bitScanReverse(TWord word, WordSize_<NUM_BITS>)
{
    return _bitScanReverseGeneric(word, WordSize_<NUM_BITS>());
}

// ----------------------------------------------------------------------------
// Function _bitScanForwardGeneric()                     [Platform independent]
// ----------------------------------------------------------------------------

// bitScanForward implementations for 64 and 32 bit values using DeBruijn sequence by Martin L�uter, Charles E. Leiserson,
// Harald Prokop and Keith H. Randall; "Using de Bruijn Sequences to Index a 1 in a Computer Word"; (1997)

// Note, the cast of word to a signed integer is necessary to fix compiler warning C4146 on Windows platforms.
template <typename TWord>
inline TWord
_bitScanForwardGeneric(TWord word, WordSize_<32>)
{
    return DeBruijnMultiplyLookup[static_cast<uint32_t>(((word & -static_cast<int32_t>(word)) * 0x077CB531U)) >> 27];
}

template <typename TWord>
inline TWord
_bitScanForwardGeneric(TWord word, WordSize_<64>)
{
    return DeBruijnMultiplyLookupBSF64[static_cast<uint64_t>((word & -static_cast<int64_t>(word)) * 0x03f79d71b4cb0a89ULL) >> 58];
}

// ----------------------------------------------------------------------------
// Function _bitScanForward()                            [Platform independent]
// ----------------------------------------------------------------------------

template <typename TWord, unsigned int NUM_BITS>
inline TWord
_bitScanForward(TWord word, WordSize_<NUM_BITS>)
{
    return _bitScanForwardGeneric(word, WordSize_<NUM_BITS>());
}

// NOTE(marehr): The Intel compiler on windows behaves the same as the visual
// studio compiler and on *nix the same as gcc. Thus, the __builtin_clz is only
// available on *nix.
#if defined(COMPILER_GCC) || defined(COMPILER_CLANG) || defined(COMPILER_LINTEL)

template <typename TWord>
inline TWord
_bitScanReverse(TWord word, WordSize_<64>)
{
    return 63 - __builtin_clzll(static_cast<unsigned long long>(word));
}

template <typename TWord>
inline TWord
_bitScanReverse(TWord word, WordSize_<32>)
{
    return 31 - __builtin_clz(static_cast<unsigned int>(word));
}


template <typename TWord>
inline TWord
_bitScanForward(TWord word, WordSize_<64>)
{
    return __builtin_ctzll(static_cast<unsigned long long>(word));
}

template <typename TWord>
inline TWord
_bitScanForward(TWord word, WordSize_<32>)
{
    return __builtin_ctz(static_cast<unsigned int>(word));
}

#elif defined(STDLIB_VS) // #if !(defined(COMPILER_GCC) || defined(COMPILER_CLANG) || defined(COMPILER_LINTEL)) && defined(STDLIB_VS)

#if (SEQAN_IS_64_BIT)

template <typename TWord>
inline TWord
_bitScanReverse(TWord word, WordSize_<64>)
{
    unsigned long index;
    _BitScanReverse64(&index, static_cast<uint64_t>(word));
    return index;
}

template <typename TWord>
inline TWord
_bitScanForward(TWord word, WordSize_<64>)
{
    unsigned long index;
    _BitScanForward64(&index, static_cast<uint64_t>(word));
    return index;
}
#else

template <typename TWord>
inline TWord
_bitScanReverse(TWord word, WordSize_<64>)
{
    unsigned long index;
    unsigned long hi = word >> 32;
    if (hi == 0u)
    {
        _BitScanReverse(&index, word);
        return index;
    }
    _BitScanReverse(&index, hi);
    return index + 32;
}

template <typename TWord>
inline TWord
_bitScanForward(TWord word, WordSize_<64>)
{
    unsigned long index;
    unsigned long lo = word & ~static_cast<unsigned long>(0);
    if (lo == 0u)
    {
        _BitScanForward(&index, word >> 32);
        return index + 32;
    }
    _BitScanForward(&index, lo);
    return index;
}
#endif  // if (SEQAN_IS_64_BIT)

template <typename TWord>
inline TWord
_bitScanReverse(TWord word, WordSize_<32>)
{
    unsigned long index;
    _BitScanReverse(&index, static_cast<unsigned long>(word));
    return index;
}

template <typename TWord>
inline TWord
_bitScanForward(TWord word, WordSize_<32>)
{
    unsigned long index;
    _BitScanForward(&index, static_cast<unsigned long>(word));
    return index;
}
#endif  // #if !(defined(COMPILER_GCC) || defined(COMPILER_CLANG) || defined(COMPILER_LINTEL)) && defined(STDLIB_VS)

// ----------------------------------------------------------------------------
// Function bitScanReverse()
// ----------------------------------------------------------------------------

/*!
 * @fn bitScanReverse
 * @headerfile <seqan/misc.h>
 * @brief Returns the index of the last set bit in the binary representation of the given value.
 * @note If <tt>val</tt> is 0 the return value is undefined.
 *
 * @signature TWord bitScanReverse(val)
 *
 * @param[in]  val The value to scan. Has to be non-zero.
 *
 * @return TWord The index of the last set bit in <tt>val</tt>, where <tt>TWord</tt> is the value of <tt>val</tt>.
 *
 * @see bitScanForward
 */

template <typename TWord>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TWord> >, TWord)
bitScanReverse(TWord word)
{
   SEQAN_ASSERT_NEQ(word, static_cast<TWord>(0));

   return _bitScanReverse(word, WordSize_<(BitsPerValue<TWord>::VALUE <= 32) ? 32 : BitsPerValue<TWord>::VALUE>());
}

// ----------------------------------------------------------------------------
// Function bitScanForward()
// ----------------------------------------------------------------------------

/*!
 * @fn bitScanForward
 * @headerfile <seqan/misc.h>
 * @brief Returns the index of the first set bit in the binary representation of the given value.
 * @note If <tt>val</tt> is 0 the return value is undefined.
 *
 * @signature TWord bitScanForward(val)
 *
 * @param[in]  val The value to scan. Has to be non-zero.
 *
 * @return TWord The index of the first set bit in <tt>val</tt>, where <tt>TWord</tt> is the value of <tt>val</tt>.
 *
 * @see bitScanReverse
 */

template <typename TWord>
inline SEQAN_FUNC_ENABLE_IF( Is<IntegerConcept<TWord> >, TWord)
bitScanForward(TWord word)
{
   SEQAN_ASSERT_NEQ(word, static_cast<TWord>(0));
   return _bitScanForward(word, WordSize_<(BitsPerValue<TWord>::VALUE <= 32) ? 32 : BitsPerValue<TWord>::VALUE>());
}

}  // namespace seqan

#endif // #ifndef SEQAN_MISC_BIT_TWIDDLING_H_
