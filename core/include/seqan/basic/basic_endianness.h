// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Basic structures and functions for endianness checks.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BASIC_BASIC_ENDIANNESS_H_
#define INCLUDE_SEQAN_BASIC_BASIC_ENDIANNESS_H_

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

template <unsigned int NUM_BITS> struct WordSize_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup ByteOrderTags Byte Order Tags
 * @brief Tags used to distinguish between endianness.
 *
 * @tag ByteOrderTags#LittleEndian
 * @headerfie <seqan/basic.h>
 * @brief Tag used for LittleEndian.
 * @signature struct ByteOrderLittleEndian_;
 *            typedef Tag<ByteOrderLittleEndian_> LittleEndian;
 *
 * @tag ByteOrderTags#BigEndian
 * @headerfie <seqan/basic.h>
 * @brief Tag used for BigEndian.
 * @signature struct ByteOrderBigEndian_;
 *            typedef Tag<ByteOrderBigEndian_> BigEndian;
 *
 * @tag ByteOrderTags#HostByteOrder
 * @headerfie <seqan/basic.h>
 * @brief Tag used to select host's byte order, which is either <tt>LittleEndian<\tt> or <tt>BigEndian<\tt>.
 */

struct ByteOrderLittleEndian_;
typedef Tag<ByteOrderLittleEndian_> LittleEndian;

struct ByteOrderBigEndian_;
typedef Tag<ByteOrderBigEndian_> BigEndian;

#if SEQAN_LITTLE_ENDIAN == 1
typedef LittleEndian HostByteOrder;
#else
typedef BigEndian HostByteOrder;
#endif

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn IsLittleEndian
 * @header <seqan/basic.h>
 * @brief Checks whether the system's byte order is little endian or not.
 *
 * @signature typename IsLittleEndian::Type;
 *
 * @return TType @link LogicalValuesTags#True @endlink if the byte order is little endian, otherwise @link LogicalValuesTags#False @endlink.
 */

#if SEQAN_LITTLE_ENDIAN == 1
struct IsLittleEndian : True{};
#else
struct IsLittleEndian : False{};
#endif

// ============================================================================
// Functions
// ============================================================================

template <typename TValue, unsigned int SIZE>
inline void _endianSwap(TValue & /*val*/, ConstUInt<SIZE>)
{
    // no-op.
}

#ifdef WINDOWS_PLATFORMS
template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<2>)
{
    val =  _byteswap_ushort(val);
}

template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<4>)
{
    val =  _byteswap_ulong(val);
}

template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<8>)
{
    val =  _byteswap_uint64(val);
}
#else  // WINDOWS_PLATFORMS
template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<2>)
{
    val = __builtin_bswap16(val);
}

template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<4>)
{
    val = __builtin_bswap32(val);
}

template <typename TValue>
inline void _endianSwap(TValue &val, ConstUInt<8>)
{
    val = __builtin_bswap64(val);
}
#endif  // WINDOWS_PLATFORMS

// ----------------------------------------------------------------------------
// Function endianSwap()
// ----------------------------------------------------------------------------

/*!
 * @fn endianSwap
 * @inlcude <seqan/basic.h>
 * @brief Swaps the byte order between the two .
 *
 * @signature   TValue ntoh(val);
 * @param   val The integral type to transform the byte order for. Must be of type @link Concept.IntegerConcept @endlink.
 *
 * @return TValue The value in host byte order.
 *
 * This function transforms the byte order of an integral type from the network to the host byte oder. The network
 * byte order is always big endian, while the host byte order depends on the used system.
 *
 * @see hton
 */

template <typename TValue, typename TFromByteOrder, typename TToByteOrder>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TValue> >, void)
endianSwap(TValue & val, TFromByteOrder /*tag*/, TToByteOrder /*tag*/)
{
    _endianSwap(val, ConstUInt<sizeof(TValue)>());
}

template <typename TValue, typename TByteOrder>
inline SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<TValue> >, void)
endianSwap(TValue & /*val*/, TByteOrder /*tag*/, TByteOrder /*tag*/)
{
    // no-op.
}

}

#endif // INCLUDE_SEQAN_BASIC_BASIC_ENDIANNESS_H_
