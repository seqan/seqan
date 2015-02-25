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
// Implements auxiliary functions for checksum operations.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_MISC_MISC_CHECKSUM_H
#define CORE_INCLUDE_SEQAN_MISC_MISC_CHECKSUM_H

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Private Functions
// ============================================================================

namespace impl
{

// ----------------------------------------------------------------------------
// Function _computeCrc()                                               [8-bit]
// ----------------------------------------------------------------------------

template <typename TCrcVal, typename TVal>
inline void _computeCrc(TCrcVal & crc, TVal val, ConstUInt<1>)
{
    crc = _mm_crc32_u8(crc, val);
}

// ----------------------------------------------------------------------------
// Function _computeCrc()                                              [16-bit]
// ----------------------------------------------------------------------------

template <typename TCrcVal, typename TVal>
inline void _computeCrc(TCrcVal & crc, TVal val, ConstUInt<2>)
{
    crc = _mm_crc32_u16(crc, val);
}

// ----------------------------------------------------------------------------
// Function _computeCrc()                                              [32-bit]
// ----------------------------------------------------------------------------

template <typename TCrcVal, typename TVal>
inline void _computeCrc(TCrcVal & crc, TVal val, ConstUInt<4>)
{
    crc = _mm_crc32_u32(crc, val);
}

// ----------------------------------------------------------------------------
// Function _computeCrc()                                              [64-bit]
// ----------------------------------------------------------------------------

template <typename TCrcVal, typename TVal>
inline void _computeCrc(TCrcVal & crc, TVal val, ConstUInt<8>)
{
    crc = _mm_crc32_u64(crc, val);
}

// ----------------------------------------------------------------------------
// Function _computeCrc()
// ----------------------------------------------------------------------------

template <typename TCrcVal, typename TVal, unsigned SIZE>
inline void _computeCrc(TCrcVal & crc, TVal /*val*/, ConstUInt<SIZE>)
{
    crc = 0;
}
}  // namespace impl

// ============================================================================
// Public Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeCrc()
// ----------------------------------------------------------------------------

/*!
 * @fn computeCrc
 * @headerfile <seqan/misc/misc_checksum.h>
 * @brief Performs cyclic redundancy check for the given sequence.
 *
 * @signature unsigned computeCrc(seq);
 *
 * @param seq The sequence to cmpute the crc value for.
 *
 * @return unsigned The crc value.
 *
 * Internally this function calls <tt>_mm_crc32_uX<\tt> which is only available for processors supporting
 * SSE4.2 or higher. Otherwise the returned value will always be <tt>0<\tt>.
 */

#ifdef __SSE4_2__
template <typename TSequence>
inline unsigned int
computeCrc(TSequence const & seq)
{
    typedef typename Iterator<TSequence const, Standard>::Type TIterator;
    typedef typename Value<TSequence>::Type TValue;
    typedef typename MakeUnsigned<TValue>::Type TUnsignedValue;

    unsigned int crc;
    for (TIterator it = begin(seq, Standard()); it != end(seq, Standard()); ++it)
        impl::computeCrc(crc, static_cast<TUnsignedValue>(*it), ConstUInt<sizeof(TUnsignedValue)>());
    return crc;
}
#else  //__SSE4_2_
template <typename TSequence>
inline unsigned int
computeCrc(TSequence const & /*ref*/)
{
    return 0;
}
#endif

}

#endif // CORE_INCLUDE_SEQAN_MISC_MISC_CHECKSUM_H
