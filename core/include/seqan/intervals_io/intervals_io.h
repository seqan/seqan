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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Reading and writing of .intervals file record plus tags.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_IO_H_
#define CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_IO_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup IntervalsFileIO Intervals File I/O
 * @brief Support for intervals file.
 *
 * The file contains format in the format <tt>CHR:POS</tt> or <tt>CHR:BEGIN-END</tt> in 1-based coordinates.
 */

// ----------------------------------------------------------------------------
// Tag Intervals
// ----------------------------------------------------------------------------

/*!
 * @tag IntervalsFileIO#Intervals
 * @headerfile <seqan/intervals_io.h>
 * @brief Tag for the intervals file format.
 */

struct Intervals_;
typedef Tag<Intervals_> Intervals;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn IntervalsFileIO#readRecord
 * @headerfile <seqan/intervals_io.h>
 * @brief Read an intervals record.
 *
 * @signature void readRecord(region, buffer, iter, tag);
 *
 * @param[out]    region The @link GenomicRegion @endlink to write interval to.
 * @param[in,out] buffer A @link CharString @endlink to use as a buffer.
 * @param[in,out] iter   A @link ForwardIteratorConcept forward iterator @endlink to read from.
 * @param[in]     tag    The tag @link IntervalsFileIO#Intervals Intervals @endlink for the format.
 */

template <typename TForwardIter>
void readRecord(GenomicRegion & record,
                CharString & buffer,
                TForwardIter & iter,
                Intervals const & /*tag*/)
{
    clear(buffer);
    // Read next line until the first whitespace (could be '\r' or '\n').
    readUntil(buffer, iter, IsWhitespace());
    // Parse buffer into GenomicRegion.
    parse(record, buffer);
    // Skip remaining line.
    skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

/*!
 * @fn IntervalsFileIO#writeRecord
 * @headerfile <seqan/intervals_io.h>
 * @brief Write out a @link GenomicRegion @endlink to an intervals file.
 *
 * @signature void writeRecord(target, buffer, region, tag);
 *
 * @param[in,out] target An @link OutputIteratorConcept output iterator @endlink to write to.
 * @param[in,out] buffer A @link CharString @endlink to use as a buffer.
 * @param[in]     region The @link GenomicRegion @endlink to write out.
 * @param[in]     tag    The tag @link IntervalsFileIO#Intervals Intervals @endlink for the format.
 */

template <typename TTarget>
void writeRecord(TTarget & target,
                 CharString & buffer,
                 GenomicRegion const & region,
                 Intervals const & /*tag*/)
{
    region.toString(buffer);
    write(target, buffer);
    write(target, "\n");
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_INTERVALS_IO_INTERVALS_IO_H_
