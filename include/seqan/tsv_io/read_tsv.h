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
// Author: David Weese <dave.weese@gmail.com>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_TSV_READ_TSV_H_
#define SEQAN_INCLUDE_SEQAN_TSV_READ_TSV_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Tsv
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Tsv
 * @headerfile <seqan/tsv_io.h>
 * @brief Variant callinf format file.
 *
 * @signature typedef Tag<Tsv_> Tsv;
 */
struct Tsv_;
typedef Tag<Tsv_> Tsv;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readHeader()                                            [TsvHeader]
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline void
readHeader(TsvHeader & header,
           TsvIOContext & context,
           TForwardIter & iter,
           Tsv const & /*tag*/)
{
    clear(header);
    CharString &buffer = context.buffer;

    clear(buffer);
    readLine(buffer, iter);
    
    strSplit(header, buffer, IsTab(), false);
}

// ----------------------------------------------------------------------------
// Function readRecord()                                            [TsvRecord]
// ----------------------------------------------------------------------------
// Read one line and split it at tabs into values.

template <typename TForwardIter>
inline void
readRecord(TsvRecord & record,
           TsvIOContext & context,
           TForwardIter & iter,
           Tsv const & /*tag*/)
{
    clear(record);
    CharString &buffer = context.buffer;

    clear(buffer);
    do
    {
        readLine(buffer, iter);
    } while (empty(buffer) && !atEnd(iter));    // skip empty lines

    strSplit(record, buffer, IsTab(), false);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_TSV_READ_TSV_H_
