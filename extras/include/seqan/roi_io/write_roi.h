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

#ifndef EXTRAS_INCLUDE_SEQAN_ROI_IO_WRITE_ROI_H_
#define EXTRAS_INCLUDE_SEQAN_ROI_IO_WRITE_ROI_H_

#include <seqan/stream.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Roi_;
typedef Tag<Roi_> Roi;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeRecord()                                           [RoiRecord]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Error handling.

// Variant without I/O Context.

template <typename TStream>
int writeRecord(TStream & out, RoiRecord const & record, Roi const & /*tag*/)
{
    streamPut(out, record.ref);
    streamWriteChar(out, '\t');
    streamPut(out, record.beginPos + 1);  // 0-based to 1-based
    streamWriteChar(out, '\t');
    streamPut(out, record.endPos);
    streamWriteChar(out, '\t');
    streamPut(out, record.name);
    streamWriteChar(out, '\t');
    streamPut(out, record.len);
    streamWriteChar(out, '\t');
    streamWriteChar(out, record.strand);
    streamWriteChar(out, '\t');
    streamPut(out, record.countMax);
    streamWriteChar(out, '\t');
    for (unsigned i = 0; i < length(record.data); ++i)
    {
        streamPut(out, record.data[i]);
        streamWriteChar(out, '\t');
    }
    
    for (unsigned i = 0; i < length(record.count); ++i)
    {
        if (i > 0)
            streamWriteChar(out, ',');
        streamPut(out, record.count[i]);
    }
    streamWriteChar(out, '\n');

    return 0;
}

// Variant with I/O Context.

template <typename TStream, typename TNameStore, typename TNameStoreCache>
int writeRecord(TStream & out, RoiRecord const & record,
                RoiIOContext<TNameStore, TNameStoreCache> const & context,
                Roi const & /*tag*/)
{
    streamPut(out, nameStore(context)[record.rID]);
    streamWriteChar(out, '\t');
    streamPut(out, record.beginPos + 1);  // 0-based to 1-based
    streamWriteChar(out, '\t');
    streamPut(out, record.endPos);
    streamWriteChar(out, '\t');
    streamPut(out, record.name);
    streamWriteChar(out, '\t');
    streamPut(out, record.len);
    streamWriteChar(out, '\t');
    streamWriteChar(out, record.strand);
    streamWriteChar(out, '\t');
    streamPut(out, record.countMax);
    streamWriteChar(out, '\t');
    for (unsigned i = 0; i < length(record.data); ++i)
    {
        streamPut(out, record.data[i]);
        streamWriteChar(out, '\t');
    }
    
    for (unsigned i = 0; i < length(record.count); ++i)
    {
        if (i > 0)
            streamWriteChar(out, ',');
        streamPut(out, record.count[i]);
    }
    streamWriteChar(out, '\n');

    return 0;
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_ROI_IO_WRITE_ROI_H_
