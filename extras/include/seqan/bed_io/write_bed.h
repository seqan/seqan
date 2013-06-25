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
// Writing of BED to files.
// ==========================================================================

// TODO(holtgrew): Add overload with BedIOContext.

#ifndef CORE_INCLUDE_SEQAN_BED_IO_WRITE_BED_H_
#define CORE_INCLUDE_SEQAN_BED_IO_WRITE_BED_H_

namespace seqan {

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
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeRecord()                                      [Bed3 BedRecord]
// ----------------------------------------------------------------------------

/**
.Function.BedRecord#writeRecord
..cat:BED I/O
..signature:int readRecord(stream, record[, context], Bed())
..summary:Write a BED record to a file.
..description:
When $context$ is given, the first field is written from the sequence with index $record.rID$ in the context's name store.
..param.stream:@Concept.StreamConcept@ object to write to.
...type:Concept.StreamConcept
..param.record:The @Class.BedRecord@ to write.
...type:Class.BedRecord
..param.context:The optional @Class.BedIOContext@ to use.
...type:Class.BedRecord
..returns:$int$ value, $0$ on success, non-$0$ value on errors.
..include:seqan/bed_io.h
*/

// Similar to readRecord() for BedRecord, we have various overloads of _writeBedRecord().

template <typename TStream>
int _writeBedRecord(TStream & out, BedRecord<Bed3> const & record, CharString const & ref)
{
    if (streamPut(out, ref) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.beginPos + 1) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.endPos) != 0)
        return 1;

    return 0;
}

template <typename TStream>
int _writeBedRecord(TStream & out, BedRecord<Bed4> const & record, CharString const & ref)
{
    if (_writeBedRecord(out, static_cast<BedRecord<Bed3> const &>(record), ref) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.name) != 0)
        return 1;

    return 0;
}

template <typename TStream>
int _writeBedRecord(TStream & out, BedRecord<Bed5> const & record, CharString const & ref)
{
    if (_writeBedRecord(out, static_cast<BedRecord<Bed4> const &>(record), ref) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.score) != 0)
        return 1;

    return 0;
}

template <typename TStream>
int _writeBedRecord(TStream & out, BedRecord<Bed6> const & record, CharString const & ref)
{
    if (_writeBedRecord(out, static_cast<BedRecord<Bed5> const &>(record), ref) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamWriteChar(out, record.strand) != 0)
        return 1;

    return 0;
}

template <typename TStream>
int _writeBedRecord(TStream & out, BedRecord<Bed12> const & record, CharString const & ref)
{
    if (_writeBedRecord(out, static_cast<BedRecord<Bed6> const &>(record), ref) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.thickBegin + 1) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.thickEnd) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.itemRgb.red) != 0)
        return 1;
    if (streamWriteChar(out, ',') != 0)
        return 1;
    if (streamPut(out, record.itemRgb.green) != 0)
        return 1;
    if (streamWriteChar(out, ',') != 0)
        return 1;
    if (streamPut(out, record.itemRgb.blue) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.blockCount) != 0)
        return 1;

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    for (unsigned i = 0; i < length(record.blockSizes); ++i)
    {
        if (i > 0 && streamWriteChar(out, ',') != 0)
            return 1;
        if (streamPut(out, record.blockSizes[i]) != 0)
            return 1;
    }

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    for (unsigned i = 0; i < length(record.blockBegins); ++i)
    {
        if (i > 0 && streamWriteChar(out, ',') != 0)
            return 1;
        if (streamPut(out, record.blockBegins[i] - 1) != 0)
            return 1;
    }

    return 0;
}

template <typename TStream, typename TRecordSpec>
int writeRecord(TStream & out, BedRecord<TRecordSpec> const & record, Bed const & /*tag*/)
{
    if (_writeBedRecord(out, record, record.ref) != 0)
        return 1;

    if (empty(record.data))
    {
        streamWriteChar(out, '\n');
        return 0;
    }

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.data) != 0)
        return 1;

    if (streamWriteChar(out, '\n') != 0)
        return 1;

    return 0;
}

template <typename TStream, typename TRecordSpec, typename TNameStore, typename TNameStoreCache>
int writeRecord(TStream & out, BedRecord<TRecordSpec> const & record,
                BedIOContext<TNameStore, TNameStoreCache> const & context, Bed const & /*tag*/)
{

    if (record.rID == BedRecord<TRecordSpec>::INVALID_REFID)
    {
        if (_writeBedRecord(out, record, record.ref) != 0)
            return 1;
    }
    else
    {
        if (_writeBedRecord(out, record, nameStore(context)[record.rID]) != 0)
            return 1;
    }

    if (empty(record.data))
    {
        streamWriteChar(out, '\n');
        return 0;
    }

    if (streamWriteChar(out, '\t') != 0)
        return 1;

    if (streamPut(out, record.data) != 0)
        return 1;

    if (streamWriteChar(out, '\n') != 0)
        return 1;

    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BED_IO_WRITE_BED_H_
