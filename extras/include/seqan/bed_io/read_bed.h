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
// Reading of BED from files.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BED_IO_READ_BED_H_
#define CORE_INCLUDE_SEQAN_BED_IO_READ_BED_H_

#include <seqan/stream.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bed_;
typedef Tag<Bed_> Bed;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BedRecord]
// ----------------------------------------------------------------------------

/**
.Function.BedRecord#readRecord
..cat:BED I/O
..signature:int readRecord(record, reader[, context], Bed())
..summary:Read a BED record from a file.
..description:
The type of the parameter $record$ decides which fields are interpreted.
The remainder of the line (excluding the line break) is written to $record.data$.
..description:
When $context$ is given, the $rID$ field is filled and the context's name store may be updated if a previously unknown reference occurs.
..param.record:@Class.BedRecord@ object to write to.
...type:Class.BedRecord
..param.reader:The @Spec.Single-Pass RecordReader@ to use.
...type:Spec.Single-Pass RecordReader
..param.context:The optional @Class.BedIOContext@ to use.
...type:Class.BedRecord
..returns:$int$ value, $0$ on success, non-$0$ value on errors.
..include:seqan/bed_io.h
*/

// We have a helper function _readBedRecordNoData() that has various
// overloads.  The one for Bed$N$ calls the one with Bed$N-1$.

// Helper function that reads first three fields of BED record.

template <typename TStream, typename TReaderSpec>
int _readBedRecordNoData(BedRecord<Bed3> & record,
                         RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
                         CharString & buffer)
{
    int res = 0;

    // Read CHROM.
    if ((res = readGraphs(record.ref, reader)) != 0)
        return res;

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read START.
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.beginPos, buffer))
        return 1;  // Cast failed.
    record.beginPos -= 1;  // 1-based to 0-based

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read END.
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.endPos, buffer))
        return 1;  // Cast failed.

    // Go over tab if any.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        if ((res = skipChar(reader, '\t')) != 0)
            return res;
    }

    return 0;
}

// Read first four fields without data.

template <typename TStream, typename TReaderSpec>
int _readBedRecordNoData(BedRecord<Bed4> & record,
                         RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
                         CharString & buffer)
{
    int res = 0;

    // Read first three fields.
    if ((res = _readBedRecordNoData(static_cast<BedRecord<Bed3> &>(record), reader, buffer)) != 0)
        return res;

    // Read NAME.
    if ((res = readGraphs(record.name, reader)) != 0)
        return res;

    // Go over tab if any.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        if ((res = skipChar(reader, '\t')) != 0)
            return res;
    }

    return 0;
}

// Read first five fields without data.

template <typename TStream, typename TReaderSpec>
int _readBedRecordNoData(BedRecord<Bed5> & record,
                         RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
                         CharString & buffer)
{
    int res = 0;

    // Read first three fields.
    if ((res = _readBedRecordNoData(static_cast<BedRecord<Bed4 > &>(record), reader, buffer)) != 0)
        return res;

    // Read SCORE.
    clear(buffer);
    if ((res = readGraphs(record.score, reader)) != 0)
        return res;

    // Go over tab if any.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        if ((res = skipChar(reader, '\t')) != 0)
            return res;
    }

    return 0;
}

// Read first six fields without data.

template <typename TStream, typename TReaderSpec>
int _readBedRecordNoData(BedRecord<Bed6> & record,
                         RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
                         CharString & buffer)
{
    int res = 0;

    // Read first three fields.
    if ((res = _readBedRecordNoData(static_cast<BedRecord<Bed5 > &>(record), reader, buffer)) != 0)
        return res;

    // Read STRAND.
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (buffer[0] != '.' && buffer[0] != '-' && buffer[0] != '+')
        return 1;
    record.strand = buffer[0];

    // Go over tab if any.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        if ((res = skipChar(reader, '\t')) != 0)
            return res;
    }

    return 0;
}

// Read first twelve fields without data.

template <typename TStream, typename TReaderSpec>
int _readBedRecordNoData(BedRecord<Bed12> & record,
                         RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
                         CharString & buffer)
{
    int res = 0;

    // Read first three fields.
    if ((res = _readBedRecordNoData(static_cast<BedRecord<Bed6 > &>(record), reader, buffer)) != 0)
        return res;

    // Read THICK BEGIN
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.thickBegin, buffer))
        return 1;  // Cast failed.
    record.thickBegin -= 1;  // 1-based to 0-based

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read THICK END
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.thickEnd, buffer))
        return 1;  // Cast failed.

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read ITEM RGB
    for (int i = 0; i < 3; ++i)
    {
        if (i != 0)  // Skip comma
        {
            if (atEnd(reader))
                return EOF_BEFORE_SUCCESS;
            if (value(reader) != ',')
                return 1;  // Should have been comma.
            goNext(reader);
        }

        clear(buffer);
        if ((res = readDigits(buffer, reader)) != 0)
            return res;
        if (i == 0 && !lexicalCast2(record.itemRgb.red, buffer))
            return 1;  // Could not cast value.
        else if (i == 1 && !lexicalCast2(record.itemRgb.green, buffer))
            return 1;  // Could not cast value.
        else if (i == 2 && !lexicalCast2(record.itemRgb.blue, buffer))
            return 1;  // Could not cast value.
    }

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // Read BLOCK COUNT
    clear(buffer);
    if ((res = readGraphs(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.blockCount, buffer))
        return 1;  // Cast failed.

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // READ BLOCK SIZES
    while (!atEnd(reader) && value(reader) != '\t')
    {
        if (!isdigit(value(reader)))
            return 1;  // Error, must be number.
        clear(buffer);
        if ((res = readDigits(buffer, reader)) != 0)
            return res;
        int tmp = 0;
        if (!lexicalCast2(tmp, buffer))
            return 1;  // Invalid number.
        appendValue(record.blockSizes, tmp);
        if (!atEnd(reader) && value(reader) == ',')
            goNext(reader);
    }
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;

    // Skip tab.
    if ((res = skipChar(reader, '\t')) != 0)
        return res;

    // READ BLOCK STARTS
    while (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n' && value(reader) != '\t')
    {
        if (!isdigit(value(reader)))
            return 1;  // Error, must be number.
        clear(buffer);
        if ((res = readDigits(buffer, reader)) != 0)
            return res;
        int tmp = 0;
        if (!lexicalCast2(tmp, buffer))
            return 1;  // Invalid number.
        tmp -= 1;  // 1-based to 0-based
        appendValue(record.blockBegins, tmp);
        if (!atEnd(reader) && value(reader) == ',')
            goNext(reader);
    }

    // Go over tab if any.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        if ((res = skipChar(reader, '\t')) != 0)
            return res;
    }

    return 0;
}

// The front-end function automatically calls the correct overload of
// _readBedRecordNoData through the type of record.

template <typename TSpec, typename TStream, typename TReaderSpec>
int readRecord(BedRecord<TSpec> & record,
               RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
               Bed const & /*tag*/)
{
    CharString buffer;
    clear(record);
    int res = 0;

    if ((res =_readBedRecordNoData(record, reader, buffer)) != 0)
        return res;

    // Read data if any, skipping over the line.
    if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
    {
        res = readLine(record.data, reader);
        if (res != 0 && res != EOF_BEFORE_SUCCESS)
            return res;
        else
            return 0;
    }

    // If there is no data then skip over line.
    if (!atEnd(reader))
        return skipLine(reader);
    else
        return 0;
}        

// This overload uses a BedIoContext object to translate the textual chromosome name into an index.

template <typename TSpec, typename TStream, typename TReaderSpec, typename TContextSpec, typename TContextSpec2>
int readRecord(BedRecord<TSpec> & record,
               RecordReader<TStream, SinglePass<TReaderSpec> > & reader,
               BedIOContext<TContextSpec, TContextSpec2> & context,
               Bed const & tag)
{
    int res = readRecord(record, reader, tag);
    if (res != 0)
        return res;

    // Translate chrom to rID using the context.  If there is no such sequence name in the context yet then we add it.
    unsigned idx = 0;
    if (!getIdByName(nameStore(context), record.ref, idx, nameStoreCache(context)))
    {
        idx = length(nameStore(context));
        appendName(nameStore(context), record.ref, nameStoreCache(context));
    }
    record.rID = idx;

    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BED_IO_READ_BED_H_
