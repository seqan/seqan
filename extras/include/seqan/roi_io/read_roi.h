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

#ifndef EXTRAS_INCLUDE_SEQAN_ROI_IO_READ_ROI_H_
#define EXTRAS_INCLUDE_SEQAN_ROI_IO_READ_ROI_H_

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
// Function readRecord()                                            [RoiRecord]
// ----------------------------------------------------------------------------

// Variant without NameStoreCache.

template <typename TStream, typename TSpec>
int readRecord(RoiRecord & record,
               RecordReader<TStream, SinglePass<TSpec> > & reader,
               Roi const & /*tag*/)
{
    CharString buffer;
    int res = 0;

    // Read reference name.
    clear(record.ref);
    if ((res = readUntilTabOrLineBreak(record.ref, reader)) != 0)
        return res;
    record.rID = RoiRecord::INVALID_REFID;

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read and parse start position.
    clear(buffer);
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.beginPos, buffer))
        return 1;
    record.beginPos -= 1;  // transform to 0-based

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read and parse end position.
    clear(buffer);
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.endPos, buffer))
        return 1;
        
    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read and parse region name.
    clear(record.name);
    if ((res = readUntilTabOrLineBreak(record.name, reader)) != 0)
        return res;

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read and parse region length.
    clear(buffer);
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.len, buffer))
        return 1;

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read strand.
    clear(buffer);
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;
    if (buffer[0] != '.' && buffer[0] != '+' && buffer[0] != '-')
        return 1;
    record.strand = buffer[0];

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Read max count.
    clear(buffer);
    if ((res = readUntilTabOrLineBreak(buffer, reader)) != 0)
        return res;
    if (!lexicalCast2(record.countMax, buffer))
        return 1;

    // Skip TAB.
    if (skipChar(reader, '\t') != 0)
        return res;

    // Buffer until the end of the line.
    seqan::CharString dataBuffer;
    if ((res = readLine(dataBuffer, reader)) != 0 && res != EOF_BEFORE_SUCCESS)
        return res;

    // Count tabs.
    int numTabs = 0;
    for (unsigned i = 0; i < length(dataBuffer); ++i)
        numTabs += (dataBuffer[i] == '\t');

    // Read data field.
    RecordReader<CharString, SinglePass<StringReader> > stringReader(dataBuffer);
    for (int tabsRead = 0; !atEnd(stringReader) && tabsRead < numTabs; ++tabsRead)
    {
        clear(buffer);
        if ((res = readGraphs(buffer, stringReader)) != 0)
            return res;
        appendValue(record.data, buffer);

        goNext(stringReader);
    }

    // TODO(holtgrew): Read additional information.

    // Individual counts.
    clear(record.count);
    clear(buffer);
    for (; !atEnd(stringReader) && value(stringReader) != '\r' && value(stringReader) != '\n'; goNext(stringReader))
    {
        if (value(stringReader) != ',')
        {
            if (!isdigit(value(stringReader)))
                return 1;  // Error parsing.
            appendValue(buffer, value(stringReader));
        }
        else
        {
            if (empty(buffer))
                continue;
            unsigned count = 0;
            if (!lexicalCast2(count, buffer))
                return 1;  // Error parsing.
            appendValue(record.count, count);
            record.countMax = std::max(record.countMax, back(record.count));
            clear(buffer);
        }
    }
    if (!empty(buffer))
    {
        unsigned count = 0;
        if (!lexicalCast2(count, buffer))
            return 1;  // Error parsing.
        appendValue(record.count, count);
        record.countMax = std::max(record.countMax, back(record.count));
    }

    return 0;
}

// Variant with NameStoreCache.

template <typename TStream, typename TSpec, typename TNameStore, typename TNameStoreCache>
int readRecord(RoiRecord & record,
               RecordReader<TStream, SinglePass<TSpec> > & reader,
               RoiIOContext<TNameStore, TNameStoreCache> & context,
               Roi const & tag)
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

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_ROI_IO_READ_ROI_H_
