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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Reading of BED from files.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BED_IO_READ_BED_H_
#define INCLUDE_SEQAN_BED_IO_READ_BED_H_

#include <seqan/stream.h>

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Bed_;
typedef Tag<Bed_> Bed;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Bed, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Bed, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Bed, T>::VALUE[1] =
{
    ".bed"     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                            [BedRecord]
// ----------------------------------------------------------------------------

// We have a helper function _readBedRecordNoData() that has various
// overloads.  The one for Bed$N$ calls the one with Bed$N-1$.

// Helper function that reads first three fields of BED record.
// NoData means the the member data (for the columns not read) is not
// filled.
template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed3> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    // Read CHROM.
    readUntil(record.ref, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    if (record.ref == "track")
    {
        skipLine(iter);
        return;
    }
    skipOne(iter);

    // Read START.
    // TODO(singer): Realy int32_t for a position ???
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    // NB: in contrast to many other text-based formats, UCSC BED uses 0-based and not 1-based coordinates.
    record.beginPos = lexicalCast<int32_t>(buffer);
    skipOne(iter);

    // Read END.
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
    record.endPos = lexicalCast<int32_t>(buffer);

    // Go over tab if any.
    if (!atEnd(iter) && IsTab()(value(iter)))
        skipOne(iter);
}

// Read first four fields without data.

template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed4> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    // Read first three fields.
    _readBedRecordNoData(static_cast<BedRecord<Bed3> &>(record), iter, buffer);

    // Read NAME.
    readUntil(record.name, iter, OrFunctor<IsTab, IsNewline>());

    // Go over tab if any.
    if (!atEnd(iter) && IsTab()(value(iter)))
        skipOne(iter);
}

// Read first five fields without data.

template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed5> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    // Read first four fields.
    _readBedRecordNoData(static_cast<BedRecord<Bed4 > &>(record), iter, buffer);

    // Read SCORE.
    readUntil(record.score, iter, OrFunctor<IsTab, IsNewline>());

    // Go over tab if any.
    if (!atEnd(iter) && IsTab()(value(iter)))
        skipOne(iter);
}

// Read first six fields without data.

template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed6> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    // Read first three fields.
    _readBedRecordNoData(static_cast<BedRecord<Bed5 > &>(record), iter, buffer);

    // Read STRAND.
    record.strand = value(iter);
    skipOne(iter, OrFunctor<OrFunctor<EqualsChar<'.'>, EqualsChar<'+'> >, EqualsChar<'-'> >());

    // Go over tab if any.
    if (!atEnd(iter) && !IsNewline()(value(iter)))
        skipOne(iter);
}

// Read first twelve fields without data.

template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed12> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    typedef AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> AssertNoNewline;

    // Read first three fields.
    _readBedRecordNoData(static_cast<BedRecord<Bed6 > &>(record), iter, buffer);

    // Read THICK BEGIN
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertNoNewline>());
    record.thickBegin = lexicalCast<int32_t>(buffer);
    skipOne(iter);

    // Read THICK END
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertNoNewline>());
    record.thickEnd = lexicalCast<int32_t>(buffer);
    skipOne(iter);

    // Read ITEM RGB
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, OrFunctor<EqualsChar<','>, AssertNoNewline> >());
    record.itemRgb.red = lexicalCast<int32_t>(buffer);
    if (value(iter) == ',')
    {
        skipOne(iter);

        clear(buffer);
        readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertNoNewline>());
        record.itemRgb.green = lexicalCast<int32_t>(buffer);
        skipOne(iter);

        clear(buffer);
        readUntil(buffer, iter, OrFunctor<IsTab, AssertNoNewline>());
        record.itemRgb.blue = lexicalCast<int32_t>(buffer);
        skipOne(iter);
    }
    else if (value(iter) == '\t' && record.itemRgb.red == 0) // allow a single '0' value for 0,0,0
    {
        skipOne(iter);
        record.itemRgb.green = 0;
        record.itemRgb.blue = 0;
    }
    else
    {
        throw ParseError("While parsing field 9 (itemRGB) in BED12 file: "
                         "value must be either 0 or a list of three integer values r,g,b.");
    }

    // Read BLOCK COUNT
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertNoNewline>());
    record.blockCount = lexicalCast<int32_t>(buffer);
    skipOne(iter);

    // READ BLOCK SIZES
    for (int i = 0; i < record.blockCount - 1; ++i)
    {
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertNoNewline>());
        appendValue(record.blockSizes, lexicalCast<int>(buffer));
        skipOne(iter);
    }
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, OrFunctor<EqualsChar<','>, AssertNoNewline> >());
    appendValue(record.blockSizes, lexicalCast<int>(buffer));
    if (value(iter) == ',') // allow trailing comma at the end of the list
        skipOne(iter);

    skipOne(iter);

    // READ BLOCK STARTS
    for (int i = 0; i < record.blockCount - 1; ++i)
    {
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertNoNewline>());
        appendValue(record.blockBegins, lexicalCast<int>(buffer));
        skipOne(iter);
    }
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, OrFunctor<EqualsChar<','>, IsNewline> >());
    appendValue(record.blockBegins, lexicalCast<int>(buffer));
    if (!atEnd(iter) && value(iter) == ',') // allow trailing comma at the end of the list
        skipOne(iter);

    if (!atEnd(iter) && value(iter) == '\t')
        skipOne(iter);

    // DO NOT skip line here, because the rest of the line (empty or not) is written
    // into record.data in the readRecord function below.
}

// The front-end function automatically calls the correct overload of
// _readBedRecordNoData through the type of record.

template <typename TSpec, typename TForwardIter>
inline void
readRecord(BedRecord<TSpec> & record,
           CharString & buffer,
           TForwardIter & iter,
           Bed const & /*tag*/)
{
    clear(record);
    _readBedRecordNoData(record, iter, buffer);
    readLine(record.data, iter);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BED_IO_READ_BED_H_
