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

/*!
 * @fn BedRecord#readRecord
 * @brief Read a BED record from a file.
 * 
 * @signature int readRecord(record, reader[, context], Bed());
 * 
 * @param[out]    record  BedRecord object to write to. Types: BedRecord
 * @param[in,out] context The optional BedIOContext to use.
 * @param[in,out] reader  The SinglePassRecordReader to use.
 * 
 * @return int Status code, 0 on success, other value on errors.
 * 
 * The type of the parameter <tt>record</tt> decides which fields are interpreted.  The remainder of the line (excluding
 * the line break) is written to <tt>record.data</tt>.
 * 
 * When <tt>context</tt> is given, the <tt>rID</tt> field is filled and the context's name store may be updated if a
 * previously unknown reference occurs.
 */

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
    // TODO(singer): Realy __int32 for a position ???
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    // TODO(singer): WHy - 1???
    // "chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0." from UCSC
    record.beginPos = lexicalCast<__int32>(buffer) - 1;
    skipOne(iter);

    // Read END.
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
    record.endPos = lexicalCast<__int32>(buffer);

    // Go over tab if any.
    if (!atEnd(iter) && value(iter) != '\r' && value(iter) != '\n')
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
    if (!atEnd(iter) && value(iter) != '\r' && value(iter) != '\n')
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
    if (!atEnd(iter) && value(iter) != '\r' && value(iter) != '\n')
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
    if (!atEnd(iter) && value(iter) != '\r' && value(iter) != '\n')
        skipOne(iter);
}

// Read first twelve fields without data.

template <typename TForwardIter>
inline void
_readBedRecordNoData(BedRecord<Bed12> & record,
                     TForwardIter & iter,
                     CharString & buffer)
{
    IsNewline isNewline;
    // Read first three fields.
    _readBedRecordNoData(static_cast<BedRecord<Bed6 > &>(record), iter, buffer);

    // Read THICK BEGIN
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.thickBegin = lexicalCast<__int32>(buffer) - 1;
    skipOne(iter);

    // Read THICK END
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.thickEnd = lexicalCast<__int32>(buffer);
    skipOne(iter);

    // Read ITEM RGB
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.itemRgb.red = lexicalCast<__int32>(buffer);
    skipOne(iter);

    clear(buffer);
    readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.itemRgb.green = lexicalCast<__int32>(buffer);
    skipOne(iter);

    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.itemRgb.blue = lexicalCast<__int32>(buffer);
    skipOne(iter);

    // Read BLOCK COUNT
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    record.blockCount = lexicalCast<__int32>(buffer);
    skipOne(iter);

    // READ BLOCK SIZES
    for (int i = 0; i < record.blockCount - 1; ++i)
    {
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
        appendValue(record.blockSizes, lexicalCast<int>(buffer));
        skipOne(iter);
    }
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
    appendValue(record.blockSizes, lexicalCast<int>(buffer));
    skipOne(iter);

    // READ BLOCK STARTS
    for (int i = 0; i < record.blockCount - 1; ++i)
    {
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<EqualsChar<','>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Bed> >());
        appendValue(record.blockBegins, lexicalCast<int>(buffer) - 1);
        skipOne(iter);
    }
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
    appendValue(record.blockBegins, lexicalCast<int>(buffer) - 1);
    // TODO(singer): Why can there be more columns than 12???
    if (isNewline(value(iter)))
    {
        skipLine(iter);
        return;
    }
    skipOne(iter);
}

// The front-end function automatically calls the correct overload of
// _readBedRecordNoData through the type of record.

template <typename TSpec, typename TForwardIter>
inline void
readRecord(BedRecord<TSpec> & record,
           TForwardIter & iter,
           Bed const & /*tag*/)
{
    CharString buffer;
    clear(record);

    _readBedRecordNoData(record, iter, buffer);

    // Read data if any, skipping over the line.
    if (!atEnd(iter) && value(iter) != '\r' && value(iter) != '\n')
    {
        readLine(record.data, iter);
        return;
    }

    // If there is no data then skip over line.
    if (!atEnd(iter))
        skipLine(iter);
}

// This overload uses a BedIoContext object to translate the textual chromosome name into an index.
template <typename TSpec, typename TForwardIter, typename TContextSpec, typename TContextSpec2>
int readRecord(BedRecord<TSpec> & record,
               TForwardIter & iter,
               BedIOContext<TContextSpec, TContextSpec2> & context,
               Bed const & tag)
{
    readRecord(record, iter, tag);

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
