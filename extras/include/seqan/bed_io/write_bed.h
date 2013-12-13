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

/*!
 * @fn BedRecord#writeRecord
 * @brief Write a BED record to a file.
 * 
 * @signature int writeRecord(stream, record[, context], Bed());
 * 
 * @param[in]     record  The BedRecord to write.
 * @param[in,out] stream  @link StreamConcept Stream @endlink object to write to.
 * @param[in,out] context The optional BedIOContext to use.
 * 
 * @return int Status code, 0 on success, other values on errors.
 * 
 * When <tt>context</tt> is given, the first field is written from the sequence with index <tt>record.rID</tt> in the
 * context's name store.
 */

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

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed3> const & record, CharString const & ref)
{
    write(target, ref);
    writeValue(target, '\t');
    // TODO(singer): Why + 1?
    // "chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0." from UCSC
    appendNumber(target, record.beginPos + 1);
    writeValue(target, '\t');
    appendNumber(target, record.endPos);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed4> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed3> const &>(record), ref);
    writeValue(target, '\t');
    write(target, record.name);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed5> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed4> const &>(record), ref);
    writeValue(target, '\t');
    write(target, record.score);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed6> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed5> const &>(record), ref);
    writeValue(target, '\t');
    writeValue(target, record.strand);
}

template <typename TTarget>
inline void
_writeBedRecord(TTarget & target, BedRecord<Bed12> const & record, CharString const & ref)
{
    _writeBedRecord(target, static_cast<BedRecord<Bed6> const &>(record), ref);
    writeValue(target, '\t');

    appendNumber(target, record.thickBegin + 1);
    writeValue(target, '\t');

    appendNumber(target, record.thickEnd);
    writeValue(target, '\t');

    appendNumber(target, record.itemRgb.red);
    writeValue(target, ',');
    appendNumber(target, record.itemRgb.green);
    writeValue(target, ',');
    appendNumber(target, record.itemRgb.blue);
    writeValue(target, '\t');

    appendNumber(target, record.blockCount);
    writeValue(target, '\t');

    for (unsigned i = 0; i < length(record.blockSizes); ++i)
    {
        if (i > 0)
            writeValue(target, ',');
        appendNumber(target, record.blockSizes[i]);
    }

    writeValue(target, '\t');

    for (unsigned i = 0; i < length(record.blockBegins); ++i)
    {
        if (i > 0)
            writeValue(target, ',');
        appendNumber(target, record.blockBegins[i] + 1);
    }
}

template <typename TTarget, typename TRecordSpec>
inline void 
writeRecord(TTarget & target, BedRecord<TRecordSpec> const & record, Bed const & /*tag*/)
{
    _writeBedRecord(target, record, record.ref);

    if (empty(record.data))
    {
        writeValue(target, '\n');
        return;
    }

    writeValue(target, '\t');
    write(target, record.data);
    writeValue(target, '\n');
}

template <typename TTarget, typename TRecordSpec, typename TNameStore, typename TNameStoreCache>
inline void 
writeRecord(TTarget & target, BedRecord<TRecordSpec> const & record,
                BedIOContext<TNameStore, TNameStoreCache> const & context, Bed const & /*tag*/)
{

    if (record.rID == BedRecord<TRecordSpec>::INVALID_REFID)
    {
        _writeBedRecord(target, record, record.ref);
    }
    else
    {
        _writeBedRecord(target, record, nameStore(context)[record.rID]);
    }

    if (empty(record.data))
    {
        writeValue(target, '\n');
        return;
    }

    writeValue(target, '\t');
    write(target, record.data);
    writeValue(target, '\n');
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BED_IO_WRITE_BED_H_
