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

// TODO(singer): get rid of internal buffers!

// Variant without NameStoreCache.

template <typename TForwardIter>
inline void readRecord(RoiRecord & record, TForwardIter & iter, Roi const & /*tag*/)
{
    typedef OrFunctor<IsTab, IsNewline> TNextEntry;

    CharString buffer;

    // Read reference name.
    clear(record.ref);
    readUntil(record.ref, iter, TNextEntry());
    record.rID = RoiRecord::INVALID_REFID;

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse start position.
    clear(buffer);
    readUntil(buffer, iter, TNextEntry());
    if (!lexicalCast(record.beginPos, buffer))
        throw BadLexicalCast(record.beginPos, buffer);
    record.beginPos -= 1;  // transform to 0-based

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse end position.
    clear(buffer);
    readUntil(buffer, iter, TNextEntry());
    if (!lexicalCast(record.endPos, buffer))
        throw BadLexicalCast(record.endPos, buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse region name.
    clear(record.name);
    readUntil(record.name, iter, TNextEntry());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read and parse region length.
    clear(buffer);
    readUntil(buffer, iter, TNextEntry());
    if (!lexicalCast(record.len, buffer))
        throw BadLexicalCast(record.len, buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read strand.
    readOne(record.strand, iter, OrFunctor<OrFunctor<EqualsChar<'.'>, EqualsChar<'+'> >, EqualsChar<'-'> >());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read max count.
    clear(buffer);
    readUntil(buffer, iter, TNextEntry());
    if (!lexicalCast(record.countMax, buffer))
        throw BadLexicalCast(record.countMax, buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Read data field.
    do
    {
        clear(buffer);
        readUntil(buffer, iter, TNextEntry());
        if (!atEnd(iter))
        {
            if(!IsNewline()(value(iter)))
                appendValue(record.data, buffer);

            skipOne(iter);
        }
        else
            break;

    }
    while (true);

    // TODO(holtgrew): Read additional information.

    // Individual counts.
    clear(record.count);
    CharString castBuffer;
    DirectionIterator<String<char>, Input>::Type castIter = begin(buffer);

    while (!atEnd(castIter))
    {
        clear(castBuffer);
        readUntil(castBuffer, castIter, OrFunctor<EqualsChar<','>, IsNewline>());
        if (!empty(castBuffer))
        {
            unsigned count = 0;
            if (!lexicalCast(count, castBuffer))
                throw BadLexicalCast(count, castBuffer);
            appendValue(record.count, count);
            record.countMax = std::max(record.countMax, back(record.count));
        }
        if (!atEnd(castIter))
            skipOne(castIter);
    }
}

// Variant with NameStoreCache.

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
inline void readRecord(RoiRecord & record,
               TForwardIter & iter,
               RoiIOContext<TNameStore, TNameStoreCache> & context,
               Roi const & tag)
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
}

}  // namespace seqan

#endif  // #ifndef EXTRAS_INCLUDE_SEQAN_ROI_IO_READ_ROI_H_
