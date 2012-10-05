// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// This file contains the MatchWriter class.
// ==========================================================================

#ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_WRITER_H_
#define SANDBOX_ESIRAGUSA_APPS_MASAI_WRITER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/find.h>

#include "tags.h"
#include "store.h"
#include "matches.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TStream, typename TDistance, typename TFormat, typename TSpec = void>
struct MatchWriter {};

template <typename TStream, typename TDistance, typename TSpec>
struct MatchWriter<TStream, TDistance, Raw, TSpec>
{
    TStream & stream;
    TFragmentStore & store;

    TReadSeqStoreSize   readsCount;

    bool                writeMatches;

    MatchWriter(TStream & stream, TFragmentStore & store, TReadSeqStoreSize readsCount, bool writeMatches) :
        stream(stream),
        store(store),
        readsCount(readsCount),
        writeMatches(writeMatches)
    {}
};

template <typename TStream, typename TDistance, typename TSpec>
struct MatchWriter<TStream, TDistance, Sam, TSpec>
{
    typedef BamIOContext<TFragmentStore::TContigNameStore>  TBamIOContext;
    typedef unsigned long                                   TWord;

    TStream & stream;
    TBamIOContext       context;
    TFragmentStore & store;

    TReadSeqStoreSize   readsCount;

    bool                writeMatches;
    bool                writeCigar;

    String<TWord>       readAligned;
    const unsigned      wordLen;

    MatchWriter(TStream & stream, TFragmentStore & store, TReadSeqStoreSize readsCount, bool writeMatches) :
        stream(stream),
        context(store.contigNameStore, store.contigNameStoreCache),
        store(store),
        readsCount(readsCount),
        writeMatches(writeMatches),
        writeCigar(true),
        wordLen(BitsPerValue<TWord>::VALUE)
    {
        resize(readAligned, (readsCount + wordLen - 1) / wordLen, (TWord)0);

        _writeHeader(store, stream, context, Sam());
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TStream, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchWriter<TStream, TDistance, Raw, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    // TODO(esiragusa):Make TMatch a template argument
    typedef Match<> TMatch;

    if (!writer.writeMatches)
        return;

    // Fill record.
    TMatch match;

    fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Write record.
    appendValue(writer.stream, match);
}

template <typename TStream, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchWriter<TStream, TDistance, Raw, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPosFwd,
                    TContigPos endPosFwd,
                    TReadId readIdFwd,
                    TErrors errorsFwd,
                    TContigPos beginPosRev,
                    TContigPos endPosRev,
                    TReadId readIdRev,
                    TErrors errorsRev)
{
    if (!writer.writeMatches)
        return;

    onMatch(writer, contigId, beginPosFwd, endPosFwd, readIdFwd, errorsFwd, false);
    onMatch(writer, contigId, beginPosRev, endPosRev, readIdRev, errorsRev, true);
}

// ============================================================================

template <typename TStream, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchWriter<TStream, TDistance, Sam, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    typedef Align<TFragmentStore::TReadSeq, ArrayGaps>  TAlign;

    TAlignedReadStoreElement        alignedRead;
    TAlignedReadStoreElement        alignedMate;
    TAlignQualityStoreElement       alignQuality;
    TAlignedReadTagStoreElement     alignedTags;

    if (!writer.writeMatches)
        return;

    // Fill aligned read.
    _fillAlignedRead(writer, alignedRead, alignQuality,
                     contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Check for secondary alignment.
    bool secondary = _checkSecondary(writer, alignedRead);

    // Align read only if cigar output is enabled.
    TAlign align;
    if (writer.writeCigar)
        _alignRead(writer, align, alignedRead, alignQuality, reverseComplemented);

    // Write aligned read.
//    _writeAlignedRead(writer.store, writer.stream, writer.context,
//                      alignedRead, alignQuality, alignedTags,
//                      alignedMate, align, secondary, Sam());

    // Fill record.
    BamAlignmentRecord record;
    _fillRecord(writer.store, record, alignedRead, alignQuality, alignedTags,
                alignedMate, align, secondary, writer.writeCigar);

    // Write record to target.
    write2(writer.stream, record, writer.context, Sam());
}

template <typename TStream, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(MatchWriter<TStream, TDistance, Sam, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPosFwd,
                    TContigPos endPosFwd,
                    TReadId readIdFwd,
                    TErrors errorsFwd,
                    TContigPos beginPosRev,
                    TContigPos endPosRev,
                    TReadId readIdRev,
                    TErrors errorsRev)
{
    typedef Align<TFragmentStore::TReadSeq, ArrayGaps>  TAlign;

    TAlignedReadStoreElement        alignedFwdMate;
    TAlignedReadStoreElement        alignedRevMate;
    TAlignQualityStoreElement       alignQualityFwdMate;
    TAlignQualityStoreElement       alignQualityRevMate;
    TAlignedReadTagStoreElement     alignedTags;

    if (!writer.writeMatches)
        return;

    // Fill aligned fwd mate.
    _fillAlignedRead(writer, alignedFwdMate, alignQualityFwdMate,
                     contigId, beginPosFwd, endPosFwd, readIdFwd, errorsFwd, false);

    // Fill aligned rev mate.
    _fillAlignedRead(writer, alignedRevMate, alignQualityRevMate,
                     contigId, beginPosRev, endPosRev, readIdRev, errorsRev, true);

    // Check for secondary alignment.
    bool secondaryFwdMate = _checkSecondary(writer, alignedFwdMate);
    bool secondaryRevMate = _checkSecondary(writer, alignedRevMate);

    // Align left mate.
    TAlign alignFwdMate;
    if (writer.writeCigar)
        _alignRead(writer, alignFwdMate, alignedFwdMate, alignQualityFwdMate, false);

    // Align right mate.
    TAlign alignRevMate;
    if (writer.writeCigar)
        _alignRead(writer, alignRevMate, alignedRevMate, alignQualityRevMate, true);

    // Fill pair information.
    alignedFwdMate.pairMatchId = alignedRevMate.readId;
    alignedRevMate.pairMatchId = alignedFwdMate.readId;

    // Write aligned left mate.
//    _writeAlignedRead(writer.store, writer.stream, writer.context,
//                      alignedFwdMate, alignQualityFwdMate, alignedTags,
//                      alignedRevMate, alignFwdMate, secondaryFwdMate, Sam());

    // Write aligned right mate.
//    _writeAlignedRead(writer.store, writer.stream, writer.context,
//                      alignedRevMate, alignQualityRevMate, alignedTags,
//                      alignedFwdMate, alignRevMate, secondaryRevMate, Sam());

    // Fill records.
    BamAlignmentRecord recordFwd;
    BamAlignmentRecord recordRev;
    _fillRecord(writer.store, recordFwd, alignedFwdMate, alignQualityFwdMate, alignedTags,
                alignedRevMate, alignFwdMate, secondaryFwdMate, writer.writeCigar);
    _fillRecord(writer.store, recordRev, alignedRevMate, alignQualityRevMate, alignedTags,
                alignedFwdMate, alignRevMate, secondaryRevMate, writer.writeCigar);

    // Write records to target.
    write2(writer.stream, recordFwd, writer.context, Sam());
    write2(writer.stream, recordRev, writer.context, Sam());
}

// ============================================================================

template <typename TStream, typename TDistance, typename TFormat, typename TSpec>
inline bool
_checkSecondary(MatchWriter<TStream, TDistance, TFormat, TSpec> & writer,
                TAlignedReadStoreElement & alignedRead)
{
    typedef MatchWriter<TStream, TDistance, TFormat, TSpec>    TMatchWriter;
    typedef typename TMatchWriter::TWord                       TWord;

    TWord mask     = (TWord)1 << (alignedRead.readId % writer.wordLen);
    bool secondary = (writer.readAligned[alignedRead.readId / writer.wordLen] & mask) != 0;

    writer.readAligned[alignedRead.readId / writer.wordLen] |= mask;

    return secondary;
}

template <typename TStream, typename TDistance, typename TFormat, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void _fillAlignedRead(MatchWriter<TStream, TDistance, TFormat, TSpec> &,
                             TAlignedReadStoreElement & alignedRead,
                             TAlignQualityStoreElement & alignQuality,
                             TContigId contigId,
                             TContigPos beginPos,
                             TContigPos endPos,
                             TReadId readId,
                             TErrors errors,
                             bool reverseComplemented)
{
    alignedRead.readId   = readId;
    alignedRead.contigId = contigId;

    if (reverseComplemented)
    {
        alignedRead.beginPos = endPos;
        alignedRead.endPos   = beginPos;
    }
    else
    {
        alignedRead.beginPos = beginPos;
        alignedRead.endPos   = endPos;
    }

    alignQuality.errors = errors;
}

template <typename TStream, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(MatchWriter<TStream, HammingDistance, TFormat, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement &,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer.store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += writer.readsCount;

    TReadSeq const & readSeq = writer.store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);
}

template <typename TStream, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(MatchWriter<TStream, EditDistance, TFormat, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement & alignQuality,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer.store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += writer.readsCount;

    TReadSeq const & readSeq = writer.store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);

    // In this case no indels are possible.
    if ((alignQuality.errors <= 1) && (length(row(align, 0)) == length(row(align, 1))))
        return;

    globalAlignment(align, Score<short, EditDistance>(),
                    (short)-alignQuality.errors, (short)alignQuality.errors,
                    NeedlemanWunsch());
}

#endif  // #ifndef SANDBOX_ESIRAGUSA_APPS_MASAI_WRITER_H_
