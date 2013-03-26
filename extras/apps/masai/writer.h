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
// This file contains the Writer class.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_WRITER_H_
#define SEQAN_EXTRAS_MASAI_WRITER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/file/file_stream.h>

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

// ----------------------------------------------------------------------------
// Class Writer
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec = void>
struct Writer {};

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
struct Writer<TGenome, TReads, Raw, TDistance, TSpec>
{
    typedef Match<>                                         TMatch;
    typedef Stream<FileStream<WriteOnly, File<>, TMatch> >	TStream;
    
    TGenome                 & genome;
    TReads                  * reads;
    TFragmentStore          & _store;
    TStream                 _stream;
    bool                    disabled;

    Writer(TGenome & genome, bool disabled = false) :
        genome(genome),
        reads(),
        _store(genome._store),
        disabled(disabled)
    {}
};

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
struct Writer<TGenome, TReads, Sam, TDistance, TSpec>
{
    typedef BamIOContext<TFragmentStore::TContigNameStore>  TBamIOContext;
    typedef unsigned long                                   TWord;
    typedef Stream<FileStream<WriteOnly, File<> > >         TStream;

    TGenome                 & genome;
    TReads                  * reads;
    TFragmentStore          & _store;
    TStream                 _stream;
    TBamIOContext           _context;
    bool                    disabled;
    bool                    _writeCigar;
    String<TWord>           _readAligned;
    const unsigned          _wordLen;

    Writer(TGenome & genome, bool disabled = false) :
        genome(genome),
        reads(),
        _store(genome._store),
        _context(_store.contigNameStore, _store.contigNameStoreCache),
        disabled(disabled),
        _writeCigar(true),
        _wordLen(BitsPerValue<TWord>::VALUE)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _fillAlignedRead()
// ----------------------------------------------------------------------------

template <typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void _fillAlignedRead(TAlignedReadStoreElement & alignedRead,
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

// ----------------------------------------------------------------------------
// Function open()                                                     [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec, typename TString>
bool open(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & writer, TString const & fileName)
{
    if (writer.disabled)
        return true;

    if (!open(writer._stream, toCString(fileName), OPEN_RDWR | OPEN_CREATE))
        return false;

    _writeHeader(writer);

    return true;
}

// ----------------------------------------------------------------------------
// Function close()                                                    [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
bool close(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & writer)
{
    if (writer.disabled)
        return true;

    return close(writer._stream);
}

// ----------------------------------------------------------------------------
// Function setReads()                                                 [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
void setReads(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & writer, TReads & reads)
{
    _clear(writer);
    writer.reads = &reads;
    _resize(writer, reads.readsCount);
}

// ----------------------------------------------------------------------------
// Function writeAlignments()                                          [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
void writeAlignments(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & /* writer */, bool /* value */)
{}

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void writeAlignments(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer, bool value)
{
    writer._writeCigar = value;
}

// ----------------------------------------------------------------------------
// Function _resize()                                                  [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec, typename TSize>
void _resize(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & /* writer */, TSize /* count */)
{}

template <typename TGenome, typename TReads, typename TDistance, typename TSpec, typename TSize>
void _resize(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer, TSize count)
{
    resize(writer._readAligned, (count + writer._wordLen - 1) / writer._wordLen, 0);
}

// ----------------------------------------------------------------------------
// Function _clear()                                                   [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
void _clear(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & /* writer */)
{}

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void _clear(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer)
{
    clear(writer._readAligned);
}

// ----------------------------------------------------------------------------
// Function _writeHeader()                                             [Writer]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TDistance, typename TSpec>
void _writeHeader(Writer<TGenome, TReads, TFormat, TDistance, TSpec> & /* writer */)
{
}

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
void _writeHeader(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer)
{
    _writeHeader(writer._store, writer._stream, writer._context, Sam());
}

// ----------------------------------------------------------------------------
// Function onMatch()                                             [Writer<Raw>]
// ----------------------------------------------------------------------------

// Single-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Raw, TDistance, TSpec> & writer,
                    TContigId contigId,
                    TContigPos beginPos,
                    TContigPos endPos,
                    TReadId readId,
                    TErrors errors,
                    bool reverseComplemented)
{
    typedef Writer<TGenome, TReads, Raw, TDistance, TSpec>  TWriter;

    if (writer.disabled)
        return;

    // Fill record.
    typename TWriter::TMatch match;
    fill(match, contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Write record.
    streamWriteChar(writer._stream, match);
}

// Paired-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Raw, TDistance, TSpec> & writer,
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
    if (writer.disabled)
        return;

    onMatch(writer, contigId, beginPosFwd, endPosFwd, readIdFwd, errorsFwd, false);
    onMatch(writer, contigId, beginPosRev, endPosRev, readIdRev, errorsRev, true);
}

// ----------------------------------------------------------------------------
// Function onMatch()                                             [Writer<Sam>]
// ----------------------------------------------------------------------------

// Single-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer,
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

    if (writer.disabled)
        return;

    // Fill aligned read.
    _fillAlignedRead(alignedRead, alignQuality,
                     contigId, beginPos, endPos, readId, errors, reverseComplemented);

    // Check for secondary alignment.
    bool secondary = _checkSecondary(writer, alignedRead);

    // Align read only if cigar output is enabled.
    TAlign align;
    if (writer._writeCigar)
        _alignRead(writer, align, alignedRead, alignQuality, reverseComplemented);

    // Write aligned read.
//    _writeAlignedRead(writer._store, writer._stream, writer.context,
//                      alignedRead, alignQuality, alignedTags,
//                      alignedMate, align, secondary, Sam());

    // Fill record.
    BamAlignmentRecord record;
    _fillRecord(writer._store, record, alignedRead, alignQuality, alignedTags,
                alignedMate, align, secondary, writer._writeCigar);

    // Write record to target.
    write2(writer._stream, record, writer._context, Sam());
}

// Paired-End
template <typename TGenome, typename TReads, typename TDistance, typename TSpec,
          typename TContigId, typename TContigPos, typename TReadId, typename TErrors>
inline void onMatch(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer,
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

    if (writer.disabled)
        return;

    // Fill aligned fwd mate.
    _fillAlignedRead(alignedFwdMate, alignQualityFwdMate,
                     contigId, beginPosFwd, endPosFwd, readIdFwd, errorsFwd, false);

    // Fill aligned rev mate.
    _fillAlignedRead(alignedRevMate, alignQualityRevMate,
                     contigId, beginPosRev, endPosRev, readIdRev, errorsRev, true);

    // Check for secondary alignment.
    bool secondaryFwdMate = _checkSecondary(writer, alignedFwdMate);
    bool secondaryRevMate = _checkSecondary(writer, alignedRevMate);

    // Align left mate.
    TAlign alignFwdMate;
    if (writer._writeCigar)
        _alignRead(writer, alignFwdMate, alignedFwdMate, alignQualityFwdMate, false);

    // Align right mate.
    TAlign alignRevMate;
    if (writer._writeCigar)
        _alignRead(writer, alignRevMate, alignedRevMate, alignQualityRevMate, true);

    // Fill pair information.
    alignedFwdMate.pairMatchId = alignedRevMate.readId;
    alignedRevMate.pairMatchId = alignedFwdMate.readId;

    // Write aligned left mate.
//    _writeAlignedRead(writer._store, writer._stream, writer.context,
//                      alignedFwdMate, alignQualityFwdMate, alignedTags,
//                      alignedRevMate, alignFwdMate, secondaryFwdMate, Sam());

    // Write aligned right mate.
//    _writeAlignedRead(writer._store, writer._stream, writer.context,
//                      alignedRevMate, alignQualityRevMate, alignedTags,
//                      alignedFwdMate, alignRevMate, secondaryRevMate, Sam());

    // Fill records.
    BamAlignmentRecord recordFwd;
    BamAlignmentRecord recordRev;
    _fillRecord(writer._store, recordFwd, alignedFwdMate, alignQualityFwdMate, alignedTags,
                alignedRevMate, alignFwdMate, secondaryFwdMate, writer._writeCigar);
    _fillRecord(writer._store, recordRev, alignedRevMate, alignQualityRevMate, alignedTags,
                alignedFwdMate, alignRevMate, secondaryRevMate, writer._writeCigar);

    // Write records to target.
    write2(writer._stream, recordFwd, writer._context, Sam());
    write2(writer._stream, recordRev, writer._context, Sam());
}

// ----------------------------------------------------------------------------
// Function _checkSecondary()                                     [Writer<Sam>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TDistance, typename TSpec>
inline bool _checkSecondary(Writer<TGenome, TReads, Sam, TDistance, TSpec> & writer,
                            TAlignedReadStoreElement & alignedRead)
{
    typedef Writer<TGenome, TReads, Sam, TDistance, TSpec>      TWriter;
    typedef typename TWriter::TWord                             TWord;

    TWord mask     = (TWord)1 << (alignedRead.readId % writer._wordLen);
    bool secondary = (writer._readAligned[alignedRead.readId / writer._wordLen] & mask) != 0;

    writer._readAligned[alignedRead.readId / writer._wordLen] |= mask;

    return secondary;
}

// ----------------------------------------------------------------------------
// Function _alignRead()                              [Writer<HammingDistance>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, TFormat, HammingDistance, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement &,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer._store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += (writer.reads)->readsCount;

    TReadSeq const & readSeq = writer._store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);
}

// ----------------------------------------------------------------------------
// Function _alignRead()                                 [Writer<EditDistance>]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TReads, typename TFormat, typename TSpec, typename TAlign>
inline void _alignRead(Writer<TGenome, TReads, TFormat, EditDistance, TSpec> & writer,
                       TAlign & align,
                       TAlignedReadStoreElement & alignedRead,
                       TAlignQualityStoreElement & alignQuality,
                       bool reverseComplemented)
{
    typedef TFragmentStore::TReadSeq    TReadSeq;

    resize(rows(align), 2);

    assignSource(row(align, 0), infix(writer._store.contigStore[alignedRead.contigId].seq,
                                      std::min(alignedRead.beginPos, alignedRead.endPos),
                                      std::max(alignedRead.beginPos, alignedRead.endPos)));

    TReadSeqStoreSize readId = alignedRead.readId;
    if (reverseComplemented)
        readId += (writer.reads)->readsCount;

    TReadSeq const & readSeq = writer._store.readSeqStore[readId];
    assignSource(row(align, 1), readSeq);

    // In this case no indels are possible.
    if ((alignQuality.errors <= 1) && (length(row(align, 0)) == length(row(align, 1))))
        return;

    globalAlignment(align, Score<short, EditDistance>(),
                    (short)-alignQuality.errors, (short)alignQuality.errors,
                    NeedlemanWunsch());
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_WRITER_H_
