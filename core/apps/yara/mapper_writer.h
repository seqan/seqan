// ==========================================================================
//                      Yara - Yet Another Read Aligner
// ==========================================================================
// Copyright (c) 2011-2014, Enrico Siragusa, FU Berlin
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
//     * Neither the name of Enrico Siragusa or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ENRICO SIRAGUSA OR THE FU BERLIN BE LIABLE
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

#ifndef APP_YARA_MAPPER_WRITER_H_
#define APP_YARA_MAPPER_WRITER_H_

using namespace seqan;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class MatchesWriter
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
struct MatchesWriter
{
    typedef typename Traits::TContigs          TContigs;
    typedef typename Traits::TReads            TReads;
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TMatches          TMatches;
    typedef typename Traits::TCigarSet         TCigarSet;
    typedef typename Traits::TOutputStream     TOutputStream;
    typedef typename Traits::TOutputContext    TOutputContext;
    typedef typename Traits::TReadsContext     TReadsContext;

    // Thread-private data.
    BamAlignmentRecord      record;
    CharString              recordBuffer;
    CharString              xa;

    // Shared-memory read-write data.
    TOutputStream &         outputStream;
    TOutputContext &        outputCtx;

    // Shared-memory read-only data.
    TMatchesSet const &     matchesSet;
    TMatches const &        primaryMatches;
    TCigarSet const &       cigarSet;
    TReadsContext const &   ctx;
    TContigs const &        contigs;
    TReads const &          reads;
    Options const &         options;

    MatchesWriter(TOutputStream & outputStream,
                  TOutputContext & outputCtx,
                  TMatchesSet const & matchesSet,
                  TMatches const & primaryMatches,
                  TCigarSet const & cigarSet,
                  TReadsContext const & ctx,
                  TContigs const & contigs,
                  TReads const & reads,
                  Options const & options) :
        outputStream(outputStream),
        outputCtx(outputCtx),
        matchesSet(matchesSet),
        primaryMatches(primaryMatches),
        cigarSet(cigarSet),
        ctx(ctx),
        contigs(contigs),
        reads(reads),
        options(options)
    {
        // Process all matches.
        iterate(primaryMatches, *this, Standard(), typename Traits::TThreading());
    }

    template <typename TIterator>
    void operator() (TIterator const & it)
    {
        _writeMatchesImpl(*this, it);
    }

    ~MatchesWriter()
    {
        if (!empty(recordBuffer))
            _writeRecordBufferImpl(*this, typename Traits::TThreading());
    }
};

// ----------------------------------------------------------------------------
// Class QualityExtractor
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this when new tokenization gets into develop.

template <typename TValue>
struct QualityExtractor : public std::unary_function<TValue, char>
{
    inline char operator()(TValue const & x) const
    {
        return '!' + static_cast<char>(getQualityValue(x));
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function fillHeader()
// ----------------------------------------------------------------------------

template <typename TContigSeqs, typename TOptions, typename TContigNames>
inline void fillHeader(BamHeader & header, TOptions const & options, TContigSeqs const & seqs, TContigNames const & names)
{
    typedef typename Iterator<TContigSeqs const, Standard>::Type    TContigSeqsIter;
    typedef typename Iterator<TContigNames const, Standard>::Type   TContigNamesIter;

    typedef BamHeader::TSequenceInfo                                TSequenceInfo;
    typedef BamHeaderRecord::TTag                                   TTag;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.4"));
    appendValue(firstRecord.tags, TTag("SO", "unsorted"));
    appendValue(header.records, firstRecord);

    // Fill sequence info header line.
    TContigSeqsIter sit = begin(seqs, Standard());
    TContigSeqsIter sitEnd = end(seqs, Standard());
    TContigNamesIter nit = begin(names, Standard());
    TContigNamesIter nitEnd = end(names, Standard());

    for (; sit != sitEnd && nit != nitEnd; ++sit, ++nit)
        appendValue(header.sequenceInfos, TSequenceInfo(*nit, length(*sit)));

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "SeqAn"));
    appendValue(pgRecord.tags, TTag("PN", "Yara"));
    appendValue(pgRecord.tags, TTag("VN", options.version));
    appendValue(pgRecord.tags, TTag("CL", options.commandLine));
    appendValue(header.records, pgRecord);
}

// ----------------------------------------------------------------------------
// Function setQual()
// ----------------------------------------------------------------------------

template <typename TString>
inline void setQual(BamAlignmentRecord & record, TString const & string)
{
    typedef typename Value<TString>::Type                               TAlphabet;
    typedef QualityExtractor<TAlphabet>                                 TQualityExtractor;
    typedef ModifiedString<TString const, ModView<TQualityExtractor> >  TQualities;

    TQualities qual(string);
    record.qual = qual;
}

// ----------------------------------------------------------------------------
// Function append*()
// ----------------------------------------------------------------------------

template <typename TErrors>
inline void appendErrors(BamAlignmentRecord & record, TErrors errors)
{
    appendTagValue(record.tags, "NM", errors, 'i');
}

template <typename TCount>
inline void appendCooptimalCount(BamAlignmentRecord & record, TCount count)
{
    appendTagValue(record.tags, "X0", count, 'i');
}

template <typename TCount>
inline void appendSuboptimalCount(BamAlignmentRecord & record, TCount count)
{
    appendTagValue(record.tags, "X1", count, 'i');
}

inline void appendType(BamAlignmentRecord & record, bool unique)
{
    appendTagValue(record.tags, "XT", unique ? 'U' : 'R', 'A');
}

template <typename TPos>
inline void appendExtraPosition(BamAlignmentRecord & record, TPos pos)
{
    appendTagValue(record.tags, "XP", pos + 1, 'i');
}

template <typename TString>
inline void appendAlignments(BamAlignmentRecord & record, TString const & xa)
{
    if (!empty(xa)) appendTagValue(record.tags, "XA", xa, 'Z');
}

// ----------------------------------------------------------------------------
// Function _writeMatchesImpl()
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TMatchIt>
inline void _writeMatchesImpl(MatchesWriter<TSpec, Traits> & me, TMatchIt const & it)
{
    typedef typename Value<TMatchIt const>::Type    TMatch;

    TMatch const & primary = value(it);

    if (isValid(primary))
        _writeMappedRead(me, position(it, me.primaryMatches), primary);
    else
        _writeUnmappedRead(me, position(it, me.primaryMatches));
}

// ----------------------------------------------------------------------------
// Function _writeUnmappedRead()
// ----------------------------------------------------------------------------
// Writes one unmapped read.

template <typename TSpec, typename Traits, typename TReadId>
inline void _writeUnmappedRead(MatchesWriter<TSpec, Traits> & me, TReadId readId)
{
    clear(me.record);
    _fillReadName(me, readId);
    _fillReadSeqQual(me, readId);
    _fillMapq(me, 0u);
    _fillMateInfo(me, readId);
    me.record.flag |= BAM_FLAG_UNMAPPED;
    _writeRecord(me);
}

// ----------------------------------------------------------------------------
// Function _writeMappedRead()
// ----------------------------------------------------------------------------
// Writes one block of matches.

template <typename TSpec, typename Traits, typename TReadId, typename TMatch>
inline void _writeMappedRead(MatchesWriter<TSpec, Traits> & me, TReadId readId, TMatch const & primary)
{
    _writeMappedReadImpl(me, readId, primary, typename Traits::TSequencing());
}

template <typename TSpec, typename Traits, typename TReadId, typename TMatch>
inline void _writeMappedReadImpl(MatchesWriter<TSpec, Traits> & me, TReadId readId, TMatch const & primary, SingleEnd)
{
    typedef typename Traits::TMatches           TMatches;
    typedef typename Size<TMatches>::Type       TSize;

    clear(me.record);
    _fillReadName(me, getReadSeqId(primary, me.reads.seqs));
    _fillReadSeqQual(me, getReadSeqId(primary, me.reads.seqs));
    _fillReadPosition(me, primary);
    _fillReadAlignment(me, primary);
    _fillMateInfo(me, readId);

    TMatches const & matches = me.matchesSet[readId];
    TSize bestCount = countBestMatches(matches);
    _fillReadInfo(me, matches, bestCount);

    if (!me.options.outputSecondary)
        _fillXa(me, matches, bestCount, 0u);

    _writeRecord(me);

    if (me.options.outputSecondary)
        _writeSecondary(me, matches, bestCount, 0u);
}

template <typename TSpec, typename Traits, typename TReadId, typename TMatch>
inline void _writeMappedReadImpl(MatchesWriter<TSpec, Traits> & me, TReadId readId, TMatch const & primary, PairedEnd)
{
    typedef typename Traits::TMatches                           TMatches;
    typedef typename Size<TMatches>::Type                       TSize;
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    clear(me.record);
    _fillReadName(me, getReadSeqId(primary, me.reads.seqs));
    _fillReadSeqQual(me, getReadSeqId(primary, me.reads.seqs));
    _fillReadPosition(me, primary);
    _fillReadAlignment(me, primary);
    _fillMateInfo(me, readId);

    if (isPaired(me.ctx, readId))
    {
        TReadId mateId = getMateId(me.reads.seqs, readId);
        TMatch const & mate = me.primaryMatches[mateId];
        _fillMatePosition(me, primary, mate);
    }

    TMatches const & matches = me.matchesSet[readId];
    TSize bestCount = countBestMatches(matches);
    _fillReadInfo(me, matches, bestCount);

    // Find the primary match in the list of matches.
    TIter it = findMatch(matches, primary);
    TSize primaryPos = position(it, matches);

    if (!me.options.outputSecondary)
        _fillXa(me, matches, bestCount, primaryPos);

    _writeRecord(me);

    if (me.options.outputSecondary)
        _writeSecondary(me, matches, bestCount, primaryPos);
}

// ----------------------------------------------------------------------------
// Function _writeSecondary()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _writeSecondary(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount, TPos primaryPos)
{
    _writeSecondaryImpl(me, matches, bestCount, primaryPos, typename Traits::TStrategy());
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _writeSecondaryImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount /* bestCount */, TPos primaryPos, All)
{
    _writeSecondary(me, prefix(matches, primaryPos));
    _writeSecondary(me, suffix(matches, primaryPos + 1));
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _writeSecondaryImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount, TPos primaryPos, Strata)
{
    if (primaryPos < bestCount)
    {
        TMatches const & cooptimal = prefix(matches, bestCount);
        _writeSecondary(me, prefix(cooptimal, primaryPos));
        _writeSecondary(me, suffix(cooptimal, primaryPos + 1));
    }
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _writeSecondary(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    TIter itEnd = end(matches, Standard());
    for (TIter it = begin(matches, Standard()); it != itEnd; ++it)
    {
        clear(me.record);
        _fillReadName(me, getReadSeqId(value(it), me.reads.seqs));
        _fillReadPosition(me, value(it));
        appendExtraPosition(me.record, getContigEnd(value(it)));
        me.record.flag |= BAM_FLAG_SECONDARY;
        _writeRecord(me);
    }
}

// ----------------------------------------------------------------------------
// Function _fillReadName()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TReadSeqId>
inline void _fillReadName(MatchesWriter<TSpec, Traits> & me, TReadSeqId readSeqId)
{
    me.record.qName = me.reads.names[getReadId(me.reads.seqs, readSeqId)];
}

// ----------------------------------------------------------------------------
// Function _fillReadSeqQual()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TReadSeqId>
inline void _fillReadSeqQual(MatchesWriter<TSpec, Traits> & me, TReadSeqId readSeqId)
{
    me.record.seq = me.reads.seqs[readSeqId];
    setQual(me.record, me.reads.seqs[readSeqId]);
}

// ----------------------------------------------------------------------------
// Function _fillReadPosition()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatch>
inline void _fillReadPosition(MatchesWriter<TSpec, Traits> & me, TMatch const & match)
{
    if (onReverseStrand(match))
        me.record.flag |= BAM_FLAG_RC;

    me.record.rID = getContigId(match);
    me.record.beginPos = getContigBegin(match);
    appendErrors(me.record, getErrors(match));
}

// ----------------------------------------------------------------------------
// Function _fillReadAlignment()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatch>
inline void _fillReadAlignment(MatchesWriter<TSpec, Traits> & me, TMatch const & match)
{
    me.record.cigar = me.cigarSet[getReadId(match)];
}

// --------------------------------------------------------------------------
// Function _fillMateInfo()
// --------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TReadId>
inline void _fillMateInfo(MatchesWriter<TSpec, Traits> & me, TReadId readId)
{
    _fillMateInfoImpl(me, readId, typename Traits::TSequencing());
}

template <typename TSpec, typename Traits, typename TReadId>
inline void _fillMateInfoImpl(MatchesWriter<TSpec, Traits> & /* me */, TReadId /* readId */, SingleEnd) {}

template <typename TSpec, typename Traits, typename TReadId>
inline void _fillMateInfoImpl(MatchesWriter<TSpec, Traits> & me, TReadId readId, PairedEnd)
{
    me.record.flag |= BAM_FLAG_NEXT_UNMAPPED;
    me.record.flag |= BAM_FLAG_MULTIPLE;

    if (isFirstMate(me.reads.seqs, readId))
        me.record.flag |= BAM_FLAG_FIRST;
    else
        me.record.flag |= BAM_FLAG_LAST;
}

// --------------------------------------------------------------------------
// Function _fillMatePosition()
// --------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatch>
inline void _fillMatePosition(MatchesWriter<TSpec, Traits> & me, TMatch const & match, TMatch const & mate)
{
    me.record.flag &= ~BAM_FLAG_NEXT_UNMAPPED;
    me.record.flag |= BAM_FLAG_ALL_PROPER;

    if (onReverseStrand(mate))
        me.record.flag |= BAM_FLAG_NEXT_RC;

    me.record.rNextId = getContigId(mate);
    me.record.pNext = getContigBegin(mate);

    if (getContigId(match) == getContigId(mate))
    {
        if (getContigBegin(match) < getContigBegin(mate))
            me.record.tLen = getContigEnd(mate) - getContigBegin(match);
        else
            me.record.tLen = getContigBegin(mate) - getContigEnd(match);
    }
}

// ----------------------------------------------------------------------------
// Function _fillMapq()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TCount>
inline void _fillMapq(MatchesWriter<TSpec, Traits> & me, TCount count)
{
    static const unsigned char MAPQ[] = { 0, 40, 3, 2, 1, 1, 1, 1, 1, 1, 0 };

    me.record.mapQ = MAPQ[_min(count, 10u)];
}

// ----------------------------------------------------------------------------
// Function _fillReadInfo()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches, typename TCount>
inline void _fillReadInfo(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount)
{
    _fillReadInfoImpl(me, matches, bestCount, typename Traits::TStrategy());
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount>
inline void _fillReadInfoImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount, All)
{
    _fillMapq(me, bestCount);
    appendCooptimalCount(me.record, bestCount);
    appendSuboptimalCount(me.record, length(matches) - bestCount);
    appendType(me.record, bestCount == 1);
    // Set number of secondary alignments and hit index.
//    appendTagValue(me.record.tags, "NH", 1, 'i');
//    appendTagValue(me.record.tags, "HI", 1, 'i');
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount>
inline void _fillReadInfoImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & /* matches */, TCount bestCount, Strata)
{
    _fillMapq(me, bestCount);
    appendCooptimalCount(me.record, bestCount);
//    appendSuboptimalCount(me.record, 0);
    appendType(me.record, bestCount == 1);
}

// ----------------------------------------------------------------------------
// Function _fillXa()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _fillXa(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount, TPos primaryPos)
{
    _fillXaImpl(me, matches, bestCount, primaryPos, typename Traits::TStrategy());
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _fillXaImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount /* bestCount */, TPos primaryPos, All)
{
    // Exclude primary match from matches list.
    clear(me.xa);
    _fillXa(me, prefix(matches, primaryPos));
    _fillXa(me, suffix(matches, primaryPos + 1));
    appendAlignments(me.record, me.xa);
}

template <typename TSpec, typename Traits, typename TMatches, typename TCount, typename TPos>
inline void _fillXaImpl(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TCount bestCount, TPos primaryPos, Strata)
{
    if (primaryPos < bestCount)
    {
        clear(me.xa);
        TMatches const & cooptimal = prefix(matches, bestCount);
        _fillXa(me, prefix(cooptimal, primaryPos));
        _fillXa(me, suffix(cooptimal, primaryPos + 1));
        appendAlignments(me.record, me.xa);
    }
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _fillXa(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    typedef typename Iterator<TMatches const, Standard>::Type   TIter;

    TIter itEnd = end(matches, Standard());
    for (TIter it = begin(matches, Standard()); it != itEnd; ++it)
    {
        append(me.xa, nameStore(me.outputCtx)[getContigId(value(it))]);
        appendValue(me.xa, ',');
        // TODO(esiragusa): convert contig begin and end to string.
//        append(me.xa, getContigBegin(value(it)) + 1);
//        appendValue(me.xa, ',');
//        append(me.xa, getContigEnd(value(it)) + 1);
//        appendValue(me.xa, ',');
        appendValue(me.xa, onForwardStrand(value(it)) ? '+' : '-');
        appendValue(me.xa, ',');
        // TODO(esiragusa): convert errors to string.
        appendValue(me.xa, '0' + getErrors(value(it)));
        appendValue(me.xa, ';');
    }
}

// ----------------------------------------------------------------------------
// Function _writeRecord()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits>
inline void _writeRecord(MatchesWriter<TSpec, Traits> & me)
{
    _writeRecordImpl(me, typename Traits::TThreading());
}

template <typename TSpec, typename Traits, typename TThreading>
inline void _writeRecordImpl(MatchesWriter<TSpec, Traits> & me, TThreading const & /* tag */)
{
    write2(me.outputStream, me.record, me.outputCtx, typename Traits::TOutputFormat());
}

template <typename TSpec, typename Traits>
inline void _writeRecordImpl(MatchesWriter<TSpec, Traits> & me, Parallel)
{
    write2(me.recordBuffer, me.record, me.outputCtx, typename Traits::TOutputFormat());

    if (length(me.recordBuffer) > Power<2, 16>::VALUE)
        _writeRecordBufferImpl(me, Parallel());
}

// ----------------------------------------------------------------------------
// Function _writeRecordBufferImpl()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TThreading>
inline void _writeRecordBufferImpl(MatchesWriter<TSpec, Traits> & /* me */, TThreading const & /* tag */) {}

template <typename TSpec, typename Traits>
inline void _writeRecordBufferImpl(MatchesWriter<TSpec, Traits> & me, Parallel)
{
    SEQAN_OMP_PRAGMA(critical(MatchesWriter_writeRecord))
    streamWriteBlock(me.outputStream, begin(me.recordBuffer, Standard()), length(me.recordBuffer));
    clear(me.recordBuffer);
}

#endif  // #ifndef APP_YARA_MAPPER_WRITER_H_
