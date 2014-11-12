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
    typedef typename Traits::TReads            TReads;
    typedef typename Traits::TMatchesSet       TMatchesSet;
    typedef typename Traits::TMatches          TMatches;
    typedef typename Traits::TCigarSet         TCigarSet;
    typedef typename Traits::TOutputFile       TOutputFile;
    typedef typename Traits::TReadsContext     TReadsContext;

    // Thread-private data.
    BamAlignmentRecord      record;
    CharString              recordBuffer;
    CharString              xa;

    // Shared-memory read-write data.
    TOutputFile &           outputFile;

    // Shared-memory read-only data.
    TMatchesSet const &     matchesSet;
    TMatches const &        primaryMatches;
    TCigarSet const &       cigarSet;
    TReadsContext const &   ctx;
    TReads const &          reads;
    Options const &         options;

    MatchesWriter(TOutputFile & outputFile,
                  TMatchesSet const & matchesSet,
                  TMatches const & primaryMatches,
                  TCigarSet const & cigarSet,
                  TReadsContext const & ctx,
                  TReads const & reads,
                  Options const & options) :
        outputFile(outputFile),
        matchesSet(matchesSet),
        primaryMatches(primaryMatches),
        cigarSet(cigarSet),
        ctx(ctx),
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function fillHeader()
// ----------------------------------------------------------------------------

template <typename TOptions>
inline void fillHeader(BamHeader & header, TOptions const & options)
{
    typedef BamHeaderRecord::TTag   TTag;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.4"));
    appendValue(firstRecord.tags, TTag("SO", "unsorted"));
    appendValue(header, firstRecord);

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "Yara"));
    appendValue(pgRecord.tags, TTag("PN", "Yara"));
    appendValue(pgRecord.tags, TTag("VN", options.version));
    appendValue(pgRecord.tags, TTag("CL", options.commandLine));
    appendValue(header, pgRecord);

    // Fill read group header line.
    BamHeaderRecord rgRecord;
    rgRecord.type = BAM_HEADER_READ_GROUP;
    appendValue(rgRecord.tags, TTag("ID", options.readGroup));
    appendValue(rgRecord.tags, TTag("SM", options.readGroup));
//    appendValue(rgRecord.tags, TTag("PI", options.libraryLength));
    appendValue(rgRecord.tags, TTag("PG", "Yara"));
    appendValue(header, rgRecord);

    // NOTE(esiragusa): Reference header line is filled by writeRecord().
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

template <typename TString>
inline void appendReadGroup(BamAlignmentRecord & record, TString const & rg)
{
    appendTagValue(record.tags, "RG", rg, 'Z');
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
    appendReadGroup(me.record, me.options.readGroup);
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
    TSize bestCount = countMatchesInBestStratum(matches);
    _fillReadInfo(me, matches, bestCount);

    if (!me.options.outputSecondary)
        _fillXa(me, matches, 0u);

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
    TSize bestCount = countMatchesInBestStratum(matches);
    _fillReadInfo(me, matches, bestCount);

    // Find the primary match in the list of matches.
    TIter it = findMatch(matches, primary);
    TSize primaryPos = position(it, matches);

    if (!me.options.outputSecondary)
        _fillXa(me, matches, primaryPos);

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
    forEach(matches, [&](typename Value<TMatches const>::Type const & match)
    {
        clear(me.record);
        _fillReadName(me, getReadSeqId(match, me.reads.seqs));
        _fillReadPosition(me, match);
        appendExtraPosition(me.record, getMember(match, ContigEnd()));
        me.record.flag |= BAM_FLAG_SECONDARY;
        _writeRecord(me);
    });
}

// ----------------------------------------------------------------------------
// Function _fillReadName()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TReadSeqId>
inline void _fillReadName(MatchesWriter<TSpec, Traits> & me, TReadSeqId readSeqId)
{
    typedef MatchesWriter<TSpec, Traits>                        TMatchesWriter;
    typedef typename TMatchesWriter::TReads                     TReads;
    typedef typename TReads::TSeqNames const                    TSeqNames;
    typedef typename Value<TSeqNames>::Type                     TSeqName;
    typedef typename DirectionIterator<TSeqName, Input>::Type   TSeqNameIt;

    TSeqName const & seqName = me.reads.names[getReadId(me.reads.seqs, readSeqId)];

    TSeqNameIt seqNameIt = begin(seqName);
    readUntil(me.record.qName, seqNameIt, OrFunctor<IsSpace, IsSlash>());

//    me.record.qName = prefix(seqName, lastOf(seqName, IsSpace()));
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

    me.record.rID = getMember(match, ContigId());
    me.record.beginPos = getMember(match, ContigBegin());
    appendErrors(me.record, getMember(match, Errors()));
}

// ----------------------------------------------------------------------------
// Function _fillReadAlignment()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatch>
inline void _fillReadAlignment(MatchesWriter<TSpec, Traits> & me, TMatch const & match)
{
    me.record.cigar = me.cigarSet[getMember(match, ReadId())];
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
    me.record.flag |= BAM_FLAG_MULTIPLE;

    if (!isMapped(me.ctx, getMateId(me.reads.seqs, readId)))
        me.record.flag |= BAM_FLAG_NEXT_UNMAPPED;

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
    me.record.flag |= BAM_FLAG_ALL_PROPER;

    if (onReverseStrand(mate))
        me.record.flag |= BAM_FLAG_NEXT_RC;

    me.record.rNextId = getMember(mate, ContigId());
    me.record.pNext = getMember(mate, ContigBegin());

    if (getMember(match, ContigId()) == getMember(mate, ContigId()))
    {
        if (getMember(match, ContigBegin()) < getMember(mate, ContigBegin()))
            me.record.tLen = getMember(mate, ContigEnd()) - getMember(match, ContigBegin());
        else
            me.record.tLen = getMember(mate, ContigBegin()) - getMember(match, ContigEnd());
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
    _fillMapq(me, bestCount);
    appendCooptimalCount(me.record, bestCount);
    appendSuboptimalCount(me.record, length(matches) - bestCount);
    appendType(me.record, bestCount == 1);
    appendReadGroup(me.record, me.options.readGroup);
    // Set number of secondary alignments and hit index.
//    appendTagValue(me.record.tags, "NH", 1, 'i');
//    appendTagValue(me.record.tags, "HI", 1, 'i');
}

// ----------------------------------------------------------------------------
// Function _fillXa()
// ----------------------------------------------------------------------------

template <typename TSpec, typename Traits, typename TMatches, typename TPos>
inline void _fillXa(MatchesWriter<TSpec, Traits> & me, TMatches const & matches, TPos primaryPos)
{
    clear(me.xa);

    // Exclude primary match from matches list.
    _fillXa(me, prefix(matches, primaryPos));
    _fillXa(me, suffix(matches, std::min(primaryPos + 1, (TPos)length(matches))));

    appendAlignments(me.record, me.xa);
}

template <typename TSpec, typename Traits, typename TMatches>
inline void _fillXa(MatchesWriter<TSpec, Traits> & me, TMatches const & matches)
{
    forEach(matches, [&](typename Value<TMatches const>::Type const & match)
    {
        append(me.xa, nameStore(context(me.outputFile))[getMember(match, ContigId())]);
        appendValue(me.xa, ',');
        appendNumber(me.xa, getMember(match, ContigBegin()) + 1);
        appendValue(me.xa, ',');
        appendNumber(me.xa, getMember(match, ContigEnd()) + 1);
        appendValue(me.xa, ',');
        appendValue(me.xa, onForwardStrand(match) ? '+' : '-');
        appendValue(me.xa, ',');
        appendNumber(me.xa, getMember(match, Errors()));
        appendValue(me.xa, ';');
    });
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
    writeRecord(me.outputFile, me.record);
}

template <typename TSpec, typename Traits>
inline void _writeRecordImpl(MatchesWriter<TSpec, Traits> & me, Parallel)
{
    write(me.recordBuffer, me.record, context(me.outputFile), me.outputFile.format);

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
    write(me.outputFile.stream, me.recordBuffer);

    clear(me.recordBuffer);
}

#endif  // #ifndef APP_YARA_MAPPER_WRITER_H_
