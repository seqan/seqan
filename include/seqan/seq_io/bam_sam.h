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
// Author: Sebastian Proft <sebastian.proft@fu-berlin.de>
// Author: Anton Komissarov <anton.komissarov@fu-berlin.de>
// Author: Temesgen Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================
// Input on BAM and SAM files.
// ==========================================================================

#include <seqan/bam_io.h>

#ifndef SEQAN_SEQ_IO_BAM_SAM_H_
#define SEQAN_SEQ_IO_BAM_SAM_H_

namespace seqan
{

// ----------------------------------------------------------------------------
// Class SamIgnoreFunctor_
// ----------------------------------------------------------------------------
template <typename TAlphabet>
struct SamIgnoreFunctor_
{
    // ignore only newline if the target alphabet is a char
    // ignore whitespace as well for all other alphabets
    typedef typename If<Or<IsSameType<TAlphabet, char>,
                           Or<IsSameType<TAlphabet, signed char>,
                                IsSameType<TAlphabet, unsigned char> > >,
                                    IsNewline, IsWhitespace >::Type Type;
};

template <typename TAlphabet>
struct SamIgnoreOrAssertFunctor_
{
    typedef typename SamIgnoreFunctor_<TAlphabet>::Type               TIgnore;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Sam>   TAsserter;

    // don't assert in case of char alphabets
    // assert being part of the alphabet for other alphabets
    typedef typename If<Or<IsSameType<TAlphabet, char>,
                            Or<IsSameType<TAlphabet, signed char>,
                                IsSameType<TAlphabet, unsigned char> > >,
                                    TIgnore, OrFunctor<TIgnore, TAsserter> >::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _skipHeader()                                             SamHeader
// ----------------------------------------------------------------------------
template <typename TFile>
inline void _skipHeader(TFile & file, Sam const & /* format */)
{
    while (nextIs(file.iter, SamHeader()))
    {
        skipLine(file.iter);
    }
}
// ----------------------------------------------------------------------------
// Function _skipHeader()                                             BamHeader
// ----------------------------------------------------------------------------
template <typename TFile>
inline void _skipHeader(TFile & file, Bam const & /* format */)
{
    BamHeader header;
    BamFileIn bamIn;
    open(bamIn, file.stream);
    file.context.bamIOContext = context(bamIn);
    readHeader(header, file.context.bamIOContext, file.iter, Bam());
}

// ----------------------------------------------------------------------------
// Function readRecord                                                     Bam
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TCharIter>
inline int32_t readBamRecord(TIdString & meta,
                             TSeqString & seq,
                             TIdString & prevQName,
                             int32_t & remainingBytes,
                             TCharIter & it)
{
    typedef typename Iterator<TSeqString, Standard>::Type             TSeqIter;

    // BamAlignmentRecordCore.
    BamAlignmentRecordCore recordCore;
    arrayCopyForward(it, it + sizeof(BamAlignmentRecordCore), reinterpret_cast<char *>(&recordCore));
    enforceLittleEndian(recordCore);
    it += sizeof(BamAlignmentRecordCore);

    clear(meta);
    clear(seq);

    remainingBytes -= sizeof(BamAlignmentRecordCore) + recordCore._l_qname +
                      recordCore._n_cigar * 4 + (recordCore._l_qseq + 1) / 2 + recordCore._l_qseq;
    SEQAN_ASSERT_GEQ(remainingBytes, 0);

    // query name.
    resize(meta, recordCore._l_qname - 1, Exact());
    arrayCopyForward(it, it + recordCore._l_qname - 1, begin(meta, Standard()));
    it += recordCore._l_qname;

    if (prevQName == meta)
        return 0;

    // skip cigar string.
    it += sizeof(uint32_t) * recordCore._n_cigar;

    // query sequence.
    resize(seq, recordCore._l_qseq, Exact());
    TSeqIter sit = begin(seq, Standard());
    TSeqIter sitEnd = sit + (recordCore._l_qseq & ~1);
    while (sit != sitEnd)
    {
        unsigned char ui = getValue(it);
        ++it;
        *sit =  Iupac(ui >> 4);
        ++sit;
        *sit =  Iupac(ui & 0x0f);
        ++sit;
    }
    if (recordCore._l_qseq & 1)
        *sit++ = Iupac((uint8_t)*it++ >> 4);
    return recordCore._l_qseq; // returns length of query sequence. Needed for reading quality
}


// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TForwardIter>
inline bool readSamRecord(TIdString & meta,
                          TSeqString & seq,
                          TIdString & prevQName,
                          TForwardIter & iter)
{
    typedef typename Value<TSeqString>::Type                            TSeqAlphabet;
    typedef typename SamIgnoreOrAssertFunctor_<TSeqAlphabet>::Type      TSeqIgnoreOrAssert;

    // fail, if we read "@" (did you miss to call readRecord(header, bamFile) first?)
    if (nextIs(iter, SamHeader()))
        SEQAN_THROW(ParseError("Unexpected SAM header encountered."));

    OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Sam> > nextEntry;

    clear(meta);
    clear(seq);

    // QNAME
    readUntil(meta, iter, nextEntry);
    skipOne(iter, IsTab());

    // if we have seen this sequence before, return
    if (meta == prevQName)
        return false;

    // Skip FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN
    for (unsigned i = 0; i < 8; ++i)
    {
        skipUntil(iter, IsTab());
        skipOne(iter, IsTab());
    }

    // SEQ
    readUntil(seq, iter, nextEntry, TSeqIgnoreOrAssert());
    // Handle case of missing sequence:  Clear seq string as documented.
    if (seq == "*")
        clear(seq);
    skipOne(iter, IsTab());
    return true;
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFile>
inline void _readRecord(TIdString & meta,
                        TSeqString & seq,
                        TQualString & qual,
                        TFile & file,
                        Sam const & /* format */)
{
    typedef typename Value<TQualString>::Type TQualAlphabet;
    typedef typename SamIgnoreOrAssertFunctor_<TQualAlphabet>::Type TQualIgnoreOrAssert;

    if (SEQAN_UNLIKELY(!(context(file).headerWasRead)))
    {
        _skipHeader(file, Sam());
        file.context.headerWasRead = true;
    }
    TIdString pId = toCString(file.context.prevId);
    clear(qual);

    if (readSamRecord(meta, seq, pId, file.iter))
    {
        // QUAL
        readUntil(qual, file.iter, OrFunctor<IsTab, IsNewline>(), TQualIgnoreOrAssert());
        clear(file.context.prevId);
        file.context.prevId = meta;
    }
    else
    {
        clear(file.context.prevId);
        file.context.prevId = meta;
        clear(meta);
    }
    skipLine(file.iter);
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFile>
inline void _readRecord(TIdString & meta,
                        TSeqString & seq,
                        TQualString & qual,
                        TFile & file,
                        Bam const & /* format */)
{
    typedef typename Iterator<CharString, Standard>::Type                             TCharIter;
    typedef typename Iterator<TQualString, Standard>::Type SEQAN_RESTRICT             TQualIter;

    if (SEQAN_UNLIKELY(!(context(file).headerWasRead)))
    {
        _skipHeader(file, Bam());
        file.context.headerWasRead = true;
    }
    TIdString pId = toCString(file.context.prevId);

    clear(qual);

    // Read size and data of the remaining block in one chunk (fastest).
    int32_t remainingBytes = _readBamRecordWithoutSize(file.context.bamIOContext.buffer, file.iter);
    TCharIter it = begin(file.context.bamIOContext.buffer, Standard());


    if (int32_t l_qseq = readBamRecord(meta, seq, pId, remainingBytes, it))
    {
        // phred quality
        resize(qual, l_qseq, Exact());
        TQualIter qitEnd = end(qual, Standard());
        for (TQualIter qit = begin(qual, Standard()); qit != qitEnd; )
        {
            *qit++ = '!' + *it++;
        }
        clear(file.context.prevId);
        file.context.prevId = meta;
    }
    else
    {
        clear(file.context.prevId);
        file.context.prevId = meta;
        clear(meta);
    }

    // skip tags
    it += remainingBytes;
}

// ----------------------------------------------------------------------------
// Function readRecord(SAM); With qualities
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TQualString, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >, void)
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           FormattedFile<Fastq, Input, TSpec> & file,
           Sam const & /* format */)
{
    clear(qual);
    _readRecord(meta, seq, qual, file, Sam());
}
// ----------------------------------------------------------------------------
// Function readRecord(SAM); Without qualities
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >, void)
readRecord(TIdString & meta,
           TSeqString & seq,
           FormattedFile<Fastq, Input, TSpec> & file,
           Sam const & /* format */)
{
    readRecord(meta, seq, context(file).buffer[2], file, Sam());
}

// ----------------------------------------------------------------------------
// Function readRecord(BAM); With qualities
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TQualString, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >, void)
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           FormattedFile<Fastq, Input, TSpec> & file,
           Bam const & /* format */)
{
    clear(qual);
    _readRecord(meta, seq, qual, file, Bam());
}
// ----------------------------------------------------------------------------
// Function readRecord(BAM); Without qualities
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TSpec>
inline SEQAN_FUNC_ENABLE_IF(Is<InputStreamConcept<typename FormattedFile<Fastq, Input, TSpec>::TStream> >, void)
readRecord(TIdString & meta,
           TSeqString & seq,
           FormattedFile<Fastq, Input, TSpec> & file,
           Bam const & /* format */)
{
    readRecord(meta, seq, context(file).buffer[2], file, Bam());
}

inline void fillHeader(BamHeader & header)
{
    typedef BamHeaderRecord::TTag   TTag;

    BamHeaderRecord seqRecord;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "2.3"));
    appendValue(firstRecord.tags, TTag("SO", "sorted"));
    appendValue(header, firstRecord);

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "SeqAn_IO"));
    appendValue(pgRecord.tags, TTag("PN", "SeqAn_IO"));
    appendValue(header, pgRecord);

    // Fill read group header line.
    BamHeaderRecord rgRecord;
    rgRecord.type = BAM_HEADER_READ_GROUP;
    appendValue(rgRecord.tags, TTag("ID", "none"));
    appendValue(rgRecord.tags, TTag("SM", "none"));
    appendValue(rgRecord.tags, TTag("PG", "SeqAn_IO"));
    appendValue(header, rgRecord);

}
// ----------------------------------------------------------------------------
// Function _writeRecord(SAM);
// ----------------------------------------------------------------------------
template <typename TFile, typename TIdString, typename TSeqString, typename TQualString>
inline void _writeRecord(TFile & file,
                         TIdString const & meta,
                         TSeqString const & seq,
                         TQualString const & qual,
                         Sam const & /* format */)
{
    if (SEQAN_UNLIKELY(!(context(file).headerWasWriten)))
    {
        BamHeader header;
        fillHeader(header);

        StringSet<CharString> contigNameStore;
        appendValue(contigNameStore, "*");
        NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
        BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

        file.context.bamIOContext = bamIOContext;
        write(file.iter, header, context(file).bamIOContext, Sam());
        file.context.headerWasWriten = true;
    }

    BamAlignmentRecord rec;
    clear(rec);
    rec.qName = meta;
    rec.seq = seq;
    rec.qual = qual;
    rec.flag = BAM_FLAG_UNMAPPED;
    
    write(file.iter, rec, context(file).bamIOContext, Sam());
}

// ----------------------------------------------------------------------------
// Function _writeRecord(BAM);
// ----------------------------------------------------------------------------
template <typename TFile, typename TIdString, typename TSeqString, typename TQualString>
inline void _writeRecord(TFile & file,
                         TIdString const & meta,
                         TSeqString const & seq,
                         TQualString const & qual,
                         Bam const & /* format */)
{
    if (SEQAN_UNLIKELY(!(context(file).headerWasWriten)))
    {
        BamHeader header;
        fillHeader(header);

        StringSet<CharString> contigNameStore;
        appendValue(contigNameStore, "*");
        NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
        BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

        file.context.bamIOContext = bamIOContext;
        write(file.iter, header, context(file).bamIOContext, Bam());
        file.context.headerWasWriten = true;
    }

    BamAlignmentRecord rec;
    clear(rec);
    rec.qName = meta;
    rec.seq = seq;
    rec.qual = qual;
    rec.flag = BAM_FLAG_UNMAPPED;
    
    write(file.iter, rec, context(file).bamIOContext, Bam());
}

// ----------------------------------------------------------------------------
// Function writeRecord(SAM); Separate Qualities
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline void writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
                        TIdString const & meta,
                        TSeqString const & seq,
                        TQualString const & qual,
                        Sam const & /* format */)
{
    _writeRecord(file, meta, seq, qual, Sam());
}

// ----------------------------------------------------------------------------
// Function writeRecord(SAM); Qualities inside seq
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(HasQualities<typename Value<TSeqString>::Type> , void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
            TIdString const & meta,
            TSeqString const & seq,
            Sam const & /* format */)
{
    typedef QualityExtractor<typename Value<TSeqString>::Type> TQualityExtractor;
    ModifiedString<TSeqString const, ModView<TQualityExtractor> > quals(seq);
    writeRecord(file, meta, seq, quals, Sam());
}

// ----------------------------------------------------------------------------
// Function writeRecord(SAM); Without qualities inside seq
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Not<HasQualities<typename Value<TSeqString>::Type> >, void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
                        TIdString const & meta,
                        TSeqString const & seq,
                        Sam const & /* format */)
{
        writeRecord(file, meta, seq, "*", Sam());
}
// ----------------------------------------------------------------------------
// Function writeRecord(BAM); Separate Qualities
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline void writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
                        TIdString const & meta,
                        TSeqString const & seq,
                        TQualString const & qual,
                        Bam const & /* format */)
{
    _writeRecord(file, meta, seq, qual, Bam());
}

// ----------------------------------------------------------------------------
// Function writeRecord(BAM); Qualities inside seq
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(HasQualities<typename Value<TSeqString>::Type> , void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
            TIdString const & meta,
            TSeqString const & seq,
            Bam const & /* format */)
{
    typedef QualityExtractor<typename Value<TSeqString>::Type> TQualityExtractor;
    ModifiedString<TSeqString const, ModView<TQualityExtractor> > quals(seq);
    CharString qual = quals;
    writeRecord(file, meta, seq, quals, Bam());
}

// ----------------------------------------------------------------------------
// Function writeRecord(BAM); Without qualities inside seq
// ----------------------------------------------------------------------------
template <typename TSpec, typename TIdString, typename TSeqString>
inline SEQAN_FUNC_ENABLE_IF(Not<HasQualities<typename Value<TSeqString>::Type> >, void)
writeRecord(FormattedFile<Fastq, Output, TSpec> & file,
            TIdString const & meta,
            TSeqString const & seq,
            Bam const & /* format */)
{
    writeRecord(file, meta, seq, "*", Bam());
}
// ----------------------------------------------------------------------------
// Function readRecord(Bam)  -> Directly from an iterator -> not Implemented
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TIterator>
inline void readRecord(TIdString & /* meta */, TSeqString & /* seq */, TIterator & /* iter */, Bam const & /* format */)
{}
// ----------------------------------------------------------------------------
// Function readRecord(Bam)  -> Directly from an iterator -> not Implemented
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TQualString, typename TIterator>
inline void readRecord(TIdString & /* meta */,
                       TSeqString & /* seq */,
                       TQualString & /* qual */,
                       TIterator & /* iter */,
                       Bam const & /* format */)
{}
// ----------------------------------------------------------------------------
// Function readRecord(Sam)  -> Directly from an iterator -> not Implemented
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TIterator>
inline void readRecord(TIdString & /* meta */, TSeqString & /* seq */, TIterator & /* iter */, Sam const & /* format */)
{}
// ----------------------------------------------------------------------------
// Function readRecord(Sam)  -> Directly from an iterator -> not Implemented
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TQualString, typename TIterator>
inline void
readRecord(TIdString & /* meta */,
           TSeqString & /* seq */,
           TQualString & /* qual */,
           TIterator & /* iter */,
           Sam const & /* format */)
{}

} // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_BAM_SAM_H_
