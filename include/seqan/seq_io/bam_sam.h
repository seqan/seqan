// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
//    typedef
//    TagList<Bam,
//    TagList<Sam
//    > >
//    BamLikeFormats;
//    typedef TagSelector<BamLikeFormats>  BamLikeFormat;

// ----------------------------------------------------------------------------
// Class SamIgnoreFunctor_
// ----------------------------------------------------------------------------

template <typename TAlphabet>
struct SamIgnoreFunctor_
{
    typedef typename If<Or<IsSameType<TAlphabet, char>,
    Or<IsSameType<TAlphabet, signed char>,
    IsSameType<TAlphabet, unsigned char> > >,
    IsNewline, // ignore only newline if the target alphabet is a char
    IsWhitespace // ignore whitespace as well for all other alphabets
    >::Type Type;
};

template <typename TAlphabet>
struct SamIgnoreOrAssertFunctor_
{
    typedef typename SamIgnoreFunctor_<TAlphabet>::Type               TIgnore;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Sam>   TAsserter;

    typedef typename If<Or<IsSameType<TAlphabet, char>,
    Or<IsSameType<TAlphabet, signed char>,
    IsSameType<TAlphabet, unsigned char> > >,
    TIgnore,   // don't assert in case of char alphabets
    OrFunctor<TIgnore, TAsserter> // assert being part of the alphabet for other alphabets
    >::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function skipHeader()                                              SamHeader
// ----------------------------------------------------------------------------
template <typename TFile>
inline void skipHeader(TFile & file)
{
    if (isEqual(file.format, Sam()))
    {
        while (nextIs(file.iter, SamHeader()))
        {
            skipLine(file.iter);
        }
    }
    else if(isEqual(file.format, Bam()))
    {
        BamHeader header;
        BamFileIn bamIn;
        open(bamIn, file.stream);
        file.context.bamIOContext = context(bamIn);
        readHeader(header, file.context.bamIOContext, file.iter, Bam());
    }
}

// ----------------------------------------------------------------------------
// Function readRecord                                                     Bam
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString,
          typename TCharIter>
inline int32_t
_readRecord(TIdString & meta,
            TSeqString & seq,
            TIdString & prevQName,
            int32_t & remainingBytes,
            TCharIter & it,
            Bam const & /* tag */)
{
    typedef typename Iterator<TSeqString, Standard>::Type             TSeqIter;

    // BamAlignmentRecordCore.
    BamAlignmentRecordCore recordCore;
    arrayCopyForward(it, it + sizeof(BamAlignmentRecordCore), reinterpret_cast<char *>(&recordCore));
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
        assignValue(sit, Iupac(ui >> 4));
        ++sit;
        assignValue(sit, Iupac(ui & 0x0f));
        ++sit;
    }
    if (recordCore._l_qseq & 1)
        *sit++ = Iupac((uint8_t)*it++ >> 4);
    return recordCore._l_qseq; // returns length of query sequence. Needed for reading quality
}


// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString,
          typename TForwardIter>
inline bool
_readRecord(TIdString & meta,
            TSeqString & seq,
            TIdString & prevQName,
            TForwardIter & iter,
            Sam const & /*tag*/)
{
    typedef typename Value<TSeqString>::Type                                TSeqAlphabet;
    typedef typename SamIgnoreOrAssertFunctor_<TSeqAlphabet>::Type          TSeqIgnoreOrAssert;

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

template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline void readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, FormattedFile<Fastq, Input, TSpec>&file, Sam)
{
    if (!(context(file).hasReadHeader))
    {
        skipHeader(file);
        file.context.hasReadHeader = true;
    }
    TIdString pId = toCString(file.context.prevId);
    if (isEqual(file.format, Sam()))
    {
        typedef typename Value<TQualString>::Type TQualAlphabet;
        typedef typename SamIgnoreOrAssertFunctor_<TQualAlphabet>::Type TQualIgnoreOrAssert;
        clear(qual);

        if (_readRecord(meta, seq, pId, file.iter, Sam()))
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
    else if (isEqual(file.format, Bam()))
    {
        typedef typename Iterator<CharString, Standard>::Type                             TCharIter;
        typedef typename Iterator<TQualString, Standard>::Type SEQAN_RESTRICT             TQualIter;

        clear(qual);

        // Read size and data of the remaining block in one chunk (fastest).
        int32_t remainingBytes = _readBamRecordWithoutSize(file.context.bamIOContext.buffer, file.iter);
        TCharIter it = begin(file.context.bamIOContext.buffer, Standard());


        if (int32_t l_qseq = _readRecord(meta, seq, pId, remainingBytes, it, Bam()))
        {
            // phred quality
            resize(qual, l_qseq, Exact());
            TQualIter qitEnd = end(qual, Standard());
            for (TQualIter qit = begin(qual, Standard()); qit != qitEnd; )
                *qit++ = '!' + *it++;

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
}

//// ----------------------------------------------------------------------------
//// Function writeRecord(BAM/SAM); Separate Qualities
//// ----------------------------------------------------------------------------
//
template <typename TSpec, typename TIdString, typename TSeqString, typename TQualString>
inline void
writeRecord(FormattedFile<Fastq, Output, TSpec>&file,
            TIdString const & meta,
            TSeqString const & seq,
            TQualString const & qual,
            Sam)
{
    if (!(context(file).hasWritenHeader))
    {
        typedef BamHeaderRecord::TTag   TTag;
        BamHeader header;

        StringSet<CharString> contigNameStore;
        appendValue(contigNameStore, "*");
        NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
        BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

        file.context.bamIOContext = bamIOContext;
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

        if (isEqual(file.format, Sam()))
        {
            write(file.iter, header, context(file).bamIOContext, Sam());
        }
        else if (isEqual(file.format, Bam()))
        {
            write(file.iter, header, context(file).bamIOContext, Bam());
        }
        file.context.hasWritenHeader = true;
    }

    BamAlignmentRecord rec;
    clear(rec);
    rec.qName = meta;
    rec.seq = seq;
    rec.qual = qual;
    rec.flag = BAM_FLAG_UNMAPPED;

    if (isEqual(file.format, Sam()))
    {
        write(file.iter, rec, context(file).bamIOContext, Sam());
    }
    else if (isEqual(file.format, Bam()))
    {
        write(file.iter, rec, context(file).bamIOContext, Bam());
    }
}

} // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_BAM_SAM_H_