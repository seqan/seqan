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
// ==========================================================================
// Input on BAM and SAM files.
// ==========================================================================

#include <seqan/bam_io.h>

#ifndef SEQAN_SEQ_IO_BAM_SAM_H_
#define SEQAN_SEQ_IO_BAM_SAM_H_

namespace seqan
{
// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector);                    SAM/BAM without quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TForwardIter,
        typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(TIdString & /*meta*/,
           TSeqString & /*seq*/,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /*context*/,
           TForwardIter & /*iter*/,
           TagSelector<> const & /*format*/)

{
    SEQAN_FAIL("BamFileIn: File format not specified.");
}

template <typename TIdString, typename TSeqString, typename TForwardIter,
        typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;
    if (isEqual(format, TFormat()))
        readRecord( meta, seq, context, iter, TFormat() );
    else
        readRecord( meta, seq, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(TagSelector);                       SAM/BAM with quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TForwardIter,
        typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(TIdString & /*meta*/,
           TSeqString & /*seq*/,
           TQualString & /*qual*/,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /*context*/,
           TForwardIter & /*iter*/,
           TagSelector<void> const & /*format*/)
{
    SEQAN_FAIL("BamFileIn: File format not specified.");
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TForwardIter,
        typename TNameStore, typename TNameStoreCache, typename TStorageSpec, typename TTagList>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           TagSelector<TTagList> const & format)
{
    typedef typename TTagList::Type TFormat;
    if (isEqual(format, TFormat()))
        readRecord( meta, seq, qual, context, iter, TFormat() );
    else
        readRecord( meta, seq, qual, context, iter, static_cast<typename TagSelector<TTagList>::Base const &>(format));
}

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn);                      SAM/BAM without quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TSpec>
inline void
readRecord(TIdString & meta, TSeqString & seq, FormattedFile<Bam, Input, TSpec> & file)
{
    readRecord(meta, seq, context(file), file.iter, file.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn);                         SAM/BAM with quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString, typename TSpec>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, FormattedFile<Bam, Input, TSpec> & fileIn)
{
    readRecord(meta, seq, qual, context(fileIn), fileIn.iter, fileIn.format);
}

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn);                          Bam without quality
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString,
        typename TCharIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline int32_t
_readRecord(TIdString & meta,
            TSeqString & seq,
            CharString & prevQName,
            BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
            int32_t & remainingBytes,
            TCharIter & it,
            Bam const & /* tag */)
{
    typedef typename Iterator<TSeqString, Standard>::Type SEQAN_RESTRICT             TSeqIter;

    // BamAlignmentRecordCore.
    BamAlignmentRecordCore recordCore;
    arrayCopyForward(it, it + sizeof(BamAlignmentRecordCore), reinterpret_cast<char*>(&recordCore));
    it += sizeof(BamAlignmentRecordCore);

    clear(meta);
    clear(seq);

    remainingBytes -= sizeof(BamAlignmentRecordCore) + recordCore._l_qname +
                      recordCore._n_cigar * 4 + (recordCore._l_qseq + 1) / 2 + recordCore._l_qseq;
    SEQAN_ASSERT_GEQ(remainingBytes, 0);

    // Translate file local rID into a global rID that is compatible with the context contigNames.
    if (recordCore.rID >= 0 && !empty(context.translateFile2GlobalRefId))
        recordCore.rID = context.translateFile2GlobalRefId[recordCore.rID];
    if (recordCore.rID >= 0)
        SEQAN_ASSERT_LT(static_cast<uint64_t>(recordCore.rID), length(contigNames(context)));

    // ... the same for rNextId
    if (recordCore.rNextId >= 0 && !empty(context.translateFile2GlobalRefId))
        recordCore.rNextId = context.translateFile2GlobalRefId[recordCore.rNextId];
    if (recordCore.rNextId >= 0)
        SEQAN_ASSERT_LT(static_cast<uint64_t>(recordCore.rNextId), length(contigNames(context)));

    // query name.
    resize(meta, recordCore._l_qname - 1, Exact());
    arrayCopyForward(it, it + recordCore._l_qname - 1, begin(meta, Standard()));
    it += recordCore._l_qname;

    if (prevQName == meta)
        return 0;

    // skip cigar string.
    it += sizeof(uint32_t)*recordCore._n_cigar;

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

template <typename TIdString, typename TSeqString,
        typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bam const & tag)
{
    typedef typename Iterator<TIdString, Standard>::Type                             TCharIter;

    //save previous sequence name
    CharString prevQName = context.buffer;

    clear(context.buffer);

    // Read size and data of the remaining block in one chunk (fastest).
    int32_t remainingBytes = _readBamRecordWithoutSize(context.buffer, iter);
    TCharIter it = begin(context.buffer, Standard());

    if ( _readRecord( meta, seq, prevQName, context, remainingBytes, it, tag ) )
    {
        clear(context.buffer);
        context.buffer=meta;
    } else
    {
        clear(context.buffer);
        context.buffer=meta;
        clear(meta);
    }

    // skip phred quality and tags
    it += remainingBytes;
}

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn)                              BAM with quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TQualString,
        typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bam const & tag )
{
    typedef typename Iterator<TIdString, Standard>::Type                              TCharIter;
    typedef typename Iterator<TQualString, Standard>::Type SEQAN_RESTRICT             TQualIter;

    CharString prevQName = context.buffer;

    clear(context.buffer);
    clear(qual);

    // Read size and data of the remaining block in one chunk (fastest).
    int32_t remainingBytes = _readBamRecordWithoutSize(context.buffer, iter);
    TCharIter it = begin(context.buffer, Standard());

    if ( int32_t l_qseq = _readRecord( meta, seq, prevQName, context, remainingBytes, it, tag ) )
    {
        // phred quality
        resize(qual, l_qseq, Exact());
        TQualIter qitEnd = end(qual, Standard());
        for (TQualIter qit = begin(qual, Standard()); qit != qitEnd;)
            *qit++ = '!' + *it++;

        // Handle case of missing quality: throw parse exception if there is no quality.
        // there is another version of readRecord that doesn't use quality
        // If qual is a sequence of 0xff (heuristic same as samtools: Only look at first byte) then we stop the program

        if (!empty(qual) && qual[0] == '\xff')
            throw ParseError("This BAM file doesn't provide PHRED quality string. "
                                     "Consider using another version of readRecord without quality");
        clear(context.buffer);
        context.buffer=meta;
    } else
    {
        clear(context.buffer);
        context.buffer=meta;
        clear(meta);
    }

    // skip tags
    it += remainingBytes;
}

// ----------------------------------------------------------------------------
// Class SamIgnoreFunctor_
// ----------------------------------------------------------------------------

template <typename TAlphabet>
struct SamIgnoreFunctor_
{
    typedef typename If< Or< IsSameType<TAlphabet, char>,
    Or< IsSameType<TAlphabet, signed char>,
    IsSameType<TAlphabet, unsigned char> > >,
    IsNewline,                     // ignore only newline if the target alphabet is a char
    IsWhitespace                   // ignore whitespace as well for all other alphabets
    >::Type Type;
};

template <typename TAlphabet>
struct SamIgnoreOrAssertFunctor_
{
    typedef typename SamIgnoreFunctor_<TAlphabet>::Type               TIgnore;
    typedef AssertFunctor<IsInAlphabet<TAlphabet>, ParseError, Sam>   TAsserter;

    typedef typename If< Or< IsSameType<TAlphabet, char>,
    Or< IsSameType<TAlphabet, signed char>,
    IsSameType<TAlphabet, unsigned char> > >,
    TIgnore,                       // don't assert in case of char alphabets
    OrFunctor<TIgnore, TAsserter>  // assert being part of the alphabet for other alphabets
    >::Type Type;
};

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn)                           SAM without quality
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString,
        typename TForwardIter>
inline bool
_readRecord(TIdString & meta,
            TSeqString & seq,
            CharString & prevQName,
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
    if ( meta == prevQName )
        return false;

    // Skip FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN
    for (unsigned i = 0; i < 8; ++i) {
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

template <typename TIdString, typename TSeqString,
        typename TNameStore, typename  TNameStoreCache, typename TStorageSpec,
        typename TForwardIter>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Sam const & tag )
{
    CharString prevQName = context.buffer;
    clear( context.buffer );

    if ( _readRecord(meta, seq, prevQName, iter, tag) )
    {
        context.buffer = meta;
    } else
    {
        context.buffer = meta;
        clear(meta);
    }
    skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function readRecord(BamFileIn)                              SAM with quality
// ----------------------------------------------------------------------------
template <typename TIdString, typename TSeqString, typename TQualString,
        typename TNameStore, typename  TNameStoreCache, typename TStorageSpec,
        typename TForwardIter>
inline void
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Sam const & tag)
{
    typedef typename Value<TQualString>::Type TQualAlphabet;
    typedef typename SamIgnoreOrAssertFunctor_<TQualAlphabet>::Type TQualIgnoreOrAssert;

    clear(qual);
    CharString prevQName = context.buffer;
    clear( context.buffer );

    if ( _readRecord( meta, seq, prevQName, iter, tag) )
    {
        // QUAL
        readUntil(qual, iter, OrFunctor<IsTab, IsNewline>(), TQualIgnoreOrAssert());

        // Handle case of missing quality: throw parse exception if there is no quality.
        // there is another version of readRecord that doesn't use quality
        if (qual == '*') {
            throw ParseError("This SAM file doesn't provide PHRED quality string. "
                                     "Consider using another version of readRecord without quality");
        }
        context.buffer = meta;
    } else
    {
        context.buffer = meta;
        clear(meta);
    }
    skipLine(iter);
}

} // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_BAM_SAM_H_