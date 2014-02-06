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
// This file contains a reimplementation of Store IO using the Bam IO module.
// ==========================================================================

#ifndef SEQAN_EXTRAS_APPS_MASAI_STORE_IO_SAM_H_
#define SEQAN_EXTRAS_APPS_MASAI_STORE_IO_SAM_H_

namespace seqan {

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function buildMDString()
// --------------------------------------------------------------------------

template <typename TMDString, typename TContigGaps, typename TReadGaps>
inline void
buildMDString(TMDString & md, TContigGaps & contigGaps, TReadGaps & readGaps)
{
    typedef typename Value<TMDString>::Type TMDChar;

    char op;
    char lastOp = ' ';
    unsigned numOps = 0;

    typename Iterator<TContigGaps>::Type contigIt = begin(contigGaps);
    typename Iterator<TReadGaps>::Type readIt = begin(readGaps);

    for (; !atEnd(contigIt) && !atEnd(readIt); goNext(contigIt), goNext(readIt))
    {
        if (isGap(contigIt))
            continue;
        if (isGap(readIt))
            op = 'D';
        else
            op = (convert<char>(*contigIt) == convert<char>(*readIt)) ? 'M' : 'R';

        // append match run
        if (lastOp != op)
        {
            if (lastOp == 'M')
                appendValue(md, numOps);
            numOps = 0;
            lastOp = op;
        }

        // append deleted/replaced reference character
        if (op != 'M')
        {
            // add ^ from non-deletion to deletion
            if (op == 'D' && lastOp != 'D')
                appendValue(md, '^');
            // add 0 from deletion to replacement
            if (op == 'R' && lastOp == 'D')
                appendValue(md, '0');

            appendValue(md, convert<TMDChar>(*contigIt));
        }

        ++numOps;
    }

    if (lastOp == 'M')
        appendValue(md, numOps);
}

// --------------------------------------------------------------------------
// Function buildCigarString()
// --------------------------------------------------------------------------

template <typename TCigar, typename TContigGaps, typename TReadGaps, typename TThresh>
inline void
buildCigarString(TCigar & cigar, TContigGaps & contigGaps, TReadGaps & readGaps, TThresh splicedGapThresh)
{
    char op;
    char lastOp = ' ';
    unsigned numOps = 0;

    typename Iterator<TContigGaps>::Type contigIt = begin(contigGaps);
    typename Iterator<TReadGaps>::Type readIt = begin(readGaps);

    for (; !atEnd(contigIt) && !atEnd(readIt); goNext(contigIt), goNext(readIt))
    {
        if (isGap(contigIt))
        {
            if (isGap(readIt))
                op = 'P';
            else if (isClipped(readIt))
                op = '?';
            else
                op = 'I';
        }
        else if (isClipped(contigIt))
        {
            op = '?';
        }
        else
        {
            if (isGap(readIt))
                op = 'D';
            else if (isClipped(readIt))
                op = 'S';
            else
                op = 'M';
        }

        if (lastOp != op)
        {
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';

            if (numOps > 0)
                appendValue(cigar, CigarElement<>(lastOp, numOps));

            numOps = 0;
            lastOp = op;
        }
        ++numOps;
    }

    if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
        lastOp = 'N';

    if (numOps > 0)
        appendValue(cigar, CigarElement<>(lastOp, numOps));
}

template <typename TCigar, typename TContigGaps, typename TReadGaps>
inline void
buildCigarString(TCigar & cigar, TContigGaps & contigGaps, TReadGaps & readGaps)
{
    buildCigarString(cigar, contigGaps, readGaps, 20);
}

// --------------------------------------------------------------------------
// Function _fillHeader()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void
_fillHeader(FragmentStore<TSpec, TConfig> & store,
            BamHeader & header)
{
    typedef FragmentStore<TSpec, TConfig>                       TFragmentStore;
    typedef typename TFragmentStore::TLibraryStore              TLibraryStore;
    typedef typename TFragmentStore::TContigStore               TContigStore;
    typedef typename TFragmentStore::TNameStore                 TNameStore;

    typedef typename Value<TContigStore>::Type                  TContig;
    typedef typename Iterator<TLibraryStore, Standard>::Type    TLibraryIter;
    typedef typename Iterator<TContigStore, Standard>::Type     TContigIter;
    typedef typename Iterator<TNameStore, Standard>::Type       TContigNameIter;
    typedef typename Id<TContig>::Type                          TId;

    typedef BamHeader::TSequenceInfo                            TSequenceInfo;
    typedef BamHeaderRecord::TTag                               TTag;

    // Fill first header line.
    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.4"));
    appendValue(firstRecord.tags, TTag("SO", "unsorted"));
    appendValue(header.records, firstRecord);

    // Fill sequence info header line.
    TContigIter it          = begin(store.contigStore, Standard());
    TContigIter itEnd       = end(store.contigStore, Standard());
    TContigNameIter nit     = begin(store.contigNameStore, Standard());
    TContigNameIter nitEnd  = end(store.contigNameStore, Standard());

    for (; it != itEnd && nit != nitEnd; ++it, ++nit)
        appendValue(header.sequenceInfos, TSequenceInfo(*nit, length((*it).seq)));

    // Fill program header line.
    BamHeaderRecord pgRecord;
    pgRecord.type = BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, TTag("ID", "SeqAn"));
    appendValue(header.records, pgRecord);

    // Fill library info header line.
    BamHeaderRecord rgRecord;
    rgRecord.type = BAM_HEADER_READ_GROUP;

    TLibraryIter lit    = begin(store.libraryStore, Standard());
    TLibraryIter litEnd = end(store.libraryStore, Standard());

    for (TId id = 0; lit != litEnd; ++lit, ++id)
    {
        appendValue(pgRecord.tags, TTag("ID", id + 1));
        appendValue(pgRecord.tags, TTag("LB", store.libraryNameStore[id]));
        appendValue(pgRecord.tags, TTag("PI", (int)store.libraryStore[id].mean));
        // Sample name needs to be included into fragment store.
        appendValue(pgRecord.tags, TTag("SM", "none"));
    }

}

// --------------------------------------------------------------------------
// Function _fillRecord()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TAlign>
inline void
_fillRecord(FragmentStore<TSpec, TConfig> & store,
            BamAlignmentRecord & record,
            typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type const & alignedRead,
            typename Value<typename FragmentStore<TSpec, TConfig>::TAlignQualityStore>::Type const & alignQuality,
            typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadTagStore>::Type const &,
            typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type const & alignedMate,
            TAlign const & align,
            bool secondaryAlignment,
            bool fillCigar = true)
{
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;

    typedef typename TFragmentStore::TContigStore                   TContigStore;
    typedef typename Value<TContigStore>::Type                      TContig;
    typedef typename TFragmentStore::TContigSeq                     TContigSeq;
    typedef typename TContig::TGapAnchors                           TContigGapAnchors;

    typedef typename TFragmentStore::TReadStore                     TReadStore;
    typedef typename Value<TReadStore>::Type                        TReadStoreElement;

    typedef typename TFragmentStore::TReadSeq                       TReadSeq;
    typedef typename Iterator<TReadSeq, Standard>::Type             TReadSeqIterator;

    typedef Iterator<CharString>::Type                              TCharStringIterator;

//    typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
//    typedef typename Value<TAlignedReadStore>::Type                 TAlignedRead;
//    typedef typename TAlignedRead::TGapAnchors                      TAlignedReadGapAnchors;

    typedef Gaps<TContigSeq, AnchorGaps<TContigGapAnchors> >        TContigGaps;
//    typedef Gaps<TReadSeq, AnchorGaps<TAlignedReadGapAnchors> >     TReadGaps;

//    TReadGaps readGaps(store.readSeqStore[alignedRead.readId], alignedRead.gaps);

    TContigGaps contigGaps(store.contigStore[alignedRead.contigId].seq, store.contigStore[alignedRead.contigId].gaps);
    setBeginPosition(contigGaps, std::min(alignedRead.beginPos, alignedRead.endPos));
    setEndPosition(contigGaps, std::max(alignedRead.beginPos, alignedRead.endPos));

    // Fill QNAME with read name.
//    record.qName = store.readNameStore[alignedRead.readId];
    for (TCharStringIterator it = begin(store.readNameStore[alignedRead.readId]); it != end(store.readNameStore[alignedRead.readId]); ++it)
    {
        if (*it == ' ' || *it == '\t' || *it == '\n' || *it == '\r')
            break;
        appendValue(record.qName, *it);
    }

    // Fill FLAG.
    record.flag = 0;

    if (alignedRead.beginPos > alignedRead.endPos)
        record.flag |= BAM_FLAG_RC;
    if (secondaryAlignment)
        record.flag |= BAM_FLAG_SECONDARY;

    signed char mateNo = getMateNo(store, alignedRead.readId);
    if (mateNo == 0)
        record.flag |= BAM_FLAG_FIRST;
    if (mateNo == 1)
        record.flag |= BAM_FLAG_LAST;

    // This read has a mate.
    if (alignedRead.pairMatchId != TReadStoreElement::INVALID_ID)
    {
        record.flag |= BAM_FLAG_MULTIPLE;

        if (alignedMate.readId != TReadStoreElement::INVALID_ID)
        {
            // Mate is mapped.
            record.flag |= BAM_FLAG_ALL_PROPER;

            if (alignedMate.beginPos > alignedMate.endPos)
                record.flag |= BAM_FLAG_NEXT_RC;
        }
        else
        {
            // Mate is unmapped (actually we should check if the mate has no match at all).
            record.flag |= BAM_FLAG_NEXT_UNMAPPED;
        }
    }

    // Fill RNAME by providing contig id.
    record.rID = alignedRead.contigId;

    // Fill POS with start position.
    record.beginPos = positionGapToSeq(contigGaps, std::min(alignedRead.beginPos, alignedRead.endPos));

    // Fill MAPQ with mapping quality.
//    if (alignQuality.score != TAlignedReadQualityElement::INVALID_SCORE)
    if (alignQuality.score > 0)
        record.mapQ = alignQuality.score;

    // Fill CIGAR using aligned read.
    if (fillCigar)
        buildCigarString(record.cigar, row(align, 0), row(align, 1));

    // Fill RNEXT by providing mate contig id.
    if (alignedMate.contigId != TContig::INVALID_ID)
        record.rNextId = alignedMate.contigId;
    else
        record.rNextId = BamAlignmentRecord::INVALID_REFID;

    // Fill PNEXT with mate start position.
    if (alignedRead.contigId == alignedMate.contigId)
        record.pNext = positionGapToSeq(contigGaps, std::min(alignedMate.beginPos, alignedMate.endPos));
    else
        record.pNext = BamAlignmentRecord::INVALID_POS;

    // Fill TLEN.
    if (alignedRead.contigId == alignedMate.contigId)
    {
        if (alignedRead.beginPos < alignedMate.beginPos)
            record.tLen = positionGapToSeq(contigGaps, std::max(alignedMate.beginPos, alignedMate.endPos) - 1) - record.beginPos + 1;
        else
            record.tLen = record.pNext - positionGapToSeq(contigGaps, std::max(alignedRead.beginPos, alignedRead.endPos) - 1) - 1;
    }
    else
        record.tLen = BamAlignmentRecord::INVALID_LEN;

    if (!secondaryAlignment)
    {
        // TODO(esiragusa):Pass read as argument, already revComplemented
        TReadSeq read = store.readSeqStore[alignedRead.readId];

        if (alignedRead.beginPos > alignedRead.endPos)
            reverseComplement(read);

        // Fill SEQ with read sequence.
        record.seq = read;

        // Fill QUAL with read quality values.
        TReadSeqIterator it    = begin(read, Standard());
        TReadSeqIterator itEnd = end(read, Standard());

        for (; it != itEnd; ++it)
            appendValue(record.qual, (char)(getQualityValue(*it) + 33));
    }

    BamTagsDict tagsDict(record.tags);

    // Fill NM with errors from quality store.
    if (alignQuality.errors != MaxValue<unsigned char>::VALUE)
        setTagValue(tagsDict, "NM", (int)alignQuality.errors);

    // TODO(esiragusa):Fill MD using aligned read.
//    CharString md;
//    buildMDString(md, row(align, 0), row(align, 1));
//    setTagValue(tagsDict, "MD", md, 'Z');

    // TODO(esiragusa):Fill other tags.
//    append(record.tags, alignedTags);
}

// --------------------------------------------------------------------------
// Sam write functions.
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Function _writeHeader()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStream, typename TContext>
inline void _writeHeader(FragmentStore<TSpec, TConfig> & store,
                         TStream & target,
                         TContext & context,
                         Sam)
{
    BamHeader header;

    // Fill header with information from fragment store.
    _fillHeader(store, header);

    // Write header to target.
    write2(target, header, context, Sam());
}

// --------------------------------------------------------------------------
// Function _writeAlignedRead()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStream, typename TContext, typename TAlign>
inline void
_writeAlignedRead(FragmentStore<TSpec, TConfig> & store,
                  TStream & target,
                  TContext & context,
                  typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type & alignedRead,
                  typename Value<typename FragmentStore<TSpec, TConfig>::TAlignQualityStore>::Type & alignQuality,
                  typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadTagStore>::Type const & alignedTags,
                  typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type const & alignedMate,
                  TAlign const & align,
                  bool secondaryAlignment,
                  Sam)
{
    BamAlignmentRecord record;

    // Fill record.
    _fillRecord(store, record, alignedRead, alignQuality, alignedTags, alignedMate, align, secondaryAlignment);

    // Write record to target.
    write2(target, record, context, Sam());
}

// --------------------------------------------------------------------------
// Function _writeUnalignedRead()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStream, typename TContext>
inline void
_writeUnalignedRead(FragmentStore<TSpec, TConfig> & store,
                    TStream & target,
                    TContext & context,
                    typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type & alignedRead,
                    typename Value<typename FragmentStore<TSpec, TConfig>::TAlignQualityStore>::Type & alignQuality,
                    typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadTagStore>::Type const & alignedTags,
                    typename Value<typename FragmentStore<TSpec, TConfig>::TAlignedReadStore>::Type const & alignedMate,
                    bool secondaryAlignment,
                    Sam)
{
    typedef FragmentStore<TSpec, TConfig>                   TFragmentStore;

    typedef typename TFragmentStore::TReadSeqStore          TReadSeqStore;
    typedef typename Value<TReadSeqStore>::Type const       TReadSeq;

    typedef Align<TReadSeq, ArrayGaps>                      TAlign;

    TAlign align;

    // Align read.
    _alignRead(store, align, alignedRead, alignQuality);

    // Write aligned read.
    _writeAlignedRead(store, target, context, alignedRead, alignQuality, alignedTags, alignedMate, align, secondaryAlignment, Sam());
}

// --------------------------------------------------------------------------
// Function _writeAlignedReads()
// --------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStream, typename TContext>
inline void _writeAlignedReads(FragmentStore<TSpec, TConfig> & store, TStream & target, TContext & context, Sam)
{
    typedef FragmentStore<TSpec, TConfig>                           TFragmentStore;

    typedef typename TFragmentStore::TReadStore                     TReadStore;
    typedef typename Value<TReadStore>::Type                        TReadStoreElement;

    typedef typename TFragmentStore::TAlignedReadStore              TAlignedReadStore;
    typedef typename Iterator<TAlignedReadStore, Standard>::Type    TAlignedReadStoreIterator;
    typedef typename Value<TAlignedReadStore>::Type                 TAlignedReadStoreElement;
    typedef typename Id<TAlignedReadStoreElement>::Type             TAlignedReadStoreElementId;

    typedef typename TFragmentStore::TAlignQualityStore             TAlignQualityStore;
    typedef typename Value<TAlignQualityStore>::Type                TAlignQualityStoreElement;

    typedef typename TFragmentStore::TAlignedReadTagStore           TAlignedReadTagStore;
    typedef typename Value<TAlignedReadTagStore>::Type              TAlignedReadTagStoreElement;

    typedef unsigned long                                           TWord;

    // Store outer library size for each pair match (indexed by pairMatchId)
    String<TAlignedReadStoreElementId> mateIndices;
    calculateMateIndices(mateIndices, store);

    // Bitset to signal wether a read was aligned at least once.
    String<TWord>   readAligned;
    const unsigned  wordLen = BitsPerValue<TWord>::VALUE;
    resize(readAligned, (length(store.readStore) + wordLen - 1) / wordLen, (TWord)0);

    // Dummy store elements to cope with missing information.
    TAlignedReadStoreElement        noAlignedMate;
    TAlignQualityStoreElement       noAlignQuality;
    TAlignedReadTagStoreElement     noAlignedTags;

    TAlignedReadStoreIterator it    = begin(store.alignedReadStore, Standard());
    TAlignedReadStoreIterator itEnd = end(store.alignedReadStore, Standard());

    for (; it != itEnd; ++it)
    {
        TAlignedReadStoreElement & alignedRead    = *it;
        TAlignQualityStoreElement & alignQuality  = noAlignQuality;
        TAlignedReadTagStoreElement & alignedTags = noAlignedTags;
        TAlignedReadStoreElement & alignedMate    = noAlignedMate;

        // Try to get quality.
        if (alignedRead.id < length(store.alignQualityStore))
            alignQuality = store.alignQualityStore[alignedRead.id];

        // Try to get tags.
        if (alignedRead.id < length(store.alignedReadTagStore))
            alignedTags = store.alignedReadTagStore[alignedRead.id];

        // Try to get a mate.
        if (alignedRead.pairMatchId != TReadStoreElement::INVALID_ID)
        {
            TAlignedReadStoreElementId mateIndex = mateIndices[2 * alignedRead.pairMatchId + getMateNo(store, alignedRead.readId)];
            alignedMate = store.alignedReadStore[mateIndex];
        }

        // Test for secondary alignment.
        TWord mask     = (TWord)1 << (alignedRead.readId % wordLen);
        bool secondary = (readAligned[alignedRead.readId / wordLen] & mask) != 0;
        readAligned[alignedRead.readId / wordLen] |= mask;

        // Write unaligned read.
        _writeUnalignedRead(store, target, context, alignedRead, alignQuality, alignedTags, alignedMate, secondary, Sam());
    }
}

// --------------------------------------------------------------------------
// Function write()
// --------------------------------------------------------------------------

///.Function.write.param.tag.type:Tag.File Format.tag.Sam

template <typename TSpec, typename TConfig, typename TStream>
inline void write(FragmentStore<TSpec, TConfig> & store, TStream & target, Sam)
{
    typedef FragmentStore<TSpec, TConfig>               TFragmentStore;
    typedef typename TFragmentStore::TContigNameStore   TContigNameStore;

    typedef BamIOContext<TContigNameStore>              TBamIOContext;

    TBamIOContext context(store.contigNameStore, store.contigNameStoreCache);

    _writeHeader(store, target, context, Sam());
    _writeAlignedReads(store, target, context, Sam());
}

}

#endif //#ifndef SEQAN_EXTRAS_APPS_MASAI_STORE_IO_SAM_H_
