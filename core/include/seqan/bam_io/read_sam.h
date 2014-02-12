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
// Code for reading SAM.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.Sam
..cat:BAM I/O
..signature:Sam
..summary:Tag for identifying the SAM format.
..include:seqan/bam_io.h
..see:Tag.Bam
*/

struct Sam_;
typedef Tag<Sam_> const Sam;

enum SamTokenizeErrors_
{
    SAM_INVALID_RECORD = 2048
};

struct SamHeader_;
typedef Tag<SamHeader_> SamHeader;

struct SamAlignment_;
typedef Tag<SamAlignment_> SamAlignment;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nextIs()                                                  SamHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline bool nextIs(TForwardIter & iter, SamHeader const & /*tag*/)
{
    if (atEnd(iter))
        return false;
    return value(iter) == '@';
}

// ----------------------------------------------------------------------------
// Function skipRecord()                                              SamHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TPass>
inline int skipRecord(TForwardIter & iter,
                      SamHeader const & tag)
{
    if (!nextIs(iter, tag))
        return SAM_INVALID_RECORD;
    skipLine(iter);
    return 0;
}

// ----------------------------------------------------------------------------
// Function skipRecord()                                           SamAlignment
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TPass>
inline int skipRecord(TForwardIter & iter,
                      SamAlignment const & tag)
{
    if (!nextIs(iter, tag))
        return SAM_INVALID_RECORD;
    skipLine(iter);
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                        BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
void readRecord(BamHeaderRecord & record,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               TForwardIter & iter,
               Sam const & /*tag*/)
{
    clear(record);

    // Make sure the first character is '@'.
    skipOne(iter, EqualsChar<'@'>());

    // Read the header tag.
    char c1;
    readOne(c1, iter);
    char c2;
    readOne(c2, iter);

    // Determine header type.
    if (c1 == 'H' && c2 == 'D')
        record.type = BAM_HEADER_FIRST;
    else if (c1 == 'S' && c2 == 'Q')
        record.type = BAM_HEADER_REFERENCE;
    else if (c1 == 'R' && c2 == 'G')
        record.type = BAM_HEADER_READ_GROUP;
    else if (c1 == 'P' && c2 == 'G')
        record.type = BAM_HEADER_PROGRAM;
    else if (c1 == 'C' && c2 == 'O')
        record.type = BAM_HEADER_COMMENT;
    else
    {
        throw std::runtime_error("Invalid Record!");
        return;
    }

    if (record.type == BAM_HEADER_COMMENT)
    {
        skipOne(iter, IsTab());
        CharString &buffer = context.buffer;
        readLine(buffer, iter);
        appendValue(record.tags, Pair<CharString>(CharString(), buffer));
    }
    else
    {
        // Read the rest of the line into the tag field of record.
        CharString key, val;
        while (!atEnd(iter) && value(iter) == '\t')
        {
            clear(key);
            clear(val);

            skipOne(iter, IsTab());
            readUntil(key, iter, EqualsChar<':'>());
            skipOne(iter);
            readUntil(val, iter, OrFunctor<IsTab, IsNewline>());

            appendValue(record.tags, Pair<CharString>(key, val));
        }
    }

<<<<<<< HEAD
    // Skip remaining line break.
    int res = skipLine(reader);
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return res;
    return 0;
=======
    // Skip remaining line break, in case of comment, we already skipped over it.
    if (record.type != BAM_HEADER_COMMENT)
        skipLine(iter);
>>>>>>> ffe8366... [FIX, FEATURE] Adjusted parts of the sam/bam module to the new seq_io module.
}

// ----------------------------------------------------------------------------
// Function readRecord()                                              BamHeader
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(headerRecord, context, recordReader, tag)
..param.recordReader:The RecordReader to read from.
...type:Class.RecordReader
...remarks:Use for SAM.
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
int readRecord(BamHeader & header,
               BamIOContext<TNameStore, TNameStoreCache> & context,
               TForwardIter & iter,
               Sam const & tag)
{
    BamHeaderRecord record;
    while (nextIs(iter, SamHeader()))
    {
        clear(record);
        readRecord(record, context, iter, tag);
        appendValue(header.records, record);

        // Get sequence information from @SQ header.
        if (record.type == BAM_HEADER_REFERENCE)
        {
            CharString sn = "unknown";
            unsigned ln = 0;
            for (unsigned i = 0; i < length(record.tags); ++i)
            {
                if (record.tags[i].i1 == "SN")
                {
                    sn = record.tags[i].i2;
                }
                else if (record.tags[i].i1 == "LN")
                {
                    if (!lexicalCast<unsigned>(ln, record.tags[i].i2))
                        ln = 0;
                }
            }
            typedef typename BamHeader::TSequenceInfo TSequenceInfo;
            appendValue(header.sequenceInfos, TSequenceInfo(sn, ln));

            // Add name to name store cache if necessary.
            unsigned unusedId = 0;
            ignoreUnusedVariableWarning(unusedId);
            if (!getIdByName(nameStore(context), sn, unusedId, nameStoreCache(context)))
                appendName(nameStore(context), sn, nameStoreCache(context));
        }
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                                     BamAlignmentRecord
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(alignmentRecord, context, recordReader, tag)
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache>
inline void readRecord(BamAlignmentRecord & record,
            BamIOContext<TNameStore, TNameStoreCache> & context,
            TForwardIter & iter,
            Sam const & /*tag*/)
{
    clear(record);
    CharString &buffer = context.buffer;

    // QNAME
    readUntil(record.qName, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    skipOne(iter, IsTab());

    // FLAG
    // TODO(holtgrew): Interpret hex and char as c-samtools -X does?
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    record.flag = lexicalCast<__uint16>(buffer);
    skipOne(iter, IsTab());

    // RNAME
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    if (buffer == "*")
    {
        record.rID = BamAlignmentRecord::INVALID_REFID;
    }
    else if (buffer == "0")
    {
        record.rID = BamAlignmentRecord::INVALID_REFID;
    }
    else if (!getIdByName(nameStore(context), buffer, record.rID, nameStoreCache(context)))
    {
        record.rID = length(nameStore(context));
        appendName(nameStore(context), buffer, nameStoreCache(context));
    }
    skipOne(iter, IsTab());

    // POS
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    if (buffer == "*")
        record.beginPos = BamAlignmentRecord::INVALID_POS;
    else if (buffer == "0")
        record.beginPos = BamAlignmentRecord::INVALID_POS;
    else
        record.beginPos = lexicalCast<__uint32>(buffer) - 1;
    skipOne(iter, IsTab());

    // MAPQ
    clear(buffer);
    if (value(iter) == '*')
    {
        record.mapQ = 255;
        skipOne(iter);
    }
    else
    {
        readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
        record.mapQ = lexicalCast<__uint16>(buffer);
    }
    skipOne(iter, IsTab());

    // CIGAR
    CigarElement<> element;
    if (value(iter) == '*')
        skipOne(iter);
    else
    {
        do
        {
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsAlpha, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
            element.count = lexicalCast<__uint32>(buffer);
            element.operation = value(iter);
            skipOne(iter);
            appendValue(record.cigar, element);
        } while (value(iter) != '\t');
    }
    skipOne(iter, IsTab());

    // RNEXT
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    if (buffer == "*")
    {
        record.rNextId = BamAlignmentRecord::INVALID_REFID;
    }
    else if (buffer == "=")
    {
        record.rNextId = record.rID;
    }
    else if (!getIdByName(nameStore(context), buffer, record.rNextId, nameStoreCache(context)))
    {
        record.rNextId = length(nameStore(context));
        appendName(nameStore(context), buffer, nameStoreCache(context));
    }
    skipOne(iter, IsTab());

    // PNEXT
    if (value(iter) == '*')
    {
        record.pNext = BamAlignmentRecord::INVALID_POS;
        skipOne(iter);
    }
    else
    {
        clear(buffer);
        readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
        if (buffer == "0")
            record.pNext = BamAlignmentRecord::INVALID_POS;
        else
            record.pNext = lexicalCast<__uint32>(buffer) - 1;
    }
    skipOne(iter, IsTab());

    // TLEN
    if (value(iter) == '*')
    {
        record.tLen = MaxValue<__int32>::VALUE;
        skipOne(iter);
    }
    else
    {
        clear(buffer);
        if (value(iter) == '-')
            readOne(buffer, iter);

        readUntil(buffer, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
        record.tLen = lexicalCast<__int32>(buffer);
    }
    skipOne(iter, IsTab());

    // SEQ
    readUntil(record.seq, iter, OrFunctor<IsTab, Asserter<NotFunctor<IsNewline>, ParseError, Sam> >());
    // Handle case of missing sequence:  Clear seq string as documented.
    if (record.seq == "*")
        clear(record.seq);
    skipOne(iter, IsTab());

    // QUAL
    IsNewline isNewline;
    readUntil(record.qual, iter, OrFunctor<IsTab, IsNewline>());
   
    // Handle case of missing quality:  Clear qual string as documented.
    if (record.qual == "*")
        clear(record.qual);

    // The following list of tags is optional.  A line break or EOF could also follow.
    if (atEnd(iter))
        return;
    if (value(iter) != '\t')
    {
        skipLine(iter);
        return;
    }
    skipOne(iter, IsTab());

    // TAGS
    clear(buffer);
    readLine(buffer, iter);
    assignTagsSamToBam(record.tags, buffer);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_SAM_H_

