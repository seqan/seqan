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
// Code for writing SAM.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

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
// Function write2()                                            BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache>
int write(TTarget & target,
           BamHeaderRecord const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & /*context*/,
           Sam const & /*tag*/)
{
    char const * headerTypes[] = {"@HD", "@SQ", "@RG", "@PG", "@CO"};
<<<<<<< HEAD
    streamPut(stream, headerTypes[header.type]);
    if (header.type == BAM_HEADER_COMMENT)
    {
        streamPut(stream, header.tags[0].i2);
=======
    write(target, headerTypes[header.type]);

    if (header.type == BAM_HEADER_COMMENT && !empty(header.tags))
    {
        writeValue(target, '\t');
        write(target, header.tags[0].i2);
>>>>>>> ffe8366... [FIX, FEATURE] Adjusted parts of the sam/bam module to the new seq_io module.
    }
    else
    {
        for (unsigned i = 0; i < length(header.tags); ++i)
        {
            writeValue(target, '\t');
            write(target, header.tags[i].i1);
            writeValue(target, ':');
            write(target, header.tags[i].i2);
        }
    }

    writeValue(target, '\n');
    return 0;
}

// ----------------------------------------------------------------------------
// Function write2()                                                  BamHeader
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache>
int write(TTarget & target,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Sam const & tag)
{
    std::set<CharString> writtenSeqInfos;

    for (unsigned i = 0; i < length(header.records); ++i)
    {
        BamHeaderRecord const & record = header.records[i];
        if (record.type == BAM_HEADER_REFERENCE)
        {
            for (unsigned i = 0; i < length(record.tags); ++i)
            {
                if (record.tags[i].i1 == "SN")
                {
                    writtenSeqInfos.insert(record.tags[i].i2);
                    break;
                }
            }
        }

        int res = write(target, record, context, tag);
        if (res != 0)
            return res;
    }

    // Write missing @SQ header records.
    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
    {
        if (writtenSeqInfos.find(header.sequenceInfos[i].i1) != writtenSeqInfos.end())
            continue;
        write(target, "@SQ\tSN:");
        write(target, header.sequenceInfos[i].i1);
        write(target, "\tLN:");
        appendNumber(target, header.sequenceInfos[i].i2);
        writeValue(target, '\n');
        return 0;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function write2()                                         BamAlignmentRecord
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache>
int write(TTarget & target,
           BamAlignmentRecord const & record,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Sam const & /*tag*/)
{
    write(target, record.qName);
    writeValue(target, '\t');

    appendNumber(target, record.flag);
    writeValue(target, '\t');

    if (record.rID == BamAlignmentRecord::INVALID_REFID)
        writeValue(target, '*');
    else
        write(target, nameStore(context)[record.rID]);

    writeValue(target, '\t');

    if (record.rID == BamAlignmentRecord::INVALID_REFID)
        writeValue(target, '*');
    else
        appendNumber(target, record.beginPos + 1);

    writeValue(target, '\t');

    appendNumber(target, static_cast<__uint16>(record.mapQ));
    writeValue(target, '\t');

    if (empty(record.cigar))
        writeValue(target, '*');
    else
        for (unsigned i = 0; i < length(record.cigar); ++i)
        {
            appendNumber(target, record.cigar[i].count);
            writeValue(target, record.cigar[i].operation);
        }

    writeValue(target, '\t');

    if (record.rNextId == BamAlignmentRecord::INVALID_REFID)
        writeValue(target, '*');
    else if (record.rID == record.rNextId)
        writeValue(target, '=');
    else
        write(target, nameStore(context)[record.rNextId]);

    writeValue(target, '\t');

    if (record.pNext == BamAlignmentRecord::INVALID_POS)
        writeValue(target, '0');
    else
        appendNumber(target, record.pNext + 1);

    writeValue(target, '\t');

    if (record.tLen == BamAlignmentRecord::INVALID_LEN)
        writeValue(target, '0');
    else
        appendNumber(target, record.tLen);

    writeValue(target, '\t');

    if (empty(record.seq))
        writeValue(target, '*');  // Case of empty seq string / "*".
    else
        write(target, record.seq);

    writeValue(target, '\t');


    if (empty(record.qual))  // Case of empty quality string / "*".
        writeValue(target, '*');
    else
        write(target, record.qual);

    if (!empty(record.tags))
    {
        writeValue(target, '\t');
        CharString buffer;
        assignTagsBamToSam(buffer, record.tags);
        write(target, buffer);
    }

    writeValue(target, '\n');

    return 0;

#undef SEQAN_PUT_TAB
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_SAM_H_
