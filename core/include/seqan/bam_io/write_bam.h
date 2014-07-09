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
// Code for writing BAM.
// ==========================================================================

// TODO(holtgrew): Add buffer to context?
// TODO(holtgrew): Rename to writeRecord from write2! Go over deprecated alias!

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_

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
// Function writeRecord()                                             BamHeader
// ----------------------------------------------------------------------------

/*!
 * @fn SamBamIO#write2
 * @brief Write a record to a SAM/BAM file.
 *
 * @signature int writeRecord(stream, record, context, tag);
 * @signature int writeRecord(stream, header, context, tag);
 *
 * @param[in,out] stream  The @link StreamConcept Stream @endlink to write to.
 * @param[out]    record  The @link BamAlignmentRecord @endlink object to write out.
 * @param[out]    header  The @link BamHeader @endlink object to write out.
 * @param[in,out] context The BamIOContext object to use.
 * @param[in]     tag     The format tag, one of <tt>Sam</tt> and <tt>Bam</tt>.
 *
 * @return int A status code, 0 on success, != 0 on failure.
 */


template <typename TTarget, typename TNameStore, typename TNameStoreCache>
void write(TTarget & target,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache> const & context,
           Bam const & /*tag*/)
{
    write(target, "BAM\1");

    // Create text of header.
    CharString headerBuffer;
    for (unsigned i = 0; i < length(header.records); ++i)
        write(headerBuffer, header.records[i], context, Sam());
    
    // Note that we do not write out a null-character to terminate the header.  This would be valid by the SAM standard
    // but the samtools do not expect this and write out the '\0' when converting from BAM to SAM.
    // appendValue(headerBuffer, '\0');

    // Write text header.
    appendRawPod(target, (__int32)length(headerBuffer));
    write(target, headerBuffer);

    // Write references.
    appendRawPod(target, (__int32)_max(length(header.sequenceInfos), length(nameStore(context))));

    for (unsigned i = 0; i < length(header.sequenceInfos); ++i)
    {
        appendRawPod(target, (__int32)(length(header.sequenceInfos[i].i1) + 1));
        write(target, header.sequenceInfos[i].i1);
        writeValue(target, '\0');
        appendRawPod(target, (__int32)header.sequenceInfos[i].i2);
    }
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                    BamAlignmentRecord
// ----------------------------------------------------------------------------

static inline int _reg2Bin(uint32_t beg, uint32_t end)
{
    --end;
    if (beg >> 14 == end >> 14)
        return 4681 + (beg >> 14);

    if (beg >> 17 == end >> 17)
        return 585 + (beg >> 17);

    if (beg >> 20 == end >> 20)
        return 73 + (beg >> 20);

    if (beg >> 23 == end >> 23)
        return 9 + (beg >> 23);

    if (beg >> 26 == end >> 26)
        return 1 + (beg >> 26);

    return 0;
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache>
void write(TTarget & target,
           BamAlignmentRecord & record,
           BamIOContext<TNameStore, TNameStoreCache> const & /*context*/,
           Bam const & /*tag*/)
{
    // TODO(singer): Once the MAP tables are gone delete the buffer.
    CharString buffer;

    // First, write record to buffer.

    // bin_mq_nl
    SEQAN_ASSERT_LT(length(record.qName) + 1u, 255u);
    record._l_qname = length(record.qName) + 1;
    unsigned l = 0;
    _getLengthInRef(record.cigar, l);
    record.bin =_reg2Bin(record.beginPos, record.beginPos + l);

    // flag_nc
    record._n_cigar = length(record.cigar);

    // l_seq
    record._l_qseq = (__int32)length(record.seq);
    appendRawPod(buffer, static_cast<BamAlignmentRecordCore &>(record));

    // read_name
    write(buffer, record.qName);
    writeValue(buffer, '\0');

    // cigar
    static __uint8 const MAP[256] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0,
        0, 0, 0, 0, 2, 0, 0, 0, 5, 1, 0, 0, 0, 0, 3, 0,
        6, 0, 0, 4, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    for (unsigned i = 0; i < length(record.cigar); ++i)
    {
        __uint32 x = record.cigar[i].count;
        x <<= 4;
        x |= MAP[static_cast<int>(record.cigar[i].operation)];
        appendRawPod(buffer, x);
    }

    // seq
    static __uint8 const MAP2[256] =
    {
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 15, 15,
        15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15,
        15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
    };
    for (size_t i = 0; i < length(record.seq); i+=2)
        writeValue(buffer, (MAP2[ordValue(record.seq[i])] << 4) | MAP2[ordValue(record.seq[i + 1])]);
    if ((length(record.seq) & 1) == 1)
        writeValue(buffer, MAP2[ordValue(back(record.seq))] << 4);

    // qual
    if (empty(record.qual))
    {
        for (unsigned i = 0; i < length(record.qual); ++i)
            writeValue(buffer, (unsigned char)(0xff));
    }
    else
    {
        for (unsigned i = 0; i < length(record.qual); ++i)
            writeValue(buffer, (char)(record.qual[i] - '!'));
    }

    // tags
    if (length(record.tags) > 0u)
        write(buffer, record.tags);

    // buffer to stream
    appendRawPod(target, (__uint32)length(buffer));
    write(target, buffer);
}
}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
