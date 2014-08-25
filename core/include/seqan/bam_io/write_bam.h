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
 *
 * @param[in,out] stream  The @link StreamConcept Stream @endlink to write to.
 * @param[out]    record  The @link BamAlignmentRecord @endlink object to write out.
 * @param[out]    header  The @link BamHeader @endlink object to write out.
 * @param[in,out] context The @link BamIOContext @endlink object to use.
 * @param[in]     tag     The format tag, one of <tt>Sam</tt> and <tt>Bam</tt>.
 *
 * @return int A status code, 0 on success, != 0 on failure.
 */


template <typename TTarget, typename TNameStore, typename TNameStoreCache>
void write(TTarget & target,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache> & context,
           Bam const & /*tag*/)
{
    write(target, "BAM\1");
    clear(context.buffer);

    // Create text of header.
    for (unsigned i = 0; i < length(header.records); ++i)
        write(context.buffer, header.records[i], context, Sam());
    
    // Note that we do not write out a null-character to terminate the header.  This would be valid by the SAM standard
    // but the samtools do not expect this and write out the '\0' when converting from BAM to SAM.
    // appendValue(context.buffer, '\0');

    // Write text header.
    appendRawPod(target, (__int32)length(context.buffer));
    write(target, context.buffer);

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
           BamIOContext<TNameStore, TNameStoreCache> & context,
           Bam const & /*tag*/)
{
    typedef typename Iterator<CharString, Standard>::Type                                   TCharIter;
    typedef typename Iterator<String<CigarElement<> > const, Standard>::Type __restrict__   TCigarIter;
    typedef typename Iterator<IupacString const, Standard>::Type __restrict__               TSeqIter;
    typedef typename Iterator<CharString const, Standard>::Type __restrict__                TQualIter;

    // First, write record to buffer.

    // set internal lengths.
    record._l_qname = length(record.qName) + 1;
    record._n_cigar = length(record.cigar);
    record._l_qseq = length(record.seq);

    resize(context.buffer, sizeof(BamAlignmentRecordCore) + record._l_qname +
           record._n_cigar * 4 + (record._l_qseq + 1) / 2 + record._l_qseq);
    TCharIter it = begin(context.buffer, Standard());

    // bin_mq_nl
    SEQAN_ASSERT_LT(length(record.qName) + 1u, 255u);
    unsigned l = 0;
    _getLengthInRef(record.cigar, l);
    record.bin =_reg2Bin(record.beginPos, record.beginPos + l);

    // BamAlignmentRecordCore.
    arrayCopyForward(reinterpret_cast<char*>(&record),
                     reinterpret_cast<char*>(&record) + sizeof(BamAlignmentRecordCore),
                     it);
    it += sizeof(BamAlignmentRecordCore);

    // read_name
    arrayCopyForward(begin(record.qName, Standard()),
                     end(record.qName, Standard()),
                     it);
    it += length(record.qName);
    *it++ = 0;

    // cigar
    static unsigned char const MAP[256] =
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
    TCigarIter citEnd = end(record.cigar, Standard());
    for (TCigarIter cit = begin(record.cigar, Standard()); cit != citEnd; ++cit)
        *it++ = ((__uint32)cit->count << 4) | MAP[(unsigned char)cit->operation];

    // seq
    TSeqIter sit = begin(record.seq, Standard());
    TSeqIter sitEnd = sit + (record._l_qseq & ~1);
    while (sit != sitEnd)
    {
        unsigned char x = (ordValue(getValue(sit++)) << 4);
        *it++ = x | ordValue(getValue(sit++));
    }
    if (record._l_qseq & 1)
        *it++ = ordValue(getValue(sit++)) << 4;

    // qual
    SEQAN_ASSERT_LEQ(length(record.qual), length(record.seq));
    TCharIter itQEnd = it + record._l_qseq;
    TQualIter qit = begin(record.qual, Standard());
    TQualIter qitEnd = end(record.qual, Standard());
    while (qit != qitEnd)
        *it++ = *qit++ - '!';
    while (it != itQEnd)
        *it++ = '\xff';         // fill with zero qualities

    // tags
    arrayCopyForward(begin(record.tags, Standard()),
                     end(record.tags, Standard()),
                     it);

    // buffer to stream
    appendRawPod(target, (__uint32)length(context.buffer));
    write(target, begin(context.buffer, Standard()), length(context.buffer));
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
