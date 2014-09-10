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
// Code for reading Bam.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @defgroup SamBamIO SAM/BAM I/O
 * @brief Tags for identifying SAM/BAM format.
 */

/*!
 * @tag SamBamIO#Bam
 * @brief Identify the BAM format.
 *
 * @tag SamBamIO#Sam
 * @brief Identify the SAM format.
 */

/**
.Tag.Bam
..cat:BAM I/O
..signature:Bam
..summary:Tag for identifying the BAM format.
..include:seqan/bam_io.h
..see:Tag.Sam
*/

struct Bam_;
typedef Tag<Bam_> Bam;


template <typename T>
struct FileFormatExtensions<Bam, T>
{
    static char const * VALUE[1];	// default is one extension
};

template <typename T>
char const * FileFormatExtensions<Bam, T>::VALUE[1] =
{
    ".bam"     // default output extension
};


template <typename T>
struct MagicHeader<Bam, T>
{
    static unsigned char const VALUE[4];
};

template <typename T>
unsigned char const MagicHeader<Bam, T>::VALUE[4] = { 'B', 'A', 'M', '\1' };  // BAM's magic header

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()                                              BamHeader
// ----------------------------------------------------------------------------

/*!
 * @fn SamBamIO#readRecord
 * @brief Read a record from a SAM/BAM file.
 *
 * @signature int readRecord(record, context, stream, tag);
 *
 * @param[out]    record  The @link BamAlignmentRecord @endlink object to read the information into.
 * @param[out]    header  The @link BamHeader @endlink object to read the header information into.
 * @param[in,out] context The @link BamIOContext @endlink object to use.
 * @param[in,out] stream  The @link StreamConcept Stream @endlink to read from.
 * @param[in]     tag     The format tag, one of <tt>Sam</tt> and <tt>Bam</tt>.
 *
 * @return int A status code, 0 on success, != 0 on failure.
 */

/**
.Function.readRecord
..signature:readRecord(headerRecord, context, stream, tag)
..param.header:@Class.BamHeader@ to read information into.
...type:Class.BamHeader
..param.context:The context to use for reading.
...type:Class.BamIOContext
..param.stream:The stream to read from (for BAM).
...remarks:BAM data can be read from any stream. For the proper decompression (from compressed BAM, the default) use @Spec.BGZF Stream@.
...type:Concept.StreamConcept
..param.tag:Format to read @Class.BamHeader@ from.
...type:Tag.Sam
...type:Tag.Bam
..include:seqan/bam_io.h
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
void readRecord(BamHeader & header,
                BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
                TForwardIter & iter,
                Bam const & /*tag*/)
{
    // Read BAM magic string.
    String<char, Array<4> > magic;
    readUntil(magic, iter, CountDownFunctor<>(4));
    if (magic != "BAM\1")
        throw std::runtime_error("Not in BAM format.");

    // Read header text, including null padding.
    __int32 lText;
    readRawByte(lText, iter);

    CharString samHeader;
    readUntil(samHeader, iter, CountDownFunctor<>(lText));

    // Truncate to first position of '\0'.
    Iterator<CharString, Rooted>::Type it = begin(samHeader);
    for (; !atEnd(it); ++it)
        if (*it == '\0')
            break;
    resize(samHeader, (it - begin(samHeader))); 

    // Parse out header records.
    BamHeaderRecord headerRecord;
    it = begin(samHeader);
    while (!atEnd(it))
    {
        clear(headerRecord);
        readRecord(headerRecord, context, it, Sam());
        appendValue(header.records, headerRecord);
    }

    // Read # reference sequences.
    __int32 nRef;
    readRawByte(nRef, iter);
    CharString name;

    clear(context.translateFile2GlobalRefId);
    resize(context.translateFile2GlobalRefId, nRef, -1);

    for (__int32 i = 0; i < nRef; ++i)
    {
        // Read length of the reference name.
        __int32 nName;
        readRawByte(nName, iter);
        clear(name);
        readUntil(name, iter, CountDownFunctor<>(nName));
        resize(name, nName - 1);
        // Read length of the reference sequence.
        __int32 lRef;
        readRawByte(lRef, iter);

        // Store sequence info.
        typedef typename BamHeader::TSequenceInfo TSequenceInfo;
        appendValue(header.sequenceInfos, TSequenceInfo(name, lRef));
        // Append contig name to name store, if not known already.
        context.translateFile2GlobalRefId[i] = getIdByName(nameStoreCache(context), name);
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                     BamAlignmentRecord
// ----------------------------------------------------------------------------

/**
.Function.readRecord
..signature:readRecord(alignmentRecord, context, stream, tag)
..param.alignmentRecord.type:Class.BamAlignmentRecord
*/

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BamAlignmentRecord & record,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Bam const & /*tag*/)
{
    typedef typename Iterator<CharString, Standard>::Type                           TCharIter;
    typedef typename Iterator<String<CigarElement<> >, Standard>::Type __restrict__ TCigarIter;
    typedef typename Iterator<IupacString, Standard>::Type __restrict__             TSeqIter;
    typedef typename Iterator<CharString, Standard>::Type __restrict__              TQualIter;

    // Read size of the remaining block.
    __int32 remainingBytes = 0;
    readRawByte(remainingBytes, iter);

    // Read remaining block in one chunk (fastest).
    clear(context.buffer);
    write(context.buffer, iter, size_t(remainingBytes));
    TCharIter it = begin(context.buffer, Standard());

    // BamAlignmentRecordCore.
    arrayCopyForward(it, it + sizeof(BamAlignmentRecordCore), reinterpret_cast<char*>(&record));
    it += sizeof(BamAlignmentRecordCore);

    remainingBytes -= sizeof(BamAlignmentRecordCore) + record._l_qname +
                      record._n_cigar * 4 + (record._l_qseq + 1) / 2 + record._l_qseq;
    SEQAN_ASSERT_GEQ(remainingBytes, 0);

    // Translate file local rID into a global rID that is compatible with the context nameStore.
    if (record.rID >= 0 && !empty(context.translateFile2GlobalRefId))
        record.rID = context.translateFile2GlobalRefId[record.rID];
    if (record.rID >= 0)
        SEQAN_ASSERT_LT(static_cast<__uint64>(record.rID), length(nameStore(context)));

    // query name.
    resize(record.qName, record._l_qname - 1, Exact());
    arrayCopyForward(it, it + record._l_qname - 1, begin(record.qName, Standard()));
    it += record._l_qname;

    // cigar string.
    resize(record.cigar, record._n_cigar, Exact());
    static char const * CIGAR_MAPPING = "MIDNSHP=X*******";
    TCigarIter cigEnd = end(record.cigar, Standard());
    for (TCigarIter cig = begin(record.cigar, Standard()); cig != cigEnd; ++cig)
    {
        unsigned opAndCnt = *reinterpret_cast<unsigned * &>(it)++;
        cig->operation = CIGAR_MAPPING[opAndCnt >> 16];
        cig->count = opAndCnt & 0xffff;
    }

    // query sequence.
    resize(record.seq, record._l_qseq, Exact());
    TSeqIter sit = begin(record.seq, Standard());
    TSeqIter sitEnd = sit + (record._l_qseq & ~1);
    while (sit != sitEnd)
    {
        unsigned char ui = getValue(it);
        ++it;
        assignValue(sit, Iupac(ui >> 4));
        ++sit;
        assignValue(sit, Iupac(ui & 0x0f));
        ++sit;
    }
    if (record._l_qseq & 1)
        *sit++ = Iupac((__uint8)*it++ >> 4);

    // phred quality
    resize(record.qual, record._l_qseq, Exact());
    // If qual is a sequence of 0xff (heuristic same as samtools: Only look at first byte) then we clear it, to get the
    // representation of '*';
    TQualIter qitEnd = end(record.qual, Standard());
    for (TQualIter qit = begin(record.qual, Standard()); qit != qitEnd;)
        *qit++ = '!' + *it++;
    if (!empty(record.qual) && record.qual[0] == '\xff')
        clear(record.qual);

    // tags
    resize(record.tags, remainingBytes, Exact());
    arrayCopyForward(it, it + remainingBytes, begin(record.tags, Standard()));
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_READ_BAM_H_
