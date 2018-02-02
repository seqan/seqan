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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_BAM_IO_TEST_WRITE_BAM_H_
#define TESTS_BAM_IO_TEST_WRITE_BAM_H_

#include <iomanip>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_io_bam_write_header)
{
    using namespace seqan;

    typedef typename BamHeaderRecord::TTag    TTag;

    // Prepare input.

    StringSet<CharString> contigNameStore;
    appendValue(contigNameStore, "REFERENCE");
    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

    BamHeader header;
    appendValue(contigLengths(bamIOContext), 10000);

    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.0"));
    appendValue(header, firstRecord);

    BamHeaderRecord seqRecord;
    seqRecord.type = BAM_HEADER_REFERENCE;
    appendValue(seqRecord.tags, TTag("SN", "REFERENCE"));
    appendValue(seqRecord.tags, TTag("LN", "10000"));
    appendValue(header, seqRecord);

    // Call code under test.
    String<char> text;
    write(text, header, bamIOContext, Bam());

    // Compare results.
    CharString bamFilename;
    append(bamFilename, getAbsolutePath("/tests/bam_io/header_uncompressed.bam"));

    String<char, MMap<> > EXPECTED;
    open(EXPECTED, toCString(bamFilename));

    /*
    String<char> EXPECTED = "\x42\x41\x4d\x01\x25\x00\x00\x00\x40\x48\x44\x09\x56\x4e\x3a\x31\x2e\x30\x0a\x40\x53\x51\x09\x53\x4e\x3a\x52\x45\x46\x45\x52\x45\x4e\x43\x45\x09\x4c\x4e\x3a\x31\x30\x30\x30\x30\x0a\x01\x00\x00\x00\x0a\x00\x00\x00\x52\x45\x46\x45\x52\x45\x4e\x43\x45\x00\x10\x27\x00\x00";
    */

    SEQAN_ASSERT_EQ(text, EXPECTED);
}

SEQAN_DEFINE_TEST(test_bam_io_bam_write_alignment)
{
    using namespace seqan;

    typedef typename BamHeaderRecord::TTag    TTag;

    // Create input.

    StringSet<CharString> contigNameStore;
    appendValue(contigNameStore, "REFERENCE");
    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

    BamHeader header;
    appendValue(contigLengths(bamIOContext), 10000);

    BamHeaderRecord firstRecord;
    firstRecord.type = BAM_HEADER_FIRST;
    appendValue(firstRecord.tags, TTag("VN", "1.0"));
    appendValue(header, firstRecord);

    BamHeaderRecord seqRecord;
    seqRecord.type = BAM_HEADER_REFERENCE;
    appendValue(seqRecord.tags, TTag("SN", "REFERENCE"));
    appendValue(seqRecord.tags, TTag("LN", "10000"));
    appendValue(header, seqRecord);

    // Call code under test.
    String<char> text;
    write(text, header, bamIOContext, Bam());

    BamAlignmentRecord record;
    record.qName = "READNAME";
    record.flag = BAM_FLAG_ALL_PROPER | BAM_FLAG_RC;
    record.rID = 0;
    record.beginPos = 30;
    record.mapQ = 8;
    appendValue(record.cigar, CigarElement<>('M', 10));
    record.rNextId = BamAlignmentRecord::INVALID_REFID;
    record.pNext = BamAlignmentRecord::INVALID_POS;
    record.tLen = BamAlignmentRecord::INVALID_LEN;
    record.seq  = "CGATCGATAA";
    record.qual = "IIIIIIIIII";

    // Call code under test.
    write(text, record, bamIOContext, Bam());

    CharString bamFilename;
    append(bamFilename, getAbsolutePath("/tests/bam_io/alignment_uncompressed.bam"));

    String<char, MMap<> > EXPECTED;
    open(EXPECTED, toCString(bamFilename));

    /*
    // Compare results.
    String<char> EXPECTED =
        "\x3c\x00\x00\x00\x00\x00\x00\x00\x1e\x00\x00\x00\x09\x08\x49\x12"
        "\x01\x00\x12\x00\x0a\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff"
        "\x00\x00\x00\x00\x52\x45\x41\x44\x4e\x41\x4d\x45\x00\xa0\x00\x00"
        "\x00\x24\x18\x24\x18\x11\x28\x28\x28\x28\x28\x28\x28\x28\x28\x28";
     */
    SEQAN_ASSERT_EQ(text, EXPECTED);
}

#endif  // TESTS_BAM_IO_TEST_WRITE_BAM_H_
