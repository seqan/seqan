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

#ifndef TESTS_BAM_IO_TEST_WRITE_SAM_H_
#define TESTS_BAM_IO_TEST_WRITE_SAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/bam_io.h>

SEQAN_DEFINE_TEST(test_bam_io_sam_write_header)
{
    using namespace seqan;

    typedef typename BamHeaderRecord::TTag    TTag;

    // Prepare input.

    StringSet<CharString> contigNameStore;
    appendValue(contigNameStore, "REF");
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
    appendValue(seqRecord.tags, TTag("SN", "REF"));
    appendValue(seqRecord.tags, TTag("LN", "10000"));
    appendValue(header, seqRecord);

    // Call code under test.
    String<char> buffer;
    write(buffer, header, bamIOContext, Sam());

    // Compare results.
    char const * EXPECTED =
            "@HD\tVN:1.0\n"
            "@SQ\tSN:REF\tLN:10000\n";

    SEQAN_ASSERT_EQ(CharString(EXPECTED), buffer);
}

SEQAN_DEFINE_TEST(test_bam_io_sam_write_alignment)
{
    using namespace seqan;

    // Create input.

    StringSet<CharString> contigNameStore;
    appendValue(contigNameStore, "REF");
    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);

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
    record.seq = "CGATCGATAA";
    record.qual = "IIIIIIIIII";

    // Call code under test.
    String<char> buffer;
    write(buffer, record, bamIOContext, Sam());

    // Compare results.
    char const * EXPECTED = "READNAME\t18\tREF\t31\t8\t10M\t*\t0\t0\tCGATCGATAA\tIIIIIIIIII\n";

    SEQAN_ASSERT_EQ(CharString(EXPECTED), buffer);
}

#endif  // TESTS_BAM_IO_TEST_WRITE_SAM_H_
