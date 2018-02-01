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
// Author: Temesgen Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_SEQ_IO_TEST_WRITE_BAM_H_
#define TESTS_SEQ_IO_TEST_WRITE_BAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

// ---------------------------------------------------------------------------
// Write Sequences with(out) qualities to BamFileIn
// ---------------------------------------------------------------------------

void testSeqIOBamFileWriteSequences(char const * extension, bool withQuals)
{
    std::string tmpPath = (std::string)SEQAN_TEMP_FILENAME() + extension;

    seqan::SeqFileOut seqFileOut(toCString(tmpPath));

    seqan::String<seqan::CharString> metas;
    seqan::String<seqan::Dna5String> seqs;
    seqan::String<seqan::CharString> quals;

    appendValue(metas, "READ0");
    appendValue(metas, "READ1");
    appendValue(metas, "READ2");

    appendValue(seqs, "AAAAAAAAAA");
    appendValue(seqs, "AAAAAAAAAA");
    appendValue(seqs, "AAAAAAAAAA");
    if (withQuals)
    {
        appendValue(quals, "!!!!!!!!!!");
        appendValue(quals, "!!!!!!!!!!");
        appendValue(quals, "!!!!!!!!!!");
    }

    if (withQuals)
    {
        writeRecords(seqFileOut, metas, seqs, quals);
    }
    else
    {
        writeRecords(seqFileOut, metas, seqs);
    }
    close(seqFileOut);

    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(tmpPath)))
    {
        std::cerr << "Could not open " << toCString(tmpPath) << "!\n";
    }

    BamHeader header;
    readHeader(header, bamFileIn);

    BamAlignmentRecord record;
    seqan::String<seqan::CharString> new_metas;
    seqan::String<seqan::Dna5String> new_seqs;
    seqan::String<seqan::CharString> new_quals;

    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);
        appendValue(new_metas, record.qName);
        appendValue(new_seqs, record.seq);
        appendValue(new_quals, record.qual);
    }

    SEQAN_ASSERT_EQ(length(new_metas), 3u);
    SEQAN_ASSERT_EQ(length(new_metas), 3u);
    SEQAN_ASSERT_EQ(length(new_metas), 3u);

    SEQAN_ASSERT_EQ(metas[0], new_metas[0]);
    SEQAN_ASSERT_EQ(metas[1], new_metas[1]);
    SEQAN_ASSERT_EQ(metas[2], new_metas[2]);

    SEQAN_ASSERT_EQ(seqs[0], new_seqs[0]);
    SEQAN_ASSERT_EQ(seqs[1], new_seqs[1]);
    SEQAN_ASSERT_EQ(seqs[2], new_seqs[2]);

    if (withQuals)
    {
        SEQAN_ASSERT_EQ(quals[0], new_quals[0]);
        SEQAN_ASSERT_EQ(quals[1], new_quals[1]);
        SEQAN_ASSERT_EQ(quals[2], new_quals[2]);
    }
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_sam_write_sequences)
{
    testSeqIOBamFileWriteSequences(".sam", false);
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_bam_write_sequences)
{
    testSeqIOBamFileWriteSequences(".bam", false);
}


SEQAN_DEFINE_TEST(test_seq_io_bam_file_sam_write_sequences_and_qualities)
{
    testSeqIOBamFileWriteSequences(".sam", true);
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_bam_write_sequences_and_qualities)
{
    testSeqIOBamFileWriteSequences(".bam", true);
}

#endif  // TESTS_SEQ_IO_TEST_WRITE_BAM_H_
