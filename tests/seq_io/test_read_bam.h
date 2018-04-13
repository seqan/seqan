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
// Author: Anton Komissarov <anton.komissarov@fu-berlin.de>
// Author: Sebastian Proft <sebastian.proft@fu-berlin.de>
// Author: Temesgen Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_SEQ_IO_TEST_READ_BAM_H_
#define TESTS_SEQ_IO_TEST_READ_BAM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

// ---------------------------------------------------------------------------
// Read Sequences without qualities from BamFileIn
// ---------------------------------------------------------------------------

void testSeqIOBamFileReadSequences(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::SeqFileIn seqFileIn(toCString(filePath));

    seqan::String<seqan::CharString> metas;
    seqan::String<seqan::Dna5String> seqs;
    while (!atEnd(seqFileIn))
    {
        resize(metas, length(metas) + 1);
        resize(seqs, length(seqs) + 1);
        readRecord(back(metas), back(seqs), seqFileIn);
    }
    SEQAN_ASSERT_EQ(length(metas), 3u);
    SEQAN_ASSERT_EQ(length(seqs), 3u);

    SEQAN_ASSERT_EQ(metas[0], "READ0");
    SEQAN_ASSERT_EQ(seqs[0], "AAAAAAAAAA");

    SEQAN_ASSERT_EQ(metas[1], "READ1");
    SEQAN_ASSERT_EQ(seqs[1], "AAAAAAAAAA");

    SEQAN_ASSERT_EQ(metas[2], "");
    SEQAN_ASSERT_EQ(seqs[2], "");
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_sam_read_sequences)
{
    testSeqIOBamFileReadSequences("/tests/seq_io/small_sequences.sam");
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_bam_read_sequences)
{
    testSeqIOBamFileReadSequences("/tests/seq_io/small_sequences.bam");
}

// ---------------------------------------------------------------------------
// Read Sequences with qualities from BamFileIn
// ---------------------------------------------------------------------------

void testSeqIOBamFileReadSequencesAndQualities(char const * pathFragment)
{
    seqan::CharString filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, pathFragment);

    seqan::SeqFileIn seqFileIn(toCString(filePath));

    seqan::String<seqan::CharString> metas;
    seqan::String<seqan::Dna5String> seqs;
    seqan::String<seqan::CharString> quals;

    while (!atEnd(seqFileIn))
    {
        resize(metas, length(metas) + 1);
        resize(seqs, length(seqs) + 1);
        resize(quals, length(quals) + 1);
        readRecord(back(metas), back(seqs), back(quals), seqFileIn);
    }
    SEQAN_ASSERT_EQ(length(metas), 3u);
    SEQAN_ASSERT_EQ(length(seqs), 3u);
    SEQAN_ASSERT_EQ(length(quals), 3u);

    SEQAN_ASSERT_EQ(metas[0], "READ0");
    SEQAN_ASSERT_EQ(seqs[0], "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(quals[0], "!!!!!!!!!!");

    SEQAN_ASSERT_EQ(metas[1], "READ1");
    SEQAN_ASSERT_EQ(seqs[1], "AAAAAAAAAA");
    SEQAN_ASSERT_EQ(quals[1], "!!!!!!!!!!");

    SEQAN_ASSERT_EQ(metas[2], "");
    SEQAN_ASSERT_EQ(seqs[2], "");
    SEQAN_ASSERT_EQ(quals[2], "");
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_sam_read_sequences_and_qualities)
{
    testSeqIOBamFileReadSequencesAndQualities("/tests/seq_io/small_sequences.sam");
}

SEQAN_DEFINE_TEST(test_seq_io_bam_file_bam_read_sequences_and_qualities)
{
    testSeqIOBamFileReadSequencesAndQualities("/tests/seq_io/small_sequences.bam");
}

#endif  // TESTS_SEQ_IO_TEST_READ_BAM_H_
