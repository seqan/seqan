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

#ifndef TESTS_SEQ_IO_TEST_TAG_SELECT_INTERSECT_H_
#define TESTS_SEQ_IO_TEST_TAG_SELECT_INTERSECT_H_

// ---------------------------------------------------------------------------
// Assign tags to an output file based on the format of input file.
// ---------------------------------------------------------------------------

void testTransferTag()
{
    // in the case of no ZLIB BAM tag is not in the TagList and tagId is shifted by 1
    int offset = 1;
#if SEQAN_HAS_ZLIB
    offset = 0;
#endif

    SeqOutFormat outFormat;
    SeqFileIn inputFile;
    seqan::CharString filePath;
    bool result = false;

    // In the case of unselected input tag
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(outFormat.tagId, -1);
    close(inputFile);

    // In the case of unavailable tag
    filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/seq_io/test_genbank.gbk");
    open(inputFile, toCString(filePath));
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, false);
    SEQAN_ASSERT_EQ(outFormat.tagId, -1);
    close(inputFile);

    // In the case of selected input tag
    filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/seq_io/test_dna.fa");
    open(inputFile, toCString(filePath));
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(outFormat.tagId, 3 - offset);
    close(inputFile);

    filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/seq_io/test_dna.fq");
    open(inputFile, toCString(filePath));
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(outFormat.tagId, 4 - offset);
    close(inputFile);

    filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/seq_io/small_sequences.sam");
    open(inputFile, toCString(filePath));
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(outFormat.tagId, 1 - offset);
    close(inputFile);

#if SEQAN_HAS_ZLIB
    filePath = SEQAN_PATH_TO_ROOT();
    append(filePath, "/tests/seq_io/small_sequences.bam");
    open(inputFile, toCString(filePath));
    result = tagSelectIntersect(outFormat, inputFile.format);
    SEQAN_ASSERT_EQ(result, true);
    SEQAN_ASSERT_EQ(outFormat.tagId, 0);
    close(inputFile);
#endif
}

SEQAN_DEFINE_TEST(test_tag_select_intersect)
{
    testTransferTag();
}

#endif  // TESTS_SEQ_IO_TEST_TAG_SELECT_INTERSECT_H_
