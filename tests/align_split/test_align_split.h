// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// Tests for the align_split module.
// ==========================================================================

#ifndef SEQAN_TESTS_ALIGN_SPLIT_TEST_ALIGN_SPLIT_H_
#define SEQAN_TESTS_ALIGN_SPLIT_TEST_ALIGN_SPLIT_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/align_split.h>

SEQAN_DEFINE_TEST(test_align_split_overlapping_reads_in_reference_align_unbanded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT";
    seqan::DnaString seqR =                         "GGCTAGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig);
    assignSource(row(alignL, 1), seqL);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig);
    assignSource(row(alignR, 1), seqR);

    int res = splitAlignment(alignL, alignR, score);
    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 24);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 24);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 46);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 46);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_contigs_in_reference_align_unbanded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCT"                                   "AGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT""TTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG";
    seqan::DnaString seqR =   "TTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGGGCT""AGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig);
    assignSource(row(alignL, 1), seqL);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig);
    assignSource(row(alignR, 1), seqR);

    int res = splitAlignment(alignL, alignR, score);
    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 68);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 68);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 90);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 90);
}

SEQAN_DEFINE_TEST(test_align_split_insertion_in_reference_align_unbanded)
{
    // Define the input sequences.  contig1 has an insertion with respect to contig1.
    seqan::DnaString contig1 = "AGCATTTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCTT";
    seqan::DnaString contig2 = "AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig1);
    assignSource(row(alignL, 1), contig2);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig1);
    assignSource(row(alignR, 1), contig2);

    int res = splitAlignment(alignL, alignR, score);
    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig1, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATTTTAGATAAGATAG"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAG"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 19);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 19);

    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCTT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 30);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 30);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 56);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 56);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_reads_in_reference_gaps_unbanded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT";
    seqan::DnaString seqR =                         "GGCTAGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig);
    assignSource(gapsVL, seqL);
    assignSource(gapsHR, contig);
    assignSource(gapsVR, seqR);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score);
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 24);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 24);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 46);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 46);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_contigs_in_reference_gaps_unbanded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCT"                                   "AGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT""TTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG";
    seqan::DnaString seqR =   "TTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGGGCT""AGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig);
    assignSource(gapsVL, seqL);
    assignSource(gapsHR, contig);
    assignSource(gapsVR, seqR);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score);
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 68);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 68);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 90);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 90);
}

SEQAN_DEFINE_TEST(test_align_split_insertion_in_reference_gaps_unbanded)
{
    // Define the input sequences.  contig1 has an insertion with respect to contig1.
    seqan::DnaString contig1 = "AGCATTTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCTT";
    seqan::DnaString contig2 = "AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig1);
    assignSource(gapsVL, contig2);
    assignSource(gapsHR, contig1);
    assignSource(gapsVR, contig2);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score);
    // clearClipping(gapsHL);
    // clearClipping(gapsVL);
    // clearClipping(gapsHR);
    // clearClipping(gapsVR);
    // std::cerr << "gapsHL\n"
    //           << gapsHL << "\n"
    //           << "gapsVL\n"
    //           << gapsVL << "\n"
    //           << "gapsHL\n"
    //           << gapsHL << "\n"
    //           << "gapsHR\n"
    //           << gapsHR << "\n";
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig1, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATTTTAGATAAGATAG"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAG"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 19);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 19);

    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCTT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 30);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 30);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 56);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 56);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_reads_in_reference_align_banded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT";
    seqan::DnaString seqR =                         "GGCTAGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig);
    assignSource(row(alignL, 1), seqL);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig);
    assignSource(row(alignR, 1), seqR);

    int res = splitAlignment(alignL, alignR, score, -10, 10);

    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;

    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 24);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 24);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 46);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 46);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_contigs_in_reference_align_banded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCT"                                   "AGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT""TTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG";
    seqan::DnaString seqR =   "TTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGGGCT""AGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig);
    assignSource(row(alignL, 1), seqL);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig);
    assignSource(row(alignR, 1), seqR);

    int res = splitAlignment(alignL, alignR, score, -10, 10);
    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 68);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 68);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 90);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 90);
}

SEQAN_DEFINE_TEST(test_align_split_insertion_in_reference_align_banded)
{
    // Define the input sequences.  contig1 has an insertion with respect to contig1.
    seqan::DnaString contig1 = "AGCATTTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCTT";
    seqan::DnaString contig2 = "AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT";

    // Variables for the alignments.
    seqan::Align<seqan::DnaString> alignL, alignR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), contig1);
    assignSource(row(alignL, 1), contig2);
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), contig1);
    assignSource(row(alignR, 1), contig2);

    int res = splitAlignment(alignL, alignR, score, -12, 0);
    // clearClipping(row(alignL, 0));
    // clearClipping(row(alignL, 1));
    // clearClipping(row(alignR, 0));
    // clearClipping(row(alignR, 1));
    // std::cerr << "alignL\n"
    //           << alignL
    //           << "--------------\n"
    //           << "alignR\n"
    //           << alignR;
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << row(alignL, 0) << row(alignR, 0);
    SEQAN_ASSERT_EQ(contig1, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATTTTAGATAAGATAG"), row(alignL, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAG"), row(alignL, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 0)), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignL, 1)), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 0)), 19);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignL, 1)), 19);

    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCTT"), row(alignR, 0));
    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCAT"), row(alignR, 1));
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 0)), 30);
    SEQAN_ASSERT_EQ(clippedBeginPosition(row(alignR, 1)), 30);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 0)), 56);
    SEQAN_ASSERT_EQ(clippedEndPosition(row(alignR, 1)), 56);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_reads_in_reference_gaps_banded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT";
    seqan::DnaString seqR =                         "GGCTAGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig);
    assignSource(gapsVL, seqL);
    assignSource(gapsHR, contig);
    assignSource(gapsVR, seqR);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score, -10, 10);
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 24);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 24);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 46);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 46);
}

SEQAN_DEFINE_TEST(test_align_split_overlapping_contigs_in_reference_gaps_banded)
{
    // Define the input sequences.
    seqan::DnaString contig = "AGCATGTTAGATAAGATAGCTGTGCT"                                   "AGTAGGCAGTCAGCGCCAT";
    seqan::DnaString seqL =   "AGCCTGTTAGATAAGATAGCTGTGGT""TTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG";
    seqan::DnaString seqR =   "TTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGGGCT""AGTAGGCAGTCAGCGACAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig);
    assignSource(gapsVL, seqL);
    assignSource(gapsHR, contig);
    assignSource(gapsVR, seqR);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score, -10, 10);
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAGCTGT"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCCTGTTAGATAAGATAGCTGT"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 23);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 23);

    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGCCAT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("GCTAGTAGGCAGTCAGCGACAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 68);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 68);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 90);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 90);
}

SEQAN_DEFINE_TEST(test_align_split_insertion_in_reference_gaps_banded)
{
    // Define the input sequences.  contig1 has an insertion with respect to contig1.
    seqan::DnaString contig1 = "AGCATTTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCTT";
    seqan::DnaString contig2 = "AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT";

    // Variables for the alignments.
    seqan::Gaps<seqan::DnaString> gapsHL, gapsVL, gapsHR, gapsVR;
    seqan::Score<int, seqan::Simple> score(0, -1, -1, -1);

    assignSource(gapsHL, contig1);
    assignSource(gapsVL, contig2);
    assignSource(gapsHR, contig1);
    assignSource(gapsVR, contig2);

    int res = splitAlignment(gapsHL, gapsVL, gapsHR, gapsVR, score, -12, 0);
    // clearClipping(gapsHL);
    // clearClipping(gapsVL);
    // clearClipping(gapsHR);
    // clearClipping(gapsVR);
    // std::cerr << "gapsHL\n"
    //           << gapsHL << "\n"
    //           << "gapsVL\n"
    //           << gapsVL << "\n"
    //           << "gapsHL\n"
    //           << gapsHL << "\n"
    //           << "gapsHR\n"
    //           << gapsHR << "\n";
    SEQAN_ASSERT_EQ(res, -2);

    std::stringstream contigS;
    contigS << gapsHL << gapsHR;
    SEQAN_ASSERT_EQ(contig1, contigS.str());

    SEQAN_ASSERT_EQ(seqan::CharString("AGCATTTTAGATAAGATAG"), gapsHL);
    SEQAN_ASSERT_EQ(seqan::CharString("AGCATGTTAGATAAGATAG"), gapsVL);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHL), 0);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVL), 0);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHL), 19);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVL), 19);

    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCTT"), gapsHR);
    SEQAN_ASSERT_EQ(seqan::CharString("CTGTGCTAGTAGGCAGTCAGCGCCAT"), gapsVR);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsHR), 30);
    SEQAN_ASSERT_EQ(clippedBeginPosition(gapsVR), 30);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsHR), 56);
    SEQAN_ASSERT_EQ(clippedEndPosition(gapsVR), 56);
}

SEQAN_DEFINE_TEST(test_align_split_issue_1679)
{
    using namespace seqan;

    DnaString refLeft  = "TTTTTTTTTTTTGAGCCGATTTTTTTT";
    DnaString refRight = "TTTTTTTTTTTTTTTTGGACCGTTTTTTTTTTTTTTTTTTTTTTT";
    DnaString read     = "GAGCCGA" "GGACCG";

    Gaps<DnaString> refGapsLeft;
    Gaps<DnaString> refGapsRight;
    Gaps<DnaString> readGapsLeft;
    Gaps<DnaString> readGapsRight;

    setSource(refGapsLeft, refLeft);
    setSource(refGapsRight, refRight);
    setSource(readGapsLeft, read);
    setSource(readGapsRight, read);

    Score<int> scoring(1, -3, -4, -5);

    int splitScore = splitAlignment(readGapsLeft, refGapsLeft, readGapsRight, refGapsRight, scoring,
                                    AlignConfig<false, true, true, false>());

    SEQAN_ASSERT_EQ(splitScore, 13);
    SEQAN_ASSERT(refGapsLeft == "TTTTTTTTTTTTGAGCCGA");
    SEQAN_ASSERT(readGapsLeft == "------------GAGCCGA");

    SEQAN_ASSERT(refGapsRight == "GGACCGTTTTTTTTTTTTTTTTTTTTTTT");
    SEQAN_ASSERT(readGapsRight == "GGACCG-----------------------");
}

#endif  // SEQAN_TESTS_ALIGN_SPLIT_TEST_ALIGN_SPLIT_H_
