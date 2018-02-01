// ==========================================================================
//                         test_evaluate_alignment.h
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

#ifndef TESTS_ALIGN_TEST_EVALUATE_ALIGNMENT_H_
#define TESTS_ALIGN_TEST_EVALUATE_ALIGNMENT_H_

SEQAN_DEFINE_TEST(test_align_compute_alignment_stats)
{
    seqan::Peptide subject =
            "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
            "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
            "PGDFGADAQGAMNKALELFRKDMASNYK";
    seqan::Peptide query =
            "MSLTKTERTIIVSMWAKISTQADTIGTETLERLFLSHPQTKTYFPHFDLHPGSA"
            "QLRAHGSKVVAAVGDAVKSIDDIGGALSKLSELHAYILRVDPVNFKLLSHCLLVTLAARF"
            "PADFTAEAHAAWDKFLSVTEKYR";

    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    setSource(row(align, 0), subject);
    setSource(row(align, 1), query);

    // Compute the alignment.
    seqan::Blosum62 scoringScheme(-1, -12);
    int score = globalAlignment(align, scoringScheme);
//     std::cerr << align;
    SEQAN_ASSERT_EQ(score, 159);

    // Compute alignment statistics.
    seqan::AlignmentStats stats;
    score = computeAlignmentStats(stats, align, scoringScheme);
    SEQAN_ASSERT_EQ(score, 159);

    SEQAN_ASSERT_EQ(stats.numGapOpens, 2u);
    SEQAN_ASSERT_EQ(stats.numGapExtensions, 9u);
    SEQAN_ASSERT_EQ(stats.numMatches, 41u);
    SEQAN_ASSERT_EQ(stats.numMismatches, 96u);
    SEQAN_ASSERT_EQ(stats.numPositiveScores, 69u);
    SEQAN_ASSERT_EQ(stats.numNegativeScores, 68u);
    SEQAN_ASSERT_EQ(stats.alignmentScore, 159);
}

#endif  // #ifndef TESTS_ALIGN_TEST_EVALUATE_ALIGNMENT_H_
