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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_FLOATING_POINT_SCORE_H_
#define CORE_TESTS_ALIGN_TEST_ALIGNMENT_FLOATING_POINT_SCORE_H_

#include <seqan/basic.h>
#include <seqan/align.h>

SEQAN_DEFINE_TEST(test_alignment_flaoting_point_score_float_linear_gaps)
{
    using namespace seqan;

    DnaString seq0 = "ACGTGACGGATCGACGGACTAGC";
    DnaString seq1 = "ACGTAGCGAGGGGATGAGACGTGAGCGACT";

    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);

    Score<float, Simple> score(2.3, -1.5, -2.1);
    float max = globalAlignment(align, score, AlignConfig<>());

    SEQAN_ASSERT_EQ(max, 20.2999935150146484375);
    SEQAN_ASSERT_EQ(row(align, 0), "ACGT---GACG-GAT-CGACG-GA-CTAGC-");
    SEQAN_ASSERT_EQ(row(align, 1), "ACGTAGCGAGGGGATGAGACGTGAGCGA-CT");
}

SEQAN_DEFINE_TEST(test_alignment_flaoting_point_score_float_affine_gaps)
{
    using namespace seqan;

    DnaString seq0 = "ACGTGACGGATCGACGGACTAGC";
    DnaString seq1 = "ACGTAGCGAGGGGATGAGACGTGAGCGACT";

    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);

    Score<float, Simple> score(2.3, -1.5, -1.1, -4.2);
    float max = globalAlignment(align, score, AlignConfig<>());

    SEQAN_ASSERT_EQ(max, 10.70000362396240234375);
    SEQAN_ASSERT_EQ(row(align, 0), "ACGT---GA-CGGAT-CGACG-GA-CTAGC");
    SEQAN_ASSERT_EQ(row(align, 1), "ACGTAGCGAGGGGATGAGACGTGAGCGACT");
}

SEQAN_DEFINE_TEST(test_alignment_flaoting_point_score_double_linear_gaps)
{
    using namespace seqan;

    DnaString seq0 = "TTTTTTTTTTGTCGGATACGTAGCACTAGC";
    DnaString seq1 = "ACGTAGCGAGGGGATGAGACCCCCCCCCCC";

    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);

    Score<double, Simple> score(2.3, -1.5, -2.1);
    double max = localAlignment(align, score);

    SEQAN_ASSERT_EQ(max, (double) 16.300000000000000710542735760100185871124267578125);
    SEQAN_ASSERT_EQ(row(align, 0), "ACGTAGC-A");
    SEQAN_ASSERT_EQ(row(align, 1), "ACGTAGCGA");
}

SEQAN_DEFINE_TEST(test_alignment_flaoting_point_score_double_affine_gaps)
{
    using namespace seqan;

    DnaString seq0 = "ACGTGACGGATCGACGGACTAGC";
    DnaString seq1 = "ACGTAGCGAGGGGATGAGACGTGAGCGACT";

    Align<DnaString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seq0);
    assignSource(row(align, 1), seq1);

    Score<double, Simple> score(2.3, -1.5, -1.1, -4.2);
    double max = globalAlignment(align, score, AlignConfig<true, true, true, true>());

    SEQAN_ASSERT_EQ(max, (double) 17.400000000000002131628207280300557613372802734375);
    SEQAN_ASSERT_EQ(row(align, 0), "-------------ACGTGACG-GATCGACGGACTAGC");
    SEQAN_ASSERT_EQ(row(align, 1), "ACGTAGCGAGGGGATGAGACGTGAGCGACT-------");
}


#endif  // #ifndef CORE_TESTS_ALIGN_TEST_ALIGNMENT_FLOATING_POINT_SCORE_H_
