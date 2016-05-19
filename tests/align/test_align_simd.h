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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
#define TESTS_ALIGN_TEST_ALIGN_SIMD_H_

#include <seqan/basic.h>
#include <seqan/align.h>

// These are interface tests.

// How can we reduce redundant test instances?

// Global, Overlap, Semi-Global
// Local, affine
// Align, Graph, Gaps, Fragments, Score

struct LocalAlignTester_
{
    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return localAlignment(align, score);
        else
            return localAlignment(align, score, lDiag, uDiag);
    }

    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return localAlignment(align, score);
        else
            return localAlignment(align, score, lDiag, uDiag);

    }
};

struct LocalAlignScoreTester_
{
    template <typename TStringsH, typename TStringsV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TStringsH const & strH, TStringsV const & strV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &)
    {
        return localAlignmentScore(strH, strV, score);
    }

    template <typename TSeqH, typename TSeqV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TSeqH const & seqH, TSeqV const & seqV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const &)
    {
        return localAlignmentScore(seqH, seqV, score);
    }
};

struct GlobalAlignTester_
{
    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignment(align, score, config);
        else
            return globalAlignment(align, score, config, lDiag, uDiag);
    }

    template <typename TAlign, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TAlign & align, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignment(align, score, config);
        else
            return globalAlignment(align, score, config, lDiag, uDiag);
    }
};

struct GlobalAlignScoreTester_
{
    template <typename TStringsH, typename TStringsV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static seqan::String<TScoreValue>
    run(TStringsH const & strH, TStringsV const & strV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
        int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignmentScore(strH, strV, score, config);
        else
            return globalAlignmentScore(strH, strV, score, config, lDiag, uDiag);
    }

    template <typename TSeqH, typename TSeqV, typename TScoreValue, typename TScoreSpec, typename TConfig>
    static TScoreValue
    gold(TSeqH const & seqH, TSeqV const & seqV, seqan::Score<TScoreValue, TScoreSpec> const & score, TConfig const & config,
         int const lDiag, int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return globalAlignmentScore(seqH, seqV, score, config);
        else
            return globalAlignmentScore(seqH, seqV, score, config, lDiag, uDiag);
    }
};

template <typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TFunctor>
void testAlignSimd(TFunctor const &,
                   seqan::DnaString const & seqH, seqan::DnaString const & seqV,
                   seqan::Score<TScoreValue, TScoreSpec> const & score,
                   TAlignConfig const & config,
                   int const lDiag = seqan::MinValue<int>::VALUE,
                   int const uDiag = seqan::MaxValue<int>::VALUE)
{
    using namespace seqan;

    StringSet<Align<DnaString> > alignments;
    for(unsigned i = 0; i < 34; ++i)
    {
        Align<DnaString> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), seqH);
        assignSource(row(align, 1), seqV);
        appendValue(alignments, align);
    }

    String<TScoreValue> scores = TFunctor::run(alignments, score, config, lDiag, uDiag);

    SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));

    Align<DnaString> goldAlign;
    resize(rows(goldAlign), 2);
    assignSource(row(goldAlign, 0), seqH);
    assignSource(row(goldAlign, 1), seqV);

    TScoreValue goldScore = TFunctor::gold(goldAlign, score, config, lDiag, uDiag);

    for(size_t i = 0; i < 34; ++i)
    {
        SEQAN_ASSERT_EQ(scores[i], goldScore);
        SEQAN_ASSERT(row(alignments[i], 0) == row(goldAlign, 0));
        SEQAN_ASSERT(row(alignments[i], 1) == row(goldAlign, 1));
    }
}

template <typename TScoreValue, typename TScoreSpec, typename TAlignConfig, typename TTester>
void testAlignSimdScore(TTester const &,
                        seqan::DnaString const & seqH, seqan::DnaString const & seqV,
                        seqan::Score<TScoreValue, TScoreSpec> const & score,
                        TAlignConfig const & config,
                        int const lDiag = seqan::MinValue<int>::VALUE,
                        int const uDiag = seqan::MaxValue<int>::VALUE)
{
    using namespace seqan;

    StringSet<DnaString> stringsH, stringsV;

    for(unsigned i = 0; i < 34; ++i)
    {
        appendValue(stringsH, seqH);
        appendValue(stringsV, seqV);
    }

    String<TScoreValue> scores = TTester::run(stringsH, stringsV, score, config, lDiag, uDiag);

    SEQAN_ASSERT_EQ(length(scores), static_cast<decltype(length(scores))>(34));


    TScoreValue goldScore = TTester::gold(seqH, seqV, score, config, lDiag, uDiag);

    for(size_t i = 0; i < 34; ++i)
        SEQAN_ASSERT_EQ(scores[i], goldScore);
}

// Problem is clearly that the result is different.
SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_linear)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps.
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in one row.
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows.
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Blosum30(-1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_global_affine)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with no overlap - Simd version
    testAlignSimd(GlobalAlignTester_(), "GGGGGGGGG", "CCCCCCCCC", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_linear)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with no overlap - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "GGGGGGGGG", "CCCCCCCCC", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_overlap_affine)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_linear)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_gaps_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_score_semi_global_affine)
{
    using namespace seqan;

    AlignConfig<true, false, false, true> alignConfig;
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "ATGT", "ATAGAT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading gaps - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAGGGGTTTT", "AAAGTT", Score<int, Simple>(2, -1, -1, -3), alignConfig);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimdScore(GlobalAlignScoreTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", Score<int, Simple>(2, -1, -1, -3), alignConfig);
}

// ----------------------------------------------------------------------------
// Global Alignments Banded.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_global_linear_banded)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    Score<int, Simple> score(2, -1, -1);
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", score, alignConfig, -3, 2);
    // Alignment with both leading and trailing gaps in one row - Simd Version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", score, alignConfig, -2, 8);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", score, alignConfig, -4, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_global_affine_banded)
{
    using namespace seqan;

    AlignConfig<> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -3);
    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -3, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", scoringScheme, alignConfig, -2, 8);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -4, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_overlap_linear_banded)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1);
    // Simple alignment without any leading or trailing gaps.
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "GGGGG", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTTTT", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -2, 2);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd(GlobalAlignTester_(), "AACGCATTTTT", "TTTACGCA", scoringScheme, alignConfig, -2, 2);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd(GlobalAlignTester_(), "AACGCA", "TTTACGCA", scoringScheme, alignConfig, -2, 4);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd(GlobalAlignTester_(), "ACGAGTGTTTGCC", "TTTTTACGA", scoringScheme, alignConfig, -5, 7);
    // Alignment that starts at first position where the band crosses the bottom of the matrix - Simd version
    testAlignSimd(GlobalAlignTester_(), "ACGA", "TTTTTACGA", scoringScheme, alignConfig, -5, 4);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_overlap_affine_banded)
{
    using namespace seqan;

    AlignConfig<true, true, true, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -3);

    // Simple alignment without any leading or trailing gaps - Simd version
    testAlignSimd(GlobalAlignTester_(), "ATGT", "ATAGAT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in one row - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAGGGGTTTT", "AAAGTT", scoringScheme, alignConfig, -2, 2);
    // Alignment with both leading and trailing gaps in different rows - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTGGG", "TTTTTTTTGGGGGGGG", scoringScheme, alignConfig, -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_semi_global_linear_banded)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1);
    // More or less simple alignment - Simd version
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTTTTG", "AATTTTTTTTTTGGGGG", scoringScheme, alignConfig, -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_semi_global_affine_banded)
{
    using namespace seqan;

    AlignConfig<false, true, false, true> alignConfig;
    Score<int, Simple> scoringScheme(2, -1, -1, -4);
    testAlignSimd(GlobalAlignTester_(), "AAAAAATTTTTTTTG", "AATTTTTTTTTTGGGGG", scoringScheme, alignConfig, -2, 2);
}

// ----------------------------------------------------------------------------
// Local Alignments
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear)
{
    testAlignSimd(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", seqan::SimpleScore(2, -1, -2), seqan::Nothing());
    testAlignSimd(LocalAlignTester_(), "GGGGGGGGG", "CCCCCCCCC", seqan::SimpleScore(2, -1, -2), seqan::Nothing());
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine)
{
    testAlignSimd(LocalAlignTester_(), "CACACTTAACTTCACAA", "GGGGCTTGAGAGCTTGGGG", seqan::SimpleScore(2, -1, -1, -3), seqan::Nothing());
}

// ----------------------------------------------------------------------------
// Local Alignments Banded
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_linear_banded)
{
    using namespace seqan;

    testAlignSimd(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", SimpleScore(2, -1, -2, -2), Nothing(), -2, 2);
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_align_local_affine_banded)
{
    using namespace seqan;

        testAlignSimd(LocalAlignTester_(), "GGGGCTTAAGCTTGGGG", "AAAACTTAGCTCTAAAA", SimpleScore(2, -1, -2, -4), Nothing(), -2, 2);
}

#endif  // #ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
