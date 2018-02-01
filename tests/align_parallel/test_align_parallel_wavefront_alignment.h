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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#include <seqan/align_parallel.h>

namespace test_align_parallel
{
struct DPTestConfig : public seqan::DPTraits::GlobalAffine
{
    using TTracebackType = seqan::TracebackOff;
};
}  // namespace test_align_parallel

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_single_global_alignment)
{
    using namespace seqan;
    // We need to be able to construct a thread pool.

    // Define the traits object based on the settings.
    // now we want to define a WavefrontAlignmentTask connecting all the components together.
    DnaString seqH = "GGTTTTGTTTGATGGAGAATTGCGCAGAGGGGTTATATCTGCGTGAGGATCTGTCACTCGGCGGTGTGGG"
                     "ATACCTCCCTGCTAAGGCGGGTTGAGTGATGTTCCCTCGGACTGGGGACCGCTGGCTTGCGAGCTATGTC"
                     "CGCTACTCTCAGTACTACACTCTCATTTGAGCCCCCGCTCAGTTTGCTAGCAGAACCCGGCACATGGTTC"
                     "GCCGATACTATGGATTTTCTAAAGAAACACTCTGTTAGGTGGTATGAGTCATGACGCACGCAGGGAGAGG"
                     "CTAAGGCTTATGCTATGCTGATCTCCGTGAATGTCTATCATTCCTCTGCAGGACCC";

    DnaString seqV = "ACAGAGCGCGTACTGTCTGACGACGTATCCGCGCGGACTAGAAGGCTGGTGCCTCGTCCAACAAATAGAT"
                     "ACAGAAATCCACCGAAGTAAAGATCTCCAATTGTGGCACCACCAGGTGGCCACCACTCTTTGAAGTGAGG"
                     "AGACTTGCTTTACGTGTTTGTTCAGCCCGAGCTTTCGCTCGCACTGGAACACTGGTGTTTCGTCCTTTCG"
                     "GACTCATCAGTCAAGGTACGCACCTTGAGACACCGGGAAACAATCGATCAATCTTTCACAGAGCAACGAG"
                     "TTCGCTACTCTTGCAAAAGATCGACTTCCTATTTCGTGGATA";

    using TDPSettings = DPSettings<Score<int, Simple>, test_align_parallel::DPTestConfig>;
    TDPSettings settings;
    settings.scoringScheme = Score<int, Simple>{2, -2, -1, -11};

    WavefrontAlignmentTask<DnaString, DnaString, TDPSettings> task{seqH, seqV, settings, 37};

    using TThreadLocal = typename WavefrontAlignmentTaskConfig<TDPSettings>::TThreadLocal;

    EnumerableThreadLocal<TThreadLocal> tls{TThreadLocal{1}};

    WavefrontTaskScheduler scheduler(1, 1);
    lockWriting(scheduler);
    waitForWriters(scheduler);

    WavefrontAlignmentExecutor<WavefrontTaskScheduler, decltype(tls)> executor{&scheduler, &tls};

    int testScore{};
    task(0, executor, [&](auto const /*id*/, auto const & score)
    {
        testScore = score;
    });

    SEQAN_ASSERT_EQ(globalAlignmentScore(seqH, seqV, settings.scoringScheme, AlignConfig<false, false, false, false>()), testScore);
    unlockWriting(scheduler);
}

SEQAN_DEFINE_TEST(test_align_parallel_wavefront_multiple_global_alignment)
{
    using namespace seqan;

    DnaString seqH = "GGTTTTGTTTGATGGAGAATTGCGCAGAGGGGTTATATCTGCGTGAGGATCTGTCACTCGGCGGTGTGGG"
                     "ATACCTCCCTGCTAAGGCGGGTTGAGTGATGTTCCCTCGGACTGGGGACCGCTGGCTTGCGAGCTATGTC"
                     "CGCTACTCTCAGTACTACACTCTCATTTGAGCCCCCGCTCAGTTTGCTAGCAGAACCCGGCACATGGTTC"
                     "GCCGATACTATGGATTTTCTAAAGAAACACTCTGTTAGGTGGTATGAGTCATGACGCACGCAGGGAGAGG"
                     "CTAAGGCTTATGCTATGCTGATCTCCGTGAATGTCTATCATTCCTCTGCAGGACCC";

    DnaString seqV = "ACAGAGCGCGTACTGTCTGACGACGTATCCGCGCGGACTAGAAGGCTGGTGCCTCGTCCAACAAATAGAT"
                     "ACAGAAATCCACCGAAGTAAAGATCTCCAATTGTGGCACCACCAGGTGGCCACCACTCTTTGAAGTGAGG"
                     "AGACTTGCTTTACGTGTTTGTTCAGCCCGAGCTTTCGCTCGCACTGGAACACTGGTGTTTCGTCCTTTCG"
                     "GACTCATCAGTCAAGGTACGCACCTTGAGACACCGGGAAACAATCGATCAATCTTTCACAGAGCAACGAG"
                     "TTCGCTACTCTTGCAAAAGATCGACTTCCTATTTCGTGGATA";

    StringSet<DnaString> setH;
    StringSet<DnaString> setV;

    for (unsigned i = 0; i < 100; ++i)
    {
        appendValue(setH, (i % 2 == 0) ? seqV : seqH);
        appendValue(setV, (i % 5 == 0) ? seqH : seqV);
    }

    ExecutionPolicy<WavefrontAlignment<>, Serial> execPolicy;
    setNumThreads(execPolicy, 4);
    setParallelAlignments(execPolicy, 8);
    setBlockSize(execPolicy, 56);

    using TDPSettings = DPSettings<Score<int, Simple>, test_align_parallel::DPTestConfig>;
    TDPSettings settings;
    settings.scoringScheme = Score<int, Simple>{2, -2, -1, -11};

    std::vector<int> alignScores(length(setH), std::numeric_limits<int>::min());

    impl::alignExecBatch(execPolicy, setH, setV, settings, [&](auto const id, auto const score)
    {
        alignScores[id] = score;
    });

    for (unsigned i = 0; i < length(setH); ++i)
    {
        SEQAN_ASSERT_EQ(globalAlignmentScore(setH[i], setV[i], settings.scoringScheme, AlignConfig<false, false, false, false>()), alignScores[i]);
    }
}

#ifdef SEQAN_SIMD_ENABLED
SEQAN_DEFINE_TEST(test_align_parallel_wavefront_multiple_global_alignment_simd)
{
    using namespace seqan;

    DnaString seqH = "GGTTTTGTTTGATGGAGAATTGCGCAGAGGGGTTATATCTGCGTGAGGATCTGTCACTCGGCGGTGTGGG"
    "ATACCTCCCTGCTAAGGCGGGTTGAGTGATGTTCCCTCGGACTGGGGACCGCTGGCTTGCGAGCTATGTC"
    "CGCTACTCTCAGTACTACACTCTCATTTGAGCCCCCGCTCAGTTTGCTAGCAGAACCCGGCACATGGTTC"
    "GCCGATACTATGGATTTTCTAAAGAAACACTCTGTTAGGTGGTATGAGTCATGACGCACGCAGGGAGAGG"
    "CTAAGGCTTATGCTATGCTGATCTCCGTGAATGTCTATCATTCCTCTGCAGGACCC";

    DnaString seqV = "ACAGAGCGCGTACTGTCTGACGACGTATCCGCGCGGACTAGAAGGCTGGTGCCTCGTCCAACAAATAGAT"
    "ACAGAAATCCACCGAAGTAAAGATCTCCAATTGTGGCACCACCAGGTGGCCACCACTCTTTGAAGTGAGG"
    "AGACTTGCTTTACGTGTTTGTTCAGCCCGAGCTTTCGCTCGCACTGGAACACTGGTGTTTCGTCCTTTCG"
    "GACTCATCAGTCAAGGTACGCACCTTGAGACACCGGGAAACAATCGATCAATCTTTCACAGAGCAACGAG"
    "TTCGCTACTCTTGCAAAAGATCGACTTCCTATTTCGTGGATA";

    StringSet<DnaString> setH;
    StringSet<DnaString> setV;

    for (unsigned i = 0; i < 100; ++i)
    {
        appendValue(setH, (i % 2 == 0) ? seqV : seqH);
        appendValue(setV, (i % 5 == 0) ? seqH : seqV);
    }

    ExecutionPolicy<WavefrontAlignment<>, Vectorial> execPolicy;
    setNumThreads(execPolicy, 2);
    setParallelAlignments(execPolicy, 4);
    setBlockSize(execPolicy, 56);

    using TDPSettings = DPSettings<Score<int, Simple>, test_align_parallel::DPTestConfig>;
    TDPSettings settings;
    settings.scoringScheme = Score<int, Simple>{2, -2, -1, -11};

    std::vector<int> alignScores(length(setH), std::numeric_limits<int>::min());
    impl::alignExecBatch(execPolicy, setH, setV, settings, [&](auto const id, auto const score)
                         {
                             alignScores[id] = score;
                         });

    for (unsigned i = 0; i < length(setH); ++i)
    {
        SEQAN_ASSERT_EQ(globalAlignmentScore(setH[i], setV[i], settings.scoringScheme, AlignConfig<false, false, false, false>()), alignScores[i]);
    }
}
#endif
