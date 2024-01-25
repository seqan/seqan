// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2024, Knut Reinert, FU Berlin
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

#include <seqan/basic.h>
#include <seqan/stream.h>

#include <seqan/align.h>

// Call the internal interface since we cannot publicly control the Traceback configuration.
template <typename TGapSequence1, typename TGapSequence2, typename TScoringScheme>
auto compute_affine_alignment(TGapSequence1 & gapSequence1,
                              TGapSequence2 & gapSequence2,
                              TScoringScheme const & scoringScheme)
{
    using TSize = typename seqan2::Size<TGapSequence1>::Type;
    using TPosition = typename seqan2::Position<TGapSequence1>::Type;
    using TTraceSegment = seqan2::TraceSegment_<TPosition, TSize>;

    // Test explicitly the complete trace functionality.
    using TTracebackPolicy = seqan2::TracebackOn<seqan2::TracebackConfig_<seqan2::CompleteTrace, seqan2::GapsLeft>>;
    using TAlignConfig = seqan2::AlignConfig2<seqan2::LocalAlignment_<>,
                                             seqan2::DPBandConfig<seqan2::BandOff>,
                                             seqan2::FreeEndGaps_<seqan2::True, seqan2::True, seqan2::True, seqan2::True>,
                                             TTracebackPolicy>;

    seqan2::String<TTraceSegment> trace;
    seqan2::DPScoutState_<seqan2::Default> dpScoutState;
    TAlignConfig config;
    auto score = seqan2::_setUpAndRunAlignment(trace,
                                              dpScoutState,
                                              seqan2::source(gapSequence1),
                                              seqan2::source(gapSequence2),
                                              scoringScheme,
                                              config,
                                              seqan2::AffineGaps());
    seqan2::_adaptTraceSegmentsTo(gapSequence1, gapSequence2, trace);
    return score;
}

SEQAN_DEFINE_TEST(test_local_alignment_with_complete_trace_produces_wrong_alignment)
{
    seqan2::String<seqan2::AminoAcid> s1 = "QDPKEPCISIYCLYLPVAGLGNLRGCCLPWM**";
    seqan2::String<seqan2::AminoAcid> s2 = "*SS*IVSRTPARTSPPGPGPIGGPH*TIQSVVATVGVYKGQGRNQRELMTRAYWEFLVQGK*LQFPIPST"
                                         "TGVQRVSRTCRPRRAHADPVSVARVRPRTSKGITDLLLLNLVWLNATCPSKKLNADRDGRVTI*QARVSF"
                                         "VIGINQTNRSTN*ERPCTTTHRIKKELSICQSSLCPGRVRFPVLSQIKPQAPLLVVPFRQFL*VSALQPY"
                                         "FPRNPKTLVSRKLPEGSSM*RPPIASWHRL*SELGRYLIVFEPLTFALD*RKHSWQMLSQSFVLRRSKNF"
                                         "TSNGTVRIAPVCPS*SLPRAPKTNKIEPRSYSIIPCTIIQAREPALNTLIFSK*TFRPPPTLSQEHRRRT"
                                         "GGKARTSSTRLATDRPPAPKIQLRAF*PQQLKYILLELELPRLLAPDLPSNGYSLKLLKCTHSNYRASNE"
                                         "SCIVIFRHYLPASGMGNLRACCLPWMW*PFLRLPLRNRTLIPRHP*EPR*ANTVPTKVDRADT*MIRRRC"
                                         "*TVRSAKLSRVTKANAPEDAAGFWSDKCT";

    seqan2::Align<seqan2::String<seqan2::AminoAcid> > correctAlignment;
    seqan2::resize(seqan2::rows(correctAlignment), 2);
    seqan2::assignSource(seqan2::row(correctAlignment, 0), s1);
    seqan2::assignSource(seqan2::row(correctAlignment, 1), s2);

    seqan2::Blosum62 scoringScheme(-1, -11);

    // Compute the known correct alignment.
    auto correctScore = seqan2::localAlignment(correctAlignment, scoringScheme);

    seqan2::Align<seqan2::String<seqan2::AminoAcid> > testAlignment;
    seqan2::resize(seqan2::rows(testAlignment), 2);
    seqan2::assignSource(seqan2::row(testAlignment, 0), s1);
    seqan2::assignSource(seqan2::row(testAlignment, 1), s2);

    // Compute the known incorrect alignment.
    auto testScore = compute_affine_alignment(seqan2::row(testAlignment, 0),
                                              seqan2::row(testAlignment, 1),
                                              scoringScheme);

    SEQAN_ASSERT_EQ(correctScore, testScore);
    SEQAN_ASSERT_EQ(correctAlignment, testAlignment);
}

SEQAN_BEGIN_TESTSUITE(test_align_simd_bugs)
{
    SEQAN_CALL_TEST(test_local_alignment_with_complete_trace_produces_wrong_alignment);
}
SEQAN_END_TESTSUITE
