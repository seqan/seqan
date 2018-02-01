// ==========================================================================
//                 test_alignment_algorithms_band_position.h
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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_BAND_POSITION_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_BAND_POSITION_H_

#include <seqan/basic.h>
#include <seqan/align.h>

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case1)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 1: LowerDiagonal <= UpperDiagonal < -length(seqV)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 2,
                                                           -static_cast<int>(length(strV)) - 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        int testScore = +std::numeric_limits<int>::min();
        SEQAN_ASSERT_EQ(score, testScore);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 2,
                                                           -static_cast<int>(length(strV)) - 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        int testScore = +std::numeric_limits<int>::min();
        SEQAN_ASSERT_EQ(score, testScore);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case2)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 2: LowerDiagonal < UpperDiagonal = -length(seqV)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1,
                                                           -static_cast<int>(length(strV))), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        int testScore = +std::numeric_limits<int>::min();
        SEQAN_ASSERT_EQ(score, testScore);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1,
                                                           -static_cast<int>(length(strV))), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "---------AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case3)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 3: LowerDiagonal < -length(seqV) < UpperDiagonal < 0
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, -3), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        int testScore = +std::numeric_limits<int>::min();
        SEQAN_ASSERT_EQ(score, testScore);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, -3),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "--------AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA-----------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case4)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 4: LowerDiagonal < -length(seqV) < UpperDiagonal = 0
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, 0), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        int testScore = +std::numeric_limits<int>::min();
        SEQAN_ASSERT_EQ(score, testScore);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, 0), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAA--CGT-GCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case5)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 5: LowerDiagonal < -length(seqV) < UpperDiagonal = length(seqH) - length(seqV)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1,
                                                           static_cast<int>(length(strH)) -
                                                           static_cast<int>(length(strV))), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1,
                                                           static_cast<int>(length(strH)) -
                                                           static_cast<int>(length(strV))), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case6)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 6: LowerDiagonal < -length(seqV) < length(seqH) - length(seqV) < UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, 6),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, 6),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case7)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 7: LowerDiagonal < -length(seqV) < length(seqH) - length(seqV) < UpperDiagonal = length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, length(strH)),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, length(strH)),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case8)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 8: LowerDiagonal < -length(seqV) < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);
        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, length(strH) + 1),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)) - 1, length(strH) + 1),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case9)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 9: -length(seqV) = LowerDiagonal <  UpperDiagonal < 0
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), -3), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), -3), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "--------AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA-----------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case10)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;
    //012345678901
    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 10: -length(seqV) = LowerDiagonal < 0 = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), 0), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), 0), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAA--CGT-GCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case11)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 11: -length(seqV) = LowerDiagonal < 0 < length(seqH) - length(seqV) = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH) - length(strV)),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH) - length(strV)),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case12)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 12: -length(seqV) = LowerDiagonal < 0 < length(seqH) - length(seqV) <  UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), 6), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), 6), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case13)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;
    //012345678901
    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 13: -length(seqV) = LowerDiagonal < 0 < length(seqH) =  UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH)),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case14)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 14: -length(seqV) = LowerDiagonal < 0 < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH) + 1),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-static_cast<int>(length(strV)), length(strH) + 1),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case15)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 15: -length(seqV) < LowerDiagonal <  0 = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, 0), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, 0), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 6);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAA--CGT-GCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case16)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 16: -length(seqV) < LowerDiagonal < 0 < length(seqH) - length(seqV) = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH) - length(strV)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH) - length(strV)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case17)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 17: -length(seqV) < LowerDiagonal < 0 < length(seqH) - length(seqV) <  UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, 6), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, 6), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case18)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 18: -length(seqV) < LowerDiagonal < 0 < length(seqH) =  UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case19)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 19: -length(seqV) < LowerDiagonal < 0 < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(-3, length(strH) + 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case20)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 20: LowerDiagonal = 0 < length(seqH) - length(seqV) = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH) - length(strV)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH) - length(strV)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case21)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 21: LowerDiagonal = 0 < length(seqH) - length(seqV) <  UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, 6), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, 6), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case22)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 22:  LowerDiagonal = 0 < length(seqH) =  UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case23)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 23:  LowerDiagonal = 0 < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(0, length(strH) + 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 15);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAA-G-GC-TTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case24)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 24: 0 < LowerDiagonal = length(seqH) - length(seqV) <  UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), 6), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;

        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), 6), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 3);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "---AAAGGCTTA");
    }
}


SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case25)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 25: 0 < LowerDiagonal = length(seqH) - length(seqV) < length(seqH) = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), length(strH)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;

        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), length(strH)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 3);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "---AAAGGCTTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case26)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 26: 0 < LowerDiagonal = length(seqH) - length(seqV) < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) - length(strV), length(strH) + 1),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 3);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "---AAAGGCTTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case27)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 27: length(seqH) - length(seqV) < LowerDiagonal < length(seqH) = UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(6, length(strH)), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(6, length(strH)), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA--------");
        SEQAN_ASSERT_EQ(ssV.str(), "-----------AAAGGCTTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case28)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);

    // Case 28: length(seqH) - length(seqV) < LowerDiagonal < length(seqH) < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(6, length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(6, length(strH) + 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 2);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA--------");
        SEQAN_ASSERT_EQ(ssV.str(), "-----------AAAGGCTTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case29)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 29: length(seqH) = LowerDiagonal < UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH), length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH), length(strH) + 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 0);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA---------");
        SEQAN_ASSERT_EQ(ssV.str(), "------------AAAGGCTTA");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case30)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 30: length(seqH) < LowerDiagonal <= UpperDiagonal
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) + 1, length(strH) + 1), TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme,
                                      DPBandConfig<BandOn>(length(strH) + 1, length(strH) + 1), TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }
}

SEQAN_DEFINE_TEST(test_alignment_algorithms_band_position_case31)
{
    using namespace seqan;

    typedef DPProfile_<GlobalAlignment_<>, LinearGaps, TracebackOn<> > TDPProfile;
    typedef DPProfile_<GlobalAlignment_<FreeEndGaps_<True, True, True, True> >, LinearGaps, TracebackOn<> > TDPProfileOverlap;
    typedef DPContext<DPCell_<int, LinearGaps>, typename TraceBitMap_<>::Type> TDPContext;
    TDPContext dpContext;
    DPScoutState_<Default> scoutState;

    Dna5String strH = "AAACGTGCTTTA";
    Dna5String strV = "AAAGGCTTA";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), strH);
    assignSource(row(align, 1), strV);

    Score<int, Simple> scoringScheme(2, -1, -1);
    // Case 31: -length(seqV) < LowerDiagonal = UpperDiagonal < length(seqH)
    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(0, 0),
                                      TDPProfile());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, std::numeric_limits<int>::min());

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "");
        SEQAN_ASSERT_EQ(ssV.str(), "");
    }

    {
        String<TraceSegment_<unsigned, unsigned> > traces;
        clearGaps(align);

        int score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(0, 0),
                                      TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);

        SEQAN_ASSERT_EQ(score, 3);

        std::stringstream ssH;
        std::stringstream ssV;

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA---");


        clearGaps(align);
        clear(traces);

        score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(-3, -3),
                                  TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);


        SEQAN_ASSERT_EQ(score, -6);

        ssH.str("");
        ssV.str("");

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "---AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------");

        clearGaps(align);
        clear(traces);

        score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(6, 6),
                                  TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);


        SEQAN_ASSERT_EQ(score, -6);

        ssH.str("");
        ssV.str("");

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA---");
        SEQAN_ASSERT_EQ(ssV.str(), "------AAAGGCTTA");

        clearGaps(align);
        clear(traces);

        score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(-9, -9),
                                  TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);


        SEQAN_ASSERT_EQ(score, 0);

        ssH.str("");
        ssV.str("");

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "---------AAACGTGCTTTA");
        SEQAN_ASSERT_EQ(ssV.str(), "AAAGGCTTA------------");

        clearGaps(align);
        clear(traces);

        score = _computeAlignment(dpContext, traces, scoutState, strH, strV, scoringScheme, DPBandConfig<BandOn>(12, 12),
                                  TDPProfileOverlap());
        _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traces);


        SEQAN_ASSERT_EQ(score, 0);

        ssH.str("");
        ssV.str("");

        ssH << row(align, 0);
        ssV << row(align, 1);

        SEQAN_ASSERT_EQ(ssH.str(), "AAACGTGCTTTA---------");
        SEQAN_ASSERT_EQ(ssV.str(), "------------AAAGGCTTA");
    }
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_ALGORITHMS_BAND_POSITION_H_
