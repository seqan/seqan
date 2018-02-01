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

#ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_BAND_H_
#define SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_BAND_H_

#include <seqan/basic.h>
#include <seqan/align.h>

void testDPBandConfigOffBandSize()
{
    using namespace seqan;

    DPBandConfig<BandOff> band;

    SEQAN_ASSERT_EQ(bandSize(band), 0u);
}

void testDPBandConfigOnBandSize()
{
    using namespace seqan;

    SEQAN_ASSERT_EQ(bandSize(DPBandConfig<BandOn>(-3, 7)), 11u);
    SEQAN_ASSERT_EQ(bandSize(DPBandConfig<BandOn>(3, 7)), 5u);
    SEQAN_ASSERT_EQ(bandSize(DPBandConfig<BandOn>(-7, -3)), 5u);
    SEQAN_ASSERT_EQ(bandSize(DPBandConfig<BandOn>(-7, 0)), 8u);
}

SEQAN_DEFINE_TEST(test_dp_band_on_constructor)
{
    using namespace seqan;
    DPBandConfig<BandOn> dpBand;

    SEQAN_ASSERT_EQ(dpBand._lowerDiagonal, 0);
    SEQAN_ASSERT_EQ(dpBand._upperDiagonal, 0);

    DPBandConfig<BandOn> dpBand2(-2, 2);
    SEQAN_ASSERT_EQ(dpBand2._lowerDiagonal, -2);
    SEQAN_ASSERT_EQ(dpBand2._upperDiagonal, 2);

    DPBandConfig<BandOn> dpBand3(dpBand2);

    SEQAN_ASSERT_EQ(dpBand3._lowerDiagonal, -2);
    SEQAN_ASSERT_EQ(dpBand3._upperDiagonal, 2);
}

SEQAN_DEFINE_TEST(test_dp_band_on_lower_diagonal)
{
    using namespace seqan;
    DPBandConfig<BandOn> dpBand;
    SEQAN_ASSERT_EQ(lowerDiagonal(dpBand), 0);

    dpBand._lowerDiagonal = -10;
    SEQAN_ASSERT_EQ(lowerDiagonal(dpBand), -10);
}

SEQAN_DEFINE_TEST(test_dp_band_on_upper_diagonal)
{
    using namespace seqan;
    DPBandConfig<BandOn> dpBand;
    SEQAN_ASSERT_EQ(upperDiagonal(dpBand), 0);

    dpBand._upperDiagonal = 10;
    SEQAN_ASSERT_EQ(upperDiagonal(dpBand), 10);
}

SEQAN_DEFINE_TEST(test_dp_band_on_set_lower_diagonal)
{
    using namespace seqan;
    DPBandConfig<BandOn> dpBand;

    setLowerDiagonal(dpBand, -10);
    SEQAN_ASSERT_EQ(dpBand._lowerDiagonal, -10);
}

SEQAN_DEFINE_TEST(test_dp_band_on_set_upper_diagonal)
{
    using namespace seqan;
    DPBandConfig<BandOn> dpBand;

    setUpperDiagonal(dpBand, 10);
    SEQAN_ASSERT_EQ(dpBand._upperDiagonal, 10);
}

SEQAN_DEFINE_TEST(test_dp_band_off_band_size)
{
    testDPBandConfigOffBandSize();
}

SEQAN_DEFINE_TEST(test_dp_band_on_band_size)
{
    testDPBandConfigOnBandSize();
}

#endif  // #ifndef SANDBOX_RMAERKER_TESTS_ALIGN2_TEST_ALIGNMENT_DP_BAND_H_
