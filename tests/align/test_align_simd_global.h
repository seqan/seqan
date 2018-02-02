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

#ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_GLOBAL_H_
#define TESTS_ALIGN_TEST_ALIGN_SIMD_GLOBAL_H_

#include "test_align_simd_base.h"

// ----------------------------------------------------------------------------
// Global Alignments.
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(SimdAlignTestCommon, Linear_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::GlobalAlignTester_(), seqan::Score<int>(2, -1, -1),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignTester_(), seqan::Blosum62(-2),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Linear_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Score<int>(2, -1, -1),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Blosum62(-2),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Affine_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::GlobalAlignTester_(), seqan::Score<int>(2, -1, -1, -3),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignTester_(), seqan::Blosum62(-2, -4),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Dynamic_Score_Matrix_Align)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimd<seqan::AminoAcid>(::impl::test_align_simd::GlobalAlignTester_(), dynScore,
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Affine_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Score<int>(2, -1, -1, -3),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Blosum62(-2, -4),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Dynamic_Score_Matrix)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimdScore<seqan::AminoAcid>(::impl::test_align_simd::GlobalAlignScoreTester_(), dynScore,
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

#endif  // #ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_GLOBAL_H_
