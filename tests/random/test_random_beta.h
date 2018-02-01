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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for the probability distribution code in seqan/random.
// ==========================================================================

#ifndef TEST_RANDOM_TEST_RANDOM_BETA_H_
#define TEST_RANDOM_TEST_RANDOM_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/random.h>

SEQAN_DEFINE_TEST(test_random_beta_constructors)
{
    using namespace seqan;

    {  // Default constructor.
        BetaDistribution<double> pdf;
        SEQAN_ASSERT_EQ(pdf.alpha(), 2.0);
        SEQAN_ASSERT_EQ(pdf.beta(), 2.0);
    }

    {  // C'tor with alpha and beta value.
        BetaDistribution<double> pdf(0.5, 0.3);
        SEQAN_ASSERT_EQ(pdf.alpha(), 0.5);
        SEQAN_ASSERT_EQ(pdf.beta(), 0.3);
    }

    {  // C'tor with param.
        using TParam = BetaDistribution<double>::param_type;
        BetaDistribution<double> pdf(TParam(0.5, 0.3));
        SEQAN_ASSERT_EQ(pdf.alpha(), 0.5);
        SEQAN_ASSERT_EQ(pdf.beta(), 0.3);
    }
}

SEQAN_DEFINE_TEST(test_random_beta_alpha)
{
    using namespace seqan;
    
    BetaDistribution<double> pdf(0.5, 0.3);
    SEQAN_ASSERT_EQ(pdf.alpha(), 0.5);
}

SEQAN_DEFINE_TEST(test_random_beta_beta)
{
    using namespace seqan;

    BetaDistribution<double> pdf(0.5, 0.3);
    SEQAN_ASSERT_EQ(pdf.beta(), 0.3);
}

SEQAN_DEFINE_TEST(test_random_beta_min)
{
    using namespace seqan;

    BetaDistribution<double> pdf(0.5, 0.3);
    SEQAN_ASSERT_EQ(pdf.min(), 0.0);
}

SEQAN_DEFINE_TEST(test_random_beta_max)
{
    using namespace seqan;

    BetaDistribution<double> pdf(0.5, 0.3);
    SEQAN_ASSERT_EQ(pdf.max(), 1.0);
}

SEQAN_DEFINE_TEST(test_random_beta_param)
{
    using namespace seqan;

    BetaDistribution<double> pdf(0.5, 0.3);
    SEQAN_ASSERT_EQ(pdf.param().alpha(), 0.5);
    SEQAN_ASSERT_EQ(pdf.param().beta(), 0.3);
}

SEQAN_DEFINE_TEST(test_random_beta_set_param)
{
    using namespace seqan;

    using TParam = BetaDistribution<double>::param_type;

    BetaDistribution<double> pdf;

    pdf.setParam(TParam(0.5, 0.3));
    SEQAN_ASSERT_EQ(pdf.alpha(), 0.5);
    SEQAN_ASSERT_EQ(pdf.beta(), 0.3);
}

SEQAN_DEFINE_TEST(test_random_beta_write)
{
    using namespace seqan;

    BetaDistribution<double> pdf(0.5, 0.3);

    std::stringstream os;
    os << pdf;
    SEQAN_ASSERT_EQ(CharString("~Beta(0.5,0.3)"), os.str());
}

SEQAN_DEFINE_TEST(test_random_beta_read)
{
    using namespace seqan;

    BetaDistribution<double> pdf;

    std::stringstream is("~Beta(0.5,0.3)");
    is >> pdf;
    SEQAN_ASSERT_EQ(pdf.alpha(), 0.5);
    SEQAN_ASSERT_EQ(pdf.beta(), 0.3);
}

SEQAN_DEFINE_TEST(test_random_beta_functor)
{
    using namespace seqan;

    // mean = 0.3, stddev = 0.2
    {
        std::mt19937 mt(42);
        BetaDistribution<double> pdf(cvtBetaDistParam(0.3, 0.2));

        double sum = 0;
        for (unsigned i = 0; i < 1000; ++i)
            sum += pdf(mt);
        SEQAN_ASSERT_IN_DELTA((sum / 1000), 0.3, 0.02);
    }
    // mean = 0.4, stddev = 0.1
    {
        std::mt19937 mt(42);
        BetaDistribution<double> pdf(cvtBetaDistParam(0.4, 0.1));

        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pdf(mt);
        SEQAN_ASSERT_IN_DELTA((sum / 10000), 0.4, 0.01);
    }
    // mean = 0.7, stddev = 0.1
    {
        std::mt19937 mt(42);
        BetaDistribution<double> pdf(cvtBetaDistParam(0.7, 0.1));

        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pdf(mt);
        SEQAN_ASSERT_IN_DELTA((sum / 10000), 0.7, 0.01);
    }
}

#endif  // TEST_RANDOM_TEST_RANDOM_DISTS_H_
