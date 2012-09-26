// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef TEST_RANDOM_TEST_RANDOM_DISTS_H_
#define TEST_RANDOM_TEST_RANDOM_DISTS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

SEQAN_DEFINE_TEST(test_random_normal_constructors)
{
    using namespace seqan;
    Pdf<Normal> pdf(1, 1);
}

SEQAN_DEFINE_TEST(test_random_normal_pick)
{
    using namespace seqan;
    
    // mean = 1.0, stddev = 1.0
    {
        Rng<MersenneTwister> mt(42);
        Pdf<Normal> pdf(1, 1);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 1), 0.02);
    }
    // mean = -3.4, stddev = 0.3
    {
        Rng<MersenneTwister> mt(42);
        Pdf<Normal> pdf(-3.4, 0.3);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 + 3.4), 0.02);
    }
    // mean = 4, stddev = 0.1
    {
        Rng<MersenneTwister> mt(42);
        Pdf<Normal> pdf(4, 0.1);
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 4), 0.02);
    }
}

SEQAN_DEFINE_TEST(test_random_geometric_fair_coin_constructors)
{
    using namespace seqan;
    Pdf<GeometricFairCoin> pdf;
}

SEQAN_DEFINE_TEST(test_random_geometric_fair_coin_pick)
{
    using namespace seqan;

    {
        Rng<MersenneTwister> mt(42);
        Pdf<GeometricFairCoin> pdf;
        
        double sum = 0;
        for (unsigned i = 0; i < 100000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_IN_DELTA(sum / 100000.0, 1, 0.01);
    }
}

SEQAN_DEFINE_TEST(test_random_lognormal_constructors)
{
    using namespace seqan;

    {
        Pdf<LogNormal> pdf(1, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._sigma, 1);
    }
    {
        Pdf<LogNormal> pdf(1, 1, MuSigma());
        SEQAN_ASSERT_EQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_EQ(pdf._normalDist._sigma, 1);
    }
    {
        Pdf<LogNormal> pdf(1, 1, MeanStdDev());
        SEQAN_ASSERT_NEQ(pdf._normalDist._mu, 1);
        SEQAN_ASSERT_NEQ(pdf._normalDist._sigma, 1);
    }
}

SEQAN_DEFINE_TEST(test_random_lognormal_pick)
{
    using namespace seqan;

    // mean = 1.0, stddev = 0.2
    {
        Rng<MersenneTwister> mt(42);
        Pdf<LogNormal> pdf(1, 0.2, MeanStdDev());
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 1), 0.02);
    }
    // mean = 0.4, stddev = 0.1
    {
        Rng<MersenneTwister> mt(42);
        Pdf<LogNormal> pdf(0.4, 0.1, MeanStdDev());
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 0.4), 0.02);
    }
    // mean = 4, stddev = 0.1
    {
        Rng<MersenneTwister> mt(42);
        Pdf<LogNormal> pdf(4, 0.1, MeanStdDev());

        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i)
            sum += pickRandomNumber(mt, pdf);
        SEQAN_ASSERT_LT(fabs(sum / 10000 - 4), 0.02);
    }
}

SEQAN_DEFINE_TEST(test_random_uniform_int_constructors)
{
    using namespace seqan;

    Pdf<Uniform<int> > pdf(-10, 10);
}

SEQAN_DEFINE_TEST(test_random_uniform_int_pick)
{
    using namespace seqan;

    Rng<MersenneTwister> mt(42);
    Pdf<Uniform<int> > pdf(-10, 10);

    unsigned gt = 0;  // Greater than 0.

    int sum = 0;
    for (unsigned i = 0; i < 100000; ++i) {
        int x = pickRandomNumber(mt, pdf);
        sum += x;
        gt += x > 0;
        SEQAN_ASSERT_GEQ(x, -10);
        SEQAN_ASSERT_LEQ(x, 10);
    }

    SEQAN_ASSERT_GT(gt, 0u);
    SEQAN_ASSERT_LT(gt, 100000u);

    SEQAN_ASSERT_LEQ(fabs(sum / 100000.0), 0.03);
}

SEQAN_DEFINE_TEST(test_random_uniform_bool_pick)
{
    using namespace seqan;
	
    Rng<MersenneTwister> mt(42);
    Pdf<Uniform<bool> > pdf;
	
    unsigned gt = 0;  // Greater than 0.
	
    int sum = 0;
    for (unsigned i = 0; i < 100000; ++i) {
        int x = pickRandomNumber(mt, pdf);
        sum += x;
        gt += x > 0;
        SEQAN_ASSERT_GEQ(x, 0);
        SEQAN_ASSERT_LEQ(x, 1);
    }
	
    SEQAN_ASSERT_GT(gt, 0u);
    SEQAN_ASSERT_LT(gt, 100000u);
	
	SEQAN_ASSERT_IN_DELTA(sum / 100000.0, 0.5, 0.01);
}

SEQAN_DEFINE_TEST(test_random_uniform_double_constructors)
{
    using namespace seqan;

    Pdf<Uniform<int> > pdf(0, 1);
}

SEQAN_DEFINE_TEST(test_random_uniform_double_pick)
{
    using namespace seqan;

    // Important case: In [0, 1]
    {
        Rng<MersenneTwister> mt(42);
        Pdf<Uniform<double> > pdf(0, 1);

        unsigned gt = 0;  // Greater than 0.5

        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i) {
            double x = pickRandomNumber(mt, pdf);
            sum += x;
            gt += x > 0.5;
            SEQAN_ASSERT_GEQ(x, 0);
            SEQAN_ASSERT_LEQ(x, 1);
        }

        SEQAN_ASSERT_GT(gt, 0u);
        SEQAN_ASSERT_LT(gt, 10000u);

        SEQAN_ASSERT_LEQ(fabs(sum / 10000.0 - 0.5), 0.02);
    }
    // Try something else...: In [-20, 30]
    {
        Rng<MersenneTwister> mt(42);
        Pdf<Uniform<double> > pdf(-20, 30);
        
        unsigned gt = 0;  // Greater than 0.5
        
        double sum = 0;
        for (unsigned i = 0; i < 10000; ++i) {
            double x = pickRandomNumber(mt, pdf);
            sum += x;
            gt += x > 0.5;
            SEQAN_ASSERT_GEQ(x, -20);
            SEQAN_ASSERT_LEQ(x, 30);
        }
        
        std::cerr << sum << std::endl;
        
        SEQAN_ASSERT_GT(gt, 0u);
        SEQAN_ASSERT_LT(gt, 50000u);
        
        SEQAN_ASSERT_LEQ(fabs(sum / 10000.0 - 5.0), 0.05);
    }
}

#endif  // TEST_RANDOM_TEST_RANDOM_DISTS_H_
