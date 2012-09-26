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
// Test the interface of the Seed class for specializations Simple Seed and
// Chained Seed.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_
#define TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds2.h>  // Include module under test.

template <typename TSeedSpec>
void testSeedsSeedBaseConstructors(TSeedSpec const &)
{
    using namespace seqan;

    // Execute default constructor.
    {
        Seed<TSeedSpec> s;
    }
    // Execute with start position and length.
    {
        Seed<TSeedSpec> s(1, 2, 3);
    }
}

template <typename TSeedSpec>
void testSeedsSeedBaseGettersSetters(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    // Check values from construction.
    SEQAN_ASSERT_EQ(1u, getBeginDim0(s));
    SEQAN_ASSERT_EQ(4u, getEndDim0(s));
    SEQAN_ASSERT_EQ(2u, getBeginDim1(s));
    SEQAN_ASSERT_EQ(5u, getEndDim1(s));
    SEQAN_ASSERT_EQ(1, getLowerDiagonal(s));
    SEQAN_ASSERT_EQ(1, getUpperDiagonal(s));
    SEQAN_ASSERT_EQ(1, getStartDiagonal(s));
    SEQAN_ASSERT_EQ(1, getEndDiagonal(s));

    // Use setters from base class.
    setLowerDiagonal(s, 42);
    SEQAN_ASSERT_EQ(42, getLowerDiagonal(s));
    setUpperDiagonal(s, 5);
    SEQAN_ASSERT_EQ(5, getUpperDiagonal(s));
}

template <typename TSeedSpec>
void testSeedsSeedBaseBasicFunctions(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    {  // test assign
        TSeed x;
        assign(x, s);
        SEQAN_ASSERT(x == s);
    }
    {  // test move
        TSeed s2(s);
        TSeed x;
        move(x, s);
        SEQAN_ASSERT(x == s2);
    }
}

template <typename TSeedSpec>
void testSeedsSeedBaseAssign(TSeedSpec const &)
{
    using namespace seqan;
    
    // Test with score.
    {
        typedef Seed<TSeedSpec, DefaultSeedConfigScore> TSeed;
        TSeed seed(0, 0, 3);
        setScore(seed, -3);

        // Via copy constructor.
        {
            TSeed seed2(seed);
            SEQAN_ASSERT(seed2 == seed);
            SEQAN_ASSERT_EQ(getScore(seed2), getScore(seed));
        }
        // Via operator=.
        {
            TSeed seed2;
            seed2 = seed;
            SEQAN_ASSERT(seed2 == seed);
            SEQAN_ASSERT_EQ(getScore(seed2), getScore(seed));
        }
        // Via assign().
        {
            TSeed seed2;
            assign(seed2, seed);
            SEQAN_ASSERT(seed2 == seed);
            SEQAN_ASSERT_EQ(getScore(seed2), getScore(seed));
        }
    }
    // Test without score.
    {
        typedef Seed<TSeedSpec, DefaultSeedConfig> TSeed;
        TSeed seed(0, 0, 3);

        // Via copy constructor.
        {
            TSeed seed2(seed);
            SEQAN_ASSERT(seed2 == seed);
        }
        // Via operator=.
        {
            TSeed seed2;
            seed2 = seed;
            SEQAN_ASSERT(seed2 == seed);
        }
        // Via assign().
        {
            TSeed seed2;
            assign(seed2, seed);
            SEQAN_ASSERT(seed2 == seed);
        }
    }
}

// Test constructors of the SimpleSeed specialization, as specified
// for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_constructors_simple)
{
    using namespace seqan;
    testSeedsSeedBaseConstructors(Simple());
}

// Test the metafunctions of the SimpleSeed specialization, as
// specified for the base class Seed (none).
SEQAN_DEFINE_TEST(test_seeds_seed_base_metafunctions_simple)
{
    using namespace seqan;
    // No metafunctions in base, intentionally left blank.
}

// Test the getters and seeters of the SimpleSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_getters_setters_simple)
{
    using namespace seqan;
    testSeedsSeedBaseGettersSetters(Simple());
}

// Test the basic functions of the SimpleSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_basic_functions_simple)
{
    using namespace seqan;
    testSeedsSeedBaseBasicFunctions(Simple());
}

// Test the assign function for the SimpleSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_base_assign_simple)
{
    using namespace seqan;
    testSeedsSeedBaseAssign(Simple());
}

// Test constructors of the ChainedSeed specialization, as specified
// for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_constructors_chained)
{
    using namespace seqan;
    testSeedsSeedBaseConstructors(ChainedSeed());
}

// Test the metafunctions of the ChainedSeed specialization, as
// specified for the base class Seed (none).
SEQAN_DEFINE_TEST(test_seeds_seed_base_metafunctions_chained)
{
    using namespace seqan;
    // No metafunctions in base, intentionally left blank.
}

// Test the getters and seeters of the ChainedSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_getters_setters_chained)
{
    using namespace seqan;
    testSeedsSeedBaseGettersSetters(ChainedSeed());
}

// Test the basic functions of the ChainedSeed specialization as
// specified for the base class Seed.
SEQAN_DEFINE_TEST(test_seeds_seed_base_basic_functions_chained)
{
    using namespace seqan;
    testSeedsSeedBaseBasicFunctions(ChainedSeed());
}

// Test the assign function for the ChainedSeed specialization.
SEQAN_DEFINE_TEST(test_seeds_seed_base_assign_chained)
{
    using namespace seqan;
    testSeedsSeedBaseAssign(ChainedSeed());
}

#endif  // TEST_SEEDS_TEST_SEEDS_SEED_BASE_H_
