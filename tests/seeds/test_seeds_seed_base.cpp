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
// Test the interface of the Seed class for specializations Simple Seed and
// Chained Seed.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

#include "seed_set_test_helpers.h"

template <typename TSeedSpec>
void testSeedsSeedBaseGettersSetters(TSeedSpec const &)
{
    using namespace seqan;

    // Define Seed type and declare a variable.
    typedef Seed<TSeedSpec> TSeed;
    TSeed s(1, 2, 3);

    // Check values from construction.
    SEQAN_ASSERT_EQ(1u, beginPositionH(s));
    SEQAN_ASSERT_EQ(4u, endPositionH(s));
    SEQAN_ASSERT_EQ(2u, beginPositionV(s));
    SEQAN_ASSERT_EQ(5u, endPositionV(s));
    SEQAN_ASSERT_EQ(-1, lowerDiagonal(s));
    SEQAN_ASSERT_EQ(-1, upperDiagonal(s));
    SEQAN_ASSERT_EQ(-1, beginDiagonal(s));
    SEQAN_ASSERT_EQ(-1, endDiagonal(s));

    // Use setters from base class.
    setLowerDiagonal(s, 5);
    SEQAN_ASSERT_EQ(5, lowerDiagonal(s));
    setUpperDiagonal(s, 42);
    SEQAN_ASSERT_EQ(42, upperDiagonal(s));
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
}

template <typename TSeedSpec>
void testSeedsSeedBaseAssign(TSeedSpec const &)
{
    using namespace seqan;
    
    typedef Seed<TSeedSpec> TSeed;
    TSeed seed(0, 0, 3);
    setScore(seed, -3);
    
    // Via copy constructor.
    {
        TSeed seed2(seed);
        SEQAN_ASSERT(seed2 == seed);
        SEQAN_ASSERT_EQ(score(seed2), score(seed));
    }
    // Via operator=.
    {
        TSeed seed2;
        seed2 = seed;
        SEQAN_ASSERT(seed2 == seed);
        SEQAN_ASSERT_EQ(score(seed2), score(seed));
    }
    // Via assign().
    {
        TSeed seed2;
        assign(seed2, seed);
        SEQAN_ASSERT(seed2 == seed);
        SEQAN_ASSERT_EQ(score(seed2), score(seed));
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

    typedef Seed<Simple> TSeed;

    typedef typename Position<TSeed const>::Type TPosition SEQAN_UNUSED_TYPEDEF;
    typedef typename Position<TSeed>::Type TPosition SEQAN_UNUSED_TYPEDEF;

    typedef typename Size<TSeed const>::Type TSize SEQAN_UNUSED_TYPEDEF;
    typedef typename Size<TSeed>::Type TSize SEQAN_UNUSED_TYPEDEF;

    typedef typename Diagonal<TSeed const >::Type TDiagonal SEQAN_UNUSED_TYPEDEF;
    typedef typename Diagonal<TSeed>::Type TDiagonal SEQAN_UNUSED_TYPEDEF;
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

SEQAN_BEGIN_TESTSUITE(test_seeds_seed_base)
{
    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_assign_simple);
    SEQAN_CALL_TEST(test_seeds_seed_base_constructors_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_metafunctions_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_getters_setters_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_basic_functions_chained);
    SEQAN_CALL_TEST(test_seeds_seed_base_assign_chained);
}
SEQAN_END_TESTSUITE
