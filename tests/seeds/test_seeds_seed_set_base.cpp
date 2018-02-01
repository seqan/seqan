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
// Test the class SeedSet.
// ==========================================================================

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/stream.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

// Test the container functions for the given Seed and SeedSet
// specialization.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetContainerFunctions(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    // Define SeedSet type and declare a variable.
    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    TSeedSet s;

    // Test length/begin/end with empty set.
    SEQAN_ASSERT(begin(s) == end(s));
    SEQAN_ASSERT_EQ(0u, length(s));
    {  // Test with const empty set.
        TSeedSet const & cs = s;
        SEQAN_ASSERT(begin(cs) == end(cs));
        SEQAN_ASSERT_EQ(0u, length(cs));
    }

    // Insert one element, test basic accessor functions.
    typedef typename Value<TSeedSet>::Type TSeed;
    addSeed(s, TSeed(1, 2, 3), Single());

    // Test length/begin/end/front/back.
    SEQAN_ASSERT_EQ(1u, length(s));
    {
        typedef typename Iterator<TSeedSet, Rooted>::Type TIterator; // TODO(holtgrew): Why explicit rooted necessary?
        TIterator it(begin(s));
        ++it;
        SEQAN_ASSERT(it == end(s));
    }
    SEQAN_ASSERT(TSeed(1, 2, 3) == *begin(s, Standard()));
    SEQAN_ASSERT(TSeed(1, 2, 3) == front(s));
    SEQAN_ASSERT(TSeed(1, 2, 3) == back(s));
    {  // Same tests with const seed set.
        TSeedSet const & cs = s;
        {
            typedef typename Iterator<TSeedSet const, Rooted>::Type TIterator; // TODO(holtgrew): Why explicit rooted necessary?
            TIterator it(begin(cs));
            ++it;
            SEQAN_ASSERT(it == end(cs));
        }
        SEQAN_ASSERT(TSeed(1, 2, 3) == *begin(cs, Standard()));
        SEQAN_ASSERT(TSeed(1, 2, 3) == front(cs));
        SEQAN_ASSERT(TSeed(1, 2, 3) == back(cs));
    }
}

// Test addSeed(..., Single) for the given Seed and SeedSet
// specialization.
//
// Case: No quality threshold.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}

// Test addSeed(..., Single) for the given Seed and SeedSet
// specialization.
//
// Case: Seed size threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 1);

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}

// Test addSeed(..., Single) for the given Seed and SeedSet
// specialization.
//
// Case: Seed size threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdNotReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 4);

    TSeed seed(3, 3, 3);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Single) for the given Seed and SeedSet
// specialization.
//
// Case: Seed score threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdReachedScore(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, -1);

    TSeed seed(3, 3, 3);
    setScore(seed, 1);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(seed, front(set));
}


// Test addSeed(..., Single) for the given Seed and SeedSet
// specialization.
//
// Case: Seed score threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSingleThresholdNotReachedScore(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, -1);

    TSeed seed(3, 3, 3);
    setScore(seed, -2);

    addSeed(set, seed, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Merge) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Merging possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(3, 3, 3), Single());
    // Add seed with maximal diagonal distance 1, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(2, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);

    // std::cout << front(set) << std::endl;
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(2u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(2u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionV(front(set)));
}


// Test addSeed(..., Merge) with the given Seed and SeedSet
// Specialization.  Case: Seed right of added; Merging possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeRightMergingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(2, 2, 3), Single());
    // Add seed with maximal diagonal distance 1, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(3, 3, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(2u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(2u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionV(front(set)));
}


// Test addSeed(..., Merge) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Merging impossible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingImpossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Merging not possible because of diagonal distance
        TSeedSet set;
        addSeed(set, TSeed(0, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(3, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 2, 3), front(set));
    }
    {  // Merging not possible because non-overlapping
        TSeedSet set;
        addSeed(set, TSeed(0, 0, 2), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // are within this diagonal distance but do not overlap.
        bool ret = addSeed(set, TSeed(3, 2, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 2), front(set));
    }
}


// Test addSeed(..., Merge) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Merging possible; Length
// quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 5);

    // Add a low-scoring seed.
    addSeed(set, TSeed(1, 1, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to merge into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);
    // The seed is merged into the first one but the first one does
    // not exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Merge) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Merging possible; Length
// quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 4);

    // Add a low-quality seed.
    addSeed(set, TSeed(1, 1, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to merge into a high-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);
    // The seed is merged into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(4u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(4u, endPositionV(front(set)));
}


// Test addSeed(..., Merge) with the given Seed and SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Score quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, -1);

    // Add a low-quality seed.
    TSeed s1(1, 1, 3);
    setScore(s1, -2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to merge into another low-quality seed.
    TSeed s2(0, 0, 2);
    setScore(s2, -1);
    bool ret = addSeed(set, s2, 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);
    // The seed is merged with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Merge) with the given Seed and SeedSet Specialization.
// Seeds have scores.  Case: Seed left of added; Merging possible;
// Score quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, -2);

    // Add a low-quality seed.
    TSeed s1(1, 1, 4);
    setScore(s1, -3);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to merge into a high-quality seed.
    TSeed s2(0, 0, 2);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), Merge());
    SEQAN_ASSERT(ret);
    // The seed is merged with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionV(front(set)));
    SEQAN_ASSERT_EQ(-2, score(front(set)));
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(4, 5, 3), Single());
    // Add seed with maximal distance 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(1, 1, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(7u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(1u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(8u, endPositionV(front(set)));
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Case: Seed right of added; Chaining possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainRightChainingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    addSeed(set, TSeed(1, 1, 3), Single());
    // Add seed with maximal distance 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(4, 5, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(7u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(1u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(8u, endPositionV(front(set)));
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining impossible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingImpossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Chaining not possible because of distance
        TSeedSet set;
        addSeed(set, TSeed(0, 0, 3), Single());
        // Add seed with maximal distance 1, the two seeds do not
        // overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(5, 5, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 3), front(set));
    }
    {  // Chaining not possible because overlapping
        TSeedSet set;
        addSeed(set, TSeed(1, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, the two seeds
        // are within this distance but do not overlap.
        bool ret = addSeed(set, TSeed(0, 0, 3), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(set));
    }
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible;
// Length quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 6);

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);
    // The seed is chained into the first one but the result does not
    // exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible;
// Length quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 5);

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, Nothing(), Score<int, Simple>()/*TODO(holtgrew): unnecessary */, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);
    // The seed is chained into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionV(front(set)));
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Seeds have scores.  Case: Seed left of added;
// Chaining possible; Score quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, 2);
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 1);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into another low-quality seed.
    TSeed s2(4, 4, 2);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, Nothing(), scoringScheme, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);
    // The seed is chained with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., SimpleChain) with the given Seed and SeedSet
// Specialization.  Seeds have scores.  Case: Seed left of added;
// Chaining possible; Score quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, 3);
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a seed of sufficiently high quality.
    TSeed s2(4, 4, 2);
    setScore(s2, 2);
    bool ret = addSeed(set, s2, 1, Nothing(), scoringScheme, Nothing(), Nothing(), SimpleChain());
    SEQAN_ASSERT(ret);
    // The seed is chained with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionV(front(set)));
    SEQAN_ASSERT_EQ(3, score(front(set)));
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    TSeedSet set;
    addSeed(set, TSeed(4, 5, 3), Single());
    // Add seed with maximal distance 1, bandwidth 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(1, 1, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(7u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(1u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(8u, endPositionV(front(set)));
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Case: Seed right of added; Chaining possible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosRightChainingPossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    TSeedSet set;
    addSeed(set, TSeed(1, 1, 3), Single());
    // Add seed with maximal distance 1, bandwidth 2, the two seeds are within this distance.
    bool ret = addSeed(set, TSeed(4, 5, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary!*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);

    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(1u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(7u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(1u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(8u, endPositionV(front(set)));
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining impossible; No
// quality threshold required.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingImpossibleNoThreshold(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    {  // Chaining not possible because of distance
        DnaString sequence0 = "CCCCCCCCCC";
        DnaString sequence1 = "CCCCCCCCCC";

        TSeedSet set;
        addSeed(set, TSeed(0, 0, 3), Single());
        // Add seed with maximal distance 1, bandwidth 2, the two
        // seeds do not overlap but are are not within this distance.
        bool ret = addSeed(set, TSeed(5, 5, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(0, 0, 3), front(set));
    }
    {  // Chaining not possible because overlapping
        DnaString sequence0 = "CCCCCCCCCC";
        DnaString sequence1 = "CCCCCCCCCC";

        TSeedSet set;
        addSeed(set, TSeed(1, 2, 3), Single());
        // Add seed with maximal diagonal distance 1, bandwidth 2, the
        // two seeds are within this distance but do not overlap.
        bool ret = addSeed(set, TSeed(0, 0, 3), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
        SEQAN_ASSERT_NOT(ret);

        SEQAN_ASSERT_EQ(1u, length(set));
        SEQAN_ASSERT_EQ(TSeed(1, 2, 3), front(set));
    }
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible;
// Length quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 6);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary*/, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);
    // The seed is chained into the first one but the result does not
    // exceed the quality threshold yet.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Case: Seed left of added; Chaining possible;
// Length quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedLength(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinSeedSize(set, 5);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    addSeed(set, TSeed(2, 2, 3), Single());
    // The seed is added to the set but not included in the set of
    // seeds above the quality threshold.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a low-quality seed.
    bool ret = addSeed(set, TSeed(0, 0, 2), 1, 2, Score<int, Simple>()/*TODO(holtgrew): unnecessary */, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);
    // The seed is chained into the first one and exceeds the quality
    // threshold.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(5u, endPositionV(front(set)));
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Seeds have scores.  Case: Seed left of added;
// Chaining possible; Score quality threshold not reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, 4);
    Score<int, Simple> scoringScheme(1, -1, -1);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 1);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into another low-quality seed.
    TSeed s2(4, 4, 2);
    setScore(s2, 1);
    bool ret = addSeed(set, s2, 1, 2, scoringScheme, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);
    // The seed is chained with the existing low-quality seed into one
    // of low quality.
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.
}


// Test addSeed(..., Chaos) with the given Seed and SeedSet
// Specialization.  Seeds have scores.  Case: Seed left of added;
// Chaining possible; Score quality threshold reached.
template <typename TSeedSpec, typename TSeedSetSpec>
void testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedScored(TSeedSpec const &, TSeedSetSpec const &)
{
    using namespace seqan;

    typedef SeedSet<Seed<TSeedSpec>, TSeedSetSpec> TSeedSet;
    typedef typename Value<TSeedSet>::Type TSeed;

    TSeedSet set;
    setMinScore(set, 5);
    Score<int, Simple> scoringScheme(1, -1, -1);

    DnaString sequence0 = "CCCCCCCCCC";
    DnaString sequence1 = "CCCCCCCCCC";

    // Add a low-quality seed.
    TSeed s1(0, 0, 3);
    setScore(s1, 2);
    addSeed(set, s1, Single());
    SEQAN_ASSERT_EQ(1u, length(set));  // TODO(holtgrew): Build-in thresholds.

    // Add a seed to chain into a seed of sufficiently high quality.
    TSeed s2(4, 4, 2);
    setScore(s2, 2);
    bool ret = addSeed(set, s2, 1, 2, scoringScheme, sequence0, sequence1, Chaos());
    SEQAN_ASSERT(ret);
    // The seed is chained with the existing low-quality seed into one
    // of sufficiently high quality.
    SEQAN_ASSERT_EQ(1u, length(set));
    SEQAN_ASSERT_EQ(0u, beginPositionH(front(set)));
    SEQAN_ASSERT_EQ(0u, beginPositionV(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionH(front(set)));
    SEQAN_ASSERT_EQ(6u, endPositionV(front(set)));
    SEQAN_ASSERT_EQ(5, score(front(set)));
}

// Test container functions for specialization Simple Seed and
// Unordered SeedSet.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_functions_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetContainerFunctions(Simple(), Unordered());
}

// Test addSeed(..., Single) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: No threshold.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleNoThreshold(Simple(), Unordered());
}

// Test addSeed(..., Single) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedLength(Simple(), Unordered());
}

// Test addSeed(..., Single) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedLength(Simple(), Unordered());
}


// Test addSeed(..., Single) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedScore(Simple(), Unordered());
}


// Test addSeed(..., Single) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedScore(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeRightMergingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingImpossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedLength(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedLength(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(Simple(), Unordered());
}


// Test addSeed(..., Merge) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(Simple(), Unordered());
}

// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_right_chaining_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainRightChainingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_impossible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingImpossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedLength(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedLength(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedScored(Simple(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedScored(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosRightChainingPossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingImpossibleNoThreshold(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedLength(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedLength(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedScored(Simple(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Simple Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_simple_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedScored(Simple(), Unordered());
}


// Test container functions for specialization Chained Seed and
// Unordered SeedSet.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_container_functions_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetContainerFunctions(ChainedSeed(), Unordered());
}


// Test addSeed(..., Single) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: No threshold.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Single) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Single) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Single) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdReachedScore(ChainedSeed(), Unordered());
}


// Test addSeed(..., Single) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Size threshold, threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSingleThresholdNotReachedScore(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeRightMergingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingImpossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Merging not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdNotReachedScored(ChainedSeed(), Unordered());
}


// Test addSeed(..., Merge) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Merging is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedMergeLeftMergingPossibleThresholdReachedScored(ChainedSeed(), Unordered());
}

// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_right_chaining_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainRightChainingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_impossible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingImpossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdNotReachedScored(ChainedSeed(), Unordered());
}


// Test addSeed(..., SimpleChain) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedSimpleChainLeftChainingPossibleThresholdReachedScored(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is right of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosRightChainingPossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining is not possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingImpossibleNoThreshold(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.
//
// Case: Seed in set is left of added;  Chaining not possible;  Length quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedLength(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold not reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdNotReachedScored(ChainedSeed(), Unordered());
}


// Test addSeed(..., Chaos) for specialization Chained Seed and
// Unordered SeedSet.  Seeds have scores.
//
// Case: Seed in set is left of added;  Chaining is possible;  Quality threshold reached.
SEQAN_DEFINE_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_chained_unordered)
{
    using namespace seqan;
    testSeedsSeedSetAddSeedChaosLeftChainingPossibleThresholdReachedScored(ChainedSeed(), Unordered());
}

template <typename TSeed, typename TSeedSet>
void testSeedsSeedSetBaseClear()
{
    using namespace seqan;
    TSeedSet seedSet;

    addSeed(seedSet, TSeed(0, 0, 4), Single());
    addSeed(seedSet, TSeed(5, 5, 6), Single());

    SEQAN_ASSERT_EQ(length(seedSet), 2u);
    clear(seedSet);
    SEQAN_ASSERT_EQ(length(seedSet), 0u);
}

SEQAN_DEFINE_TEST(test_seeds_seed_set_base_clear_simple)
{
    using namespace seqan;
    testSeedsSeedSetBaseClear<Seed<Simple>, SeedSet<Seed<Simple> > >();
}

SEQAN_DEFINE_TEST(test_seeds_seed_set_base_clear_chained)
{
    using namespace seqan;
    testSeedsSeedSetBaseClear<Seed<ChainedSeed, DefaultSeedConfig>, SeedSet<Seed<ChainedSeed, DefaultSeedConfig> > >();
}

SEQAN_BEGIN_TESTSUITE(test_seeds_seed_set_base)
{
    // Tests for unordered seed sets and simple seeds.
    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_functions_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_clear_simple);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_clear_chained);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_right_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_simple_chain_left_chaining_possible_threshold_reached_scored_simple_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_simple_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_simple_unordered);

    // Tests for unordered seed sets and chained seeds.

    SEQAN_CALL_TEST(test_seeds_seed_set_base_container_functions_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_reached_score_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_single_threshold_not_reached_score_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_right_merging_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_merge_left_merging_possible_threshold_reached_scored_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_right_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chained_chain_left_chaining_possible_threshold_reached_scored_chained_unordered);

    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_right_chaining_possible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_impossible_no_threshold_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_length_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_not_reached_scored_chained_unordered);
    SEQAN_CALL_TEST(test_seeds_seed_set_base_add_seed_chaos_left_chaining_possible_threshold_reached_scored_chained_unordered);
}
SEQAN_END_TESTSUITE
