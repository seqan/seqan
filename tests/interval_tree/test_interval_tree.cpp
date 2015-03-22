// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#include <seqan/basic.h>
#include <seqan/file.h>

#include <seqan/interval_tree.h>

SEQAN_DEFINE_TEST(test_interval_tree_empty)
{
    seqan::IntervalTree<seqan::IntervalWithCargo<unsigned, int>, seqan::Static> intervalTree;

    SEQAN_ASSERT_EQ(0u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) == end(intervalTree, seqan::Standard()));
}

struct SmallIntervalTreeFixture
{
    typedef seqan::IntervalWithCargo<unsigned, int>                              TIntervalWithCargo;
    typedef seqan::IntervalTree<TIntervalWithCargo, seqan::Static>               TIntervalTree;
    typedef seqan::Iterator<TIntervalTree const, seqan::Standard>::Type TIterator;
    typedef std::vector<TIterator>                                      TResult;

    std::vector<TIntervalWithCargo> intervals;

    SmallIntervalTreeFixture()
    {
        intervals.push_back(TIntervalWithCargo(0, 1, 4));
        intervals.push_back(TIntervalWithCargo(1, 5, 9));
        intervals.push_back(TIntervalWithCargo(2, 4, 8));
        intervals.push_back(TIntervalWithCargo(3, 5, 7));
        intervals.push_back(TIntervalWithCargo(4, 16, 20));
        intervals.push_back(TIntervalWithCargo(5, 11, 16));
        intervals.push_back(TIntervalWithCargo(6, 30, 67));
    }
};

SEQAN_DEFINE_TEST(test_interval_tree_small_construction)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::IntervalWithCargo<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(7u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));
}

SEQAN_DEFINE_TEST(test_interval_tree_small_find_overlapping_with_point)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::IntervalWithCargo<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(7u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));

    // overlaps(0) => {}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithPoint(intervalTree, 0, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
    // overlaps(5) => {1,2,3}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithPoint(intervalTree, 5, result);
        SEQAN_ASSERT_EQ(length(result), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[0]), 2u);
        SEQAN_ASSERT_EQ(cargo(*result[1]), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[2]), 1u);
    }
    // overlaps(67) => {}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithPoint(intervalTree, 67, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
}

struct BugIntervalTreeFixture
{
    typedef seqan::IntervalWithCargo<unsigned, int>                              TIntervalWithCargo;
    typedef seqan::IntervalTree<TIntervalWithCargo, seqan::Static>               TIntervalTree;
    typedef seqan::Iterator<TIntervalTree const, seqan::Standard>::Type TIterator;
    typedef std::vector<TIterator>                                      TResult;

    std::vector<TIntervalWithCargo> intervals;

    BugIntervalTreeFixture()
    {
        intervals.push_back(TIntervalWithCargo(0, 2985989, 2985990));
        intervals.push_back(TIntervalWithCargo(1, 3102611, 3102757));
        intervals.push_back(TIntervalWithCargo(2, 3102881, 3102984));
        intervals.push_back(TIntervalWithCargo(3, 3103125, 3103127));
    }
};

SEQAN_DEFINE_TEST(test_interval_tree_small_find_overlapping_with_interval)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::IntervalWithCargo<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(7u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));

    // overlaps(0, 1) => {}
    {
        BugIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 0, 1, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
    // overlaps(5, 6) => {1,2,3}
    {
        BugIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 5, 6, result);
        SEQAN_ASSERT_EQ(length(result), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[0]), 2u);
        SEQAN_ASSERT_EQ(cargo(*result[1]), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[2]), 1u);
    }
    // overlaps(67, 68) => {}
    {
        BugIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 67, 68, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
}

SEQAN_DEFINE_TEST(test_interval_tree_small_find_overlapping_with_interval_bug)
{
    BugIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::IntervalWithCargo<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(4u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));

    // overlaps(0, 1) => {}
    {
        BugIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 3102981, 3103082, result);
        SEQAN_ASSERT_EQ(length(result), 1u);
         SEQAN_ASSERT_EQ(cargo(*result[0]), 2u);
    }
}

SEQAN_BEGIN_TESTSUITE(test_interval_tree)
{
    SEQAN_CALL_TEST(test_interval_tree_empty);
    SEQAN_CALL_TEST(test_interval_tree_small_construction);
    SEQAN_CALL_TEST(test_interval_tree_small_find_overlapping_with_point);
    SEQAN_CALL_TEST(test_interval_tree_small_find_overlapping_with_interval);
    SEQAN_CALL_TEST(test_interval_tree_small_find_overlapping_with_interval_bug);
}
SEQAN_END_TESTSUITE
