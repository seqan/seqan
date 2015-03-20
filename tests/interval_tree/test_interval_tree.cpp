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
    seqan::IntervalTree<seqan::Interval<unsigned, int>, seqan::Static> intervalTree;

    SEQAN_ASSERT_EQ(0u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) == end(intervalTree, seqan::Standard()));
}

struct SmallIntervalTreeFixture
{
    typedef seqan::Interval<unsigned, int>                              TInterval;
    typedef seqan::IntervalTree<TInterval, seqan::Static>               TIntervalTree;
    typedef seqan::Iterator<TIntervalTree const, seqan::Standard>::Type TIterator;
    typedef std::vector<TIterator>                                      TResult;

    std::vector<TInterval> intervals;

    SmallIntervalTreeFixture()
    {
        intervals.push_back(TInterval(0, 1, 4));
        intervals.push_back(TInterval(1, 5, 9));
        intervals.push_back(TInterval(2, 4, 8));
        intervals.push_back(TInterval(3, 5, 7));
        intervals.push_back(TInterval(4, 16, 20));
        intervals.push_back(TInterval(5, 11, 16));
        intervals.push_back(TInterval(6, 30, 67));
    }
};

SEQAN_DEFINE_TEST(test_interval_tree_small_construction)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::Interval<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(7u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));
}

SEQAN_DEFINE_TEST(test_interval_tree_small_find_overlapping_with_point)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::Interval<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

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

SEQAN_DEFINE_TEST(test_interval_tree_small_find_overlapping_with_interval)
{
    SmallIntervalTreeFixture fixture;
    seqan::IntervalTree<seqan::Interval<unsigned, int>, seqan::Static> intervalTree(fixture.intervals);

    SEQAN_ASSERT_EQ(7u, length(intervalTree));
    SEQAN_ASSERT(begin(intervalTree, seqan::Standard()) != end(intervalTree, seqan::Standard()));

    // overlaps(0, 1) => {}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 0, 1, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
    // overlaps(5, 6) => {1,2,3}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 5, 6, result);
        SEQAN_ASSERT_EQ(length(result), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[0]), 2u);
        SEQAN_ASSERT_EQ(cargo(*result[1]), 3u);
        SEQAN_ASSERT_EQ(cargo(*result[2]), 1u);
    }
    // overlaps(67, 68) => {}
    {
        SmallIntervalTreeFixture::TResult result;
        findOverlappingWithInterval(intervalTree, 67, 68, result);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }
}

SEQAN_BEGIN_TESTSUITE(test_interval_tree)
{
    SEQAN_CALL_TEST(test_interval_tree_empty);
    SEQAN_CALL_TEST(test_interval_tree_small_construction);
    SEQAN_CALL_TEST(test_interval_tree_small_find_overlapping_with_point);
    SEQAN_CALL_TEST(test_interval_tree_small_find_overlapping_with_interval);
}
SEQAN_END_TESTSUITE
