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
// Author: Anne-Katrin Emde <emde@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Tests for misc/misc_interval_tree.h.
// ==========================================================================

#ifndef SEQAN_HEADER_TEST_MISC_INTERVAL_TREE_H
#define SEQAN_HEADER_TEST_MISC_INTERVAL_TREE_H

// TODO(holtgrew): Are the static_casts<> to TValue necessary?
// TODO(holtgrew): Split up large tests into one setup and multiple test functions.

#include <seqan/misc/interval_tree.h>  // Header under test.

#include <seqan/basic.h>

namespace seqan {




// Build an IntervalTree from randomly generated intervals.
// Query IntervalTree with randomly generated points and compare result with naive
// algorithm (walk through intervals and collect intervals containing the point)
template<typename TValue, typename TConstructSpec, typename TStoreSpec>
void IntervalTreeTest_Random(TValue minValue, TValue maxValue) {

    // Define some types to test.
    typedef int                                     TCargo;
    typedef IntervalTree<TValue,TCargo>             TIntervalTree;

    String<TValue> intervalBegins;
    String<TValue> intervalEnds;
    String<TCargo> intervalCargos;

    //generate random intervals
    unsigned numIntervals = rand() % 1000 +1;
    for(unsigned i = 0; i < numIntervals; ++i)
    {
        TValue iBegin = (TValue)((double)rand()/((double)RAND_MAX+1) * (maxValue-minValue)) + minValue;
        TValue iEnd = (TValue)((double)rand()/((double)RAND_MAX+1) * (maxValue-minValue)) + minValue;
        if(iEnd < iBegin)
        {
             TValue tmp = iEnd;
             iEnd = iBegin;
             iBegin = tmp;
        }
        appendValue(intervalBegins,iBegin);
        appendValue(intervalEnds,iEnd);
        appendValue(intervalCargos,i);
    }

    // create interval tree
    TIntervalTree itree(begin(intervalBegins),begin(intervalEnds),begin(intervalCargos),numIntervals);

    // generate query points and test
    unsigned numQueries = rand() % 10 + 1;
    for(unsigned i = 0; i < numQueries; ++i)
    {
        TValue query = (TValue)((double)rand()/((double)RAND_MAX+1) * (maxValue-minValue)) + minValue;
        String<TCargo> itreeResult;
        findIntervals(itreeResult,itree,query);
        String<TCargo> naiveResult;
        for(unsigned j = 0; j < numIntervals; ++j)
        {
            if(intervalBegins[j] <= query && intervalEnds[j] > query)
                appendValue(naiveResult,intervalCargos[j]);
        }
        std::sort(begin(itreeResult),end(itreeResult));
        std::sort(begin(naiveResult),end(naiveResult));
        SEQAN_ASSERT_EQ(length(itreeResult),length(naiveResult));
        for(unsigned j = 0; j < length(naiveResult); ++j)
            SEQAN_ASSERT_EQ(itreeResult[j],naiveResult[j]);
    }
}

// Build an IntervalTree from constant data.  The perform some queries
// and check the results.
template<typename TValue, typename TConstructSpec, typename TStoreSpec>
void IntervalTreeTest_GraphMap() {

    // Define some types to test.
    typedef int                                     TCargo;
    typedef IntervalAndCargo<TValue, TCargo>        TInterval;
    typedef IntervalTreeNode<TInterval, TStoreSpec> TNode;
    typedef Graph<>                                 TGraph;
    typedef String<TNode>                           TPropertyMap;

    // Constant test fixture data.
    const size_t kValuesLen = 10;
    const double kValues[kValuesLen][2] = {
        { 5.0, 20.0}, { 6.0, 13.0}, { 8.0, 10.0}, { 9.0, 18.0}, {15.0, 23.0},
        {21.0, 25.0}, {29.0, 42.0}, {20.0, 36.0}, {38.0, 42.0}, {36.0, 48.0}
    };

    // Fill intervals with the intervals from kValues and a cargo that
    // is equal to its index.
    String<TInterval> intervals;
    resize(intervals, kValuesLen);
    for (size_t i = 0; i < kValuesLen; ++i) {
        intervals[i].i1 = static_cast<TValue>(kValues[i][0]);
        intervals[i].i2 = static_cast<TValue>(kValues[i][1]);
        intervals[i].cargo = i;
    }

    // Construct interval tree graph and its corresponding property map.
    TGraph g;
    TPropertyMap pm;
    createIntervalTree(g, pm, intervals, static_cast<TValue>(26.0),
                       TConstructSpec());


    // Query for 37.0 and check for the expected result.
    {
        TValue query = static_cast<TValue>(37.0);
        String<TCargo> result;
        findIntervals(result, g, pm, query);

        SEQAN_ASSERT_EQ(length(result), 2u);
        SEQAN_ASSERT((result[0] == 9 && result[1] == 6) ||
                          (result[0] == 6 && result[1] == 9));
    }
}




// Build an IntervalTree from constant data.  The, perform extensive
// checks on its structure
template <typename TValue>
void IntervalTreeTest_TreeStructure() {
    typedef int TCargo;
    typedef IntervalAndCargo<TValue, TCargo> TInterval;

    typedef IntervalTreeNode<TInterval, StoreIntervals> TNode;
    typedef Graph< > TGraph;
    typedef String<TNode> TPropertyMap;


    // Constant test fixture data.
    const size_t kValuesLen = 8;
    const double kValues[kValuesLen][2] = {
        { 1.0,  5.0}, { 3.0,  8.0}, { 2.0, 10.0}, {21.0, 23.0}, {22.0, 29.0},
        {23.0, 25.0}, {25.0, 31.0}, {26.0, 35.0}
    };

    // Fill intervals with the intervals from kValues and a cargo that
    // is equal to its index.
    String<TInterval> intervals;
    resize(intervals, kValuesLen);
    for (size_t i = 0; i < kValuesLen; ++i) {
        intervals[i].i1 = static_cast<TValue>(kValues[i][0]);
        intervals[i].i2 = static_cast<TValue>(kValues[i][1]);
        intervals[i].cargo = i;
    }

    // Create interval tree from the intervals.
    TGraph g;
    TPropertyMap pm;
    createIntervalTree(g, pm, intervals, ComputeCenter());

    // Perform extensives checks on the interval tree's structure.
    {
        SEQAN_ASSERT_EQ(length(pm), 6u);
        SEQAN_ASSERT_EQ(pm[0].center, static_cast<TValue>(18.0));
        SEQAN_ASSERT_EQ(pm[1].center, static_cast<TValue>(5.5));
        SEQAN_ASSERT_EQ(pm[2].center, static_cast<TValue>(3.0));
        SEQAN_ASSERT_EQ(pm[3].center, static_cast<TValue>(28.0));
        SEQAN_ASSERT_EQ(pm[4].center, static_cast<TValue>(23.0));
        SEQAN_ASSERT_EQ(pm[5].center, static_cast<TValue>(22.0));

        SEQAN_ASSERT_EQ(length(pm[0].list1), 0u);
        SEQAN_ASSERT_EQ(length(pm[1].list2), 2u);
        SEQAN_ASSERT_EQ(length(pm[2].list1), 1u);
        SEQAN_ASSERT_EQ(length(pm[3].list2), 3u);
        SEQAN_ASSERT_EQ(length(pm[4].list1), 1u);
        SEQAN_ASSERT_EQ(length(pm[5].list2), 1u);

        SEQAN_ASSERT_EQ(leftBoundary(pm[1].list1[0]), static_cast<TValue>(2.0));
        SEQAN_ASSERT_EQ(leftBoundary(pm[1].list1[1]), static_cast<TValue>(3.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[1].list2[0]), static_cast<TValue>(10.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[1].list2[1]), static_cast<TValue>(8.0));
        SEQAN_ASSERT_EQ(cargo(pm[1].list1[0]), 2);
        SEQAN_ASSERT_EQ(cargo(pm[1].list2[0]), 2);
        SEQAN_ASSERT_EQ(cargo(pm[1].list1[1]), 1);
        SEQAN_ASSERT_EQ(cargo(pm[1].list2[1]), 1);

        SEQAN_ASSERT_EQ(getLeftBoundary(pm[1].list1[0]), static_cast<TValue>(2.0));
        SEQAN_ASSERT_EQ(getLeftBoundary(pm[1].list1[1]), static_cast<TValue>(3.0));
        SEQAN_ASSERT_EQ(getRightBoundary(pm[1].list2[0]), static_cast<TValue>(10.0));
        SEQAN_ASSERT_EQ(getRightBoundary(pm[1].list2[1]), static_cast<TValue>(8.0));
        SEQAN_ASSERT_EQ(getCargo(pm[1].list1[0]), 2);
        SEQAN_ASSERT_EQ(getCargo(pm[1].list2[0]), 2);
        SEQAN_ASSERT_EQ(getCargo(pm[1].list1[1]), 1);
        SEQAN_ASSERT_EQ(getCargo(pm[1].list2[1]), 1);

        SEQAN_ASSERT_EQ(leftBoundary(pm[2].list1[0]), static_cast<TValue>(1.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[2].list2[0]), static_cast<TValue>(5.0));
        SEQAN_ASSERT_EQ(cargo(pm[2].list1[0]), 0);
        SEQAN_ASSERT_EQ(getCargo(pm[2].list2[0]), 0);

        SEQAN_ASSERT_EQ(leftBoundary(pm[3].list1[0]), static_cast<TValue>(22.0));
        SEQAN_ASSERT_EQ(leftBoundary(pm[3].list1[1]), static_cast<TValue>(25.0));
        SEQAN_ASSERT_EQ(leftBoundary(pm[3].list1[2]), static_cast<TValue>(26.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[3].list2[0]), static_cast<TValue>(35.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[3].list2[1]), static_cast<TValue>(31.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[3].list2[2]), static_cast<TValue>(29.0));

        SEQAN_ASSERT_EQ(leftBoundary(pm[4].list1[0]), static_cast<TValue>(23.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[4].list2[0]), static_cast<TValue>(25.0));

        SEQAN_ASSERT_EQ(leftBoundary(pm[5].list1[0]), static_cast<TValue>(21.0));
        SEQAN_ASSERT_EQ(rightBoundary(pm[5].list2[0]), static_cast<TValue>(23.0));
    }

}




// Build an IntervalTree from constant data.  Then perform some interval queries
// and check the results.
template<typename TValue, typename TConstructSpec>
void IntervalTreeTest_FindIntervalsIntervals() {

    // Define some types to test.
    typedef int                                     TCargo;
    typedef IntervalAndCargo<TValue, TCargo>        TInterval;
    typedef IntervalTree<TValue, TCargo>            TIntervalTree;

    // Constant test fixture data.
    const size_t kValuesLen = 10;
    const double kValues[kValuesLen][2] = {
        { 5.0, 20.0}, { 6.0, 13.0}, { 8.0, 10.0}, { 9.0, 18.0}, {15.0, 23.0},
        {21.0, 25.0}, {29.0, 42.0}, {20.0, 36.0}, {38.0, 42.0}, {36.0, 48.0}
    };

    // Fill intervals with the intervals from kValues and a cargo that
    // is equal to its index.
    String<TInterval> intervals;
    resize(intervals, kValuesLen);
    for (size_t i = 0; i < kValuesLen; ++i) {
        intervals[i].i1 = static_cast<TValue>(kValues[i][0]);
        intervals[i].i2 = static_cast<TValue>(kValues[i][1]);
        intervals[i].cargo = i;
    }

    // Construct interval tree graph and its corresponding property map.
    TIntervalTree itree(intervals,TConstructSpec());

    // Query for interval [10.0,20.0) and check for the expected result.
    {
        String<TCargo> result;
        findIntervals(result, itree, (TValue)10.0, (TValue)20.0);

        // intervals 0, 1, 3, and 4 should be returned
        SEQAN_ASSERT_EQ(length(result), 4u);

        TCargo check = 0;
        TCargo found = 100;
        for(unsigned i = 0; i < length(result); ++i)
        {
            if(result[i]==check){ found = check; break;}
        }
        SEQAN_ASSERT(found == check);

        check = 1;
        found = 100;
        for(unsigned i = 0; i < length(result); ++i)
        {
            if(result[i]==check){ found = check; break;}
        }
        SEQAN_ASSERT(found == check);

        check = 3;
        found = 100;
        for(unsigned i = 0; i < length(result); ++i)
        {
            if(result[i]==check){ found = check; break;}
        }
        SEQAN_ASSERT(found == check);

        check = 4;
        found = 100;
        for(unsigned i = 0; i < length(result); ++i)
        {
            if(result[i]==check){ found = check; break;}
        }
        SEQAN_ASSERT(found == check);

    }



}


// Build an IntervalTree from constant data.  The, perform extensive
// checks on its structure and also test some query results against
// expected results.
template <typename TValue>
void IntervalTreeTest_FindNoInterval() {
    typedef int TCargo;
    typedef IntervalAndCargo<TValue, TCargo> TInterval;

    typedef IntervalTreeNode<TInterval, StoreIntervals> TNode;
    typedef Graph< > TGraph;
    typedef String<TNode> TPropertyMap;


    // Constant test fixture data.
    const size_t kValuesLen = 8;
    const double kValues[kValuesLen][2] = {
        { 1.0,  5.0}, { 3.0,  8.0}, { 2.0, 10.0}, {21.0, 23.0}, {22.0, 29.0},
        {23.0, 25.0}, {25.0, 31.0}, {26.0, 35.0}
    };

    // Fill intervals with the intervals from kValues and a cargo that
    // is equal to its index.
    String<TInterval> intervals;
    resize(intervals, kValuesLen);
    for (size_t i = 0; i < kValuesLen; ++i) {
        intervals[i].i1 = static_cast<TValue>(kValues[i][0]);
        intervals[i].i2 = static_cast<TValue>(kValues[i][1]);
        intervals[i].cargo = i;
    }

    // Create interval tree from the intervals.
    TGraph g;
    TPropertyMap pm;
    createIntervalTree(g, pm, intervals,ComputeCenter());

    // There should be no interval containing 13.0.
    {
        TValue query = static_cast<TValue>(13.0);
        String<TCargo> result;
        findIntervals(result, g, pm, query);
        SEQAN_ASSERT_EQ(length(result), 0u);
    }

}

// Build an IntervalTree from constant data. Test query result against
// expected result.
template <typename TValue>
void IntervalTreeTest_FindIntervalExcludeTouching() {
    typedef int TCargo;
    typedef IntervalAndCargo<TValue, TCargo> TInterval;

    typedef IntervalTreeNode<TInterval, StoreIntervals> TNode;
    typedef Graph< > TGraph;
    typedef String<TNode> TPropertyMap;


    // Constant test fixture data.
    const size_t kValuesLen = 8;
    const double kValues[kValuesLen][2] = {
        { 1.0,  5.0}, { 3.0,  8.0}, { 2.0, 10.0}, {21.0, 23.0}, {22.0, 29.0},
        {23.0, 25.0}, {25.0, 31.0}, {26.0, 35.0}
    };

    // Fill intervals with the intervals from kValues and a cargo that
    // is equal to its index.
    String<TInterval> intervals;
    resize(intervals, kValuesLen);
    for (size_t i = 0; i < kValuesLen; ++i) {
        intervals[i].i1 = static_cast<TValue>(kValues[i][0]);
        intervals[i].i2 = static_cast<TValue>(kValues[i][1]);
        intervals[i].cargo = i;
    }

    // Create interval tree from the intervals.
    TGraph g;
    TPropertyMap pm;
    createIntervalTree(g, pm, intervals,ComputeCenter());

    // Query for all interval containing 23.0, only one interval should be returned
    {
        TValue query = static_cast<TValue>(23.0);
        String<TCargo> result;
        findIntervalsExcludeTouching(result, g, pm, query);
        SEQAN_ASSERT_EQ(length(result), 1u);
        SEQAN_ASSERT_EQ(result[0], 4);
    }
}



// test interval tree using the IntervalTree class
template <typename TValue>
void IntervalTreeTest_IntervalTree() {
    typedef IntervalAndCargo<TValue, int> TInterval;
    typedef IntervalTree<TValue, int> TIntervalTree;
    typedef typename Cargo<TIntervalTree>::Type TCargo;

    // Define test fixture data.
    const size_t kIntervalCount = 5;
    const TValue kBeginValues[kIntervalCount] = {2, 5, 3, 7, 1};
    const TValue kEndValues[kIntervalCount] = {7, 8, 5, 10, 5};
    const TCargo kCargoValues[kIntervalCount] = {2, 7, 1, 5, 8};

    // Build begin/end value and intervals vectors.
    String<TValue> begins;
    String<TValue> ends;
    String<TInterval> intervals;
    String<TCargo> cargos;
    resize(begins, kIntervalCount);
    resize(ends, kIntervalCount);
    resize(intervals, kIntervalCount);
    resize(cargos, kIntervalCount);
    for (size_t i = 0; i < kIntervalCount; ++i) {
        begins[i] = kBeginValues[i];
        ends[i] = kEndValues[i];
        cargos[i] = kCargoValues[i];

        intervals[i].i1 = kBeginValues[i];
        intervals[i].i2 = kEndValues[i];
        intervals[i].cargo = i;
    }

    // Construct IntervalTree from intervals vector, perform query on
    // it.
    // TODO(holtgrew): Duplicate of test IntervalTreeTest?
    {
        TIntervalTree itree(intervals);
        TValue query = 7;
        String<TCargo> result;

        findIntervals(result, itree, query);
        SEQAN_ASSERT_EQ(length(result), 2u);
        SEQAN_ASSERT_EQ(result[0], 1);
        SEQAN_ASSERT_EQ(result[1], 3);
    }


}


// test interval tree using the IntervalTree class and using the begin/end iterators
template <typename TValue>
void IntervalTreeTest_IntervalTreeFromIterator() {
    typedef IntervalAndCargo<TValue, int> TInterval;
    typedef IntervalTree<TValue, int> TIntervalTree;
    typedef typename Cargo<TIntervalTree>::Type TCargo;

    // Define test fixture data.
    const size_t kIntervalCount = 5;
    const TValue kBeginValues[kIntervalCount] = {2, 5, 3, 7, 1};
    const TValue kEndValues[kIntervalCount] = {7, 8, 5, 10, 5};
    const TCargo kCargoValues[kIntervalCount] = {2, 7, 1, 5, 8};

    // Build begin/end value and intervals vectors.
    String<TValue> begins;
    String<TValue> ends;
    String<TInterval> intervals;
    String<TCargo> cargos;
    resize(begins, kIntervalCount);
    resize(ends, kIntervalCount);
    resize(intervals, kIntervalCount);
    resize(cargos, kIntervalCount);
    for (size_t i = 0; i < kIntervalCount; ++i) {
        begins[i] = kBeginValues[i];
        ends[i] = kEndValues[i];
        cargos[i] = kCargoValues[i];

        intervals[i].i1 = kBeginValues[i];
        intervals[i].i2 = kEndValues[i];
        intervals[i].cargo = i;
    }

    // Construct IntervalTree from begin/end iterators.  Then, perform
    // a query and check the result.
    {
        TIntervalTree itree(begin(begins), begin(ends), length(begins));
        TValue query = 7;
        String<TCargo> result;

        findIntervals(result, itree, query);
        SEQAN_ASSERT_EQ(length(result), 2u);
        SEQAN_ASSERT_EQ(result[0], 1);
        SEQAN_ASSERT_EQ(result[1], 3);
    }

}



template <typename TValue>
void IntervalTreeTest_NonFullLength() {
    typedef IntervalAndCargo<TValue, int> TInterval;
    typedef IntervalTree<TValue, int> TIntervalTree;
    typedef typename Cargo<TIntervalTree>::Type TCargo;

    // Define test fixture data.
    const size_t kIntervalCount = 5;
    const TValue kBeginValues[kIntervalCount] = {2, 5, 3, 7, 1};
    const TValue kEndValues[kIntervalCount] = {7, 8, 5, 10, 5};
    const TCargo kCargoValues[kIntervalCount] = {2, 7, 1, 5, 8};

    // Build begin/end value and intervals vectors.
    String<TValue> begins;
    String<TValue> ends;
    String<TInterval> intervals;
    String<TCargo> cargos;
    resize(begins, kIntervalCount);
    resize(ends, kIntervalCount);
    resize(intervals, kIntervalCount);
    resize(cargos, kIntervalCount);
    for (size_t i = 0; i < kIntervalCount; ++i) {
        begins[i] = kBeginValues[i];
        ends[i] = kEndValues[i];
        cargos[i] = kCargoValues[i];

        intervals[i].i1 = kBeginValues[i];
        intervals[i].i2 = kEndValues[i];
        intervals[i].cargo = i;
    }


    // Construct IntervalTree from begins/ends/cargo vectors of
    // non-full length.  Then, perform a query on it that should yield
    // no result and a query that should yield a result.
    {
        TIntervalTree itree(begin(begins), begin(ends), begin(cargos), 4);
        String<TCargo> result;

        findIntervals(result, itree, 1);
        SEQAN_ASSERT_EQ(length(result), 0u);

        findIntervals(result, itree, 4);
        SEQAN_ASSERT_EQ(length(result), 2u);
        SEQAN_ASSERT((result[0] == 2 && result[1] == 1) ||
                          (result[0] == 1 && result[1] == 2));
    }


}





template <typename TValue>
void IntervalTreeTest_AddInterval() {
    typedef IntervalAndCargo<TValue, int> TInterval;
    typedef IntervalTree<TValue, int> TIntervalTree;
    typedef typename Cargo<TIntervalTree>::Type TCargo;

    // Define test fixture data.
    const size_t kIntervalCount = 5;
    const TValue kBeginValues[kIntervalCount] = {2, 5, 3, 7, 1};
    const TValue kEndValues[kIntervalCount] = {7, 8, 5, 10, 5};
    const TCargo kCargoValues[kIntervalCount] = {2, 7, 1, 5, 8};

    // Build begin/end value and intervals vectors.
    String<TValue> begins;
    String<TValue> ends;
    String<TInterval> intervals;
    String<TCargo> cargos;
    resize(begins, kIntervalCount);
    resize(ends, kIntervalCount);
    resize(intervals, kIntervalCount);
    resize(cargos, kIntervalCount);
    for (size_t i = 0; i < kIntervalCount; ++i) {
        begins[i] = kBeginValues[i];
        ends[i] = kEndValues[i];
        cargos[i] = kCargoValues[i];

        intervals[i].i1 = kBeginValues[i];
        intervals[i].i2 = kEndValues[i];
        intervals[i].cargo = i;
    }


    // Add intervals to a tree, then perform some queries and check
    // the results.
    {
        TIntervalTree itree(begin(begins), begin(ends), begin(cargos), 4);
        String<TCargo> result;

        // Add the "missing" interval to the tree.
        addInterval(itree, intervals[4]);
        SEQAN_ASSERT_EQ(itree.interval_counter, 5u);

        // Add interval to interval tree.
        TInterval iv;
        iv.i1 = 100;
        iv.i2 = 100;
        iv.cargo = 100;
        addInterval(itree, iv);
        SEQAN_ASSERT_EQ(itree.interval_counter, 6u);

        findIntervals(result, itree, 100);
        SEQAN_ASSERT_EQ(length(result), 1u);
        SEQAN_ASSERT_EQ(result[0], 100);

        findIntervalsExcludeTouching(result, itree, 100);
        SEQAN_ASSERT_EQ(length(result), 0u);

        addInterval(itree, 44, 88, 100);
        SEQAN_ASSERT_EQ(itree.interval_counter, 7u);

        findIntervals(result, itree, 77);
        SEQAN_ASSERT_EQ(length(result), 1u);
        SEQAN_ASSERT_EQ(result[0], 100);

        addInterval(itree, 30, 50);
        SEQAN_ASSERT_EQ(itree.interval_counter, 8u);

        findIntervals(result, itree, 48);
        SEQAN_ASSERT_EQ(length(result), 2u);
        SEQAN_ASSERT((result[0] == 100 && result[1] == 7) ||
                          (result[1] == 100 && result[0] == 7));
    }
}


SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_QueryAtBoundary)
{
    typedef IntervalAndCargo<int, double> TInterval;
    typedef String<TInterval, Alloc<Exact> > TIntervalList;

    TIntervalList tmpL;
    appendValue(tmpL, TInterval(0,30,1.4));
    appendValue(tmpL, TInterval(30,40,2.2));
    appendValue(tmpL, TInterval(40,60,3.3));
    IntervalTree<int, double> tmpT(tmpL);
    String<double> resultstmp;
    findIntervals(resultstmp, tmpT, 20, 30);
    SEQAN_ASSERT_EQ(length(resultstmp), 1u);
    SEQAN_ASSERT_EQ(resultstmp[0], 1.4);
}


// Call testIntervalTree_IntervalTree with <int> parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_IntervalTree__int)
{
    IntervalTreeTest_IntervalTree<int>();
}

// Call testIntervalTree_IntervalTreeFromIterator with <int> parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_IntervalTreeFromIterator__int)
{
    IntervalTreeTest_IntervalTreeFromIterator<int>();
}

// Call testIntervalTree_NonFullLength with <int> parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_NonFullLength__int)
{
    IntervalTreeTest_NonFullLength<int>();
}

// Call testIntervalTree_AddInterval with <int> parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_AddInterval__int)
{
    IntervalTreeTest_AddInterval<int>();
}

// Call IntervalTreeTests with <int, ComputeCenter, StorePointsOnly>
// parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_TreeStructure__int)
{
    IntervalTreeTest_TreeStructure<int>();
}

// Call IntervalTreeTest__FindNoInterval with <int, ComputeCenter, StorePointsOnly>
// parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_FindNoInterval__int)
{
    IntervalTreeTest_FindNoInterval<int>();
}

// Call IntervalTreeTest_FindIntervalExcludeTouching with <int>
// parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_FindIntervalExcludeTouching__int)
{
    IntervalTreeTest_FindIntervalExcludeTouching<int>();
}

// Call IntervalTreeTest_GraphMap with <int, ComputeCenter, StoreIntervals>
// parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_GraphMap__int_ComputeCenter_StoreIntervals)
{
    IntervalTreeTest_GraphMap<int, ComputeCenter, StoreIntervals>();
}


// Call IntervalTreeTest_FindIntervalsIntervals with <int, ComputeCenter, StoreIntervals>
// parametrization.
SEQAN_DEFINE_TEST(Interval_Tree__IntervalTreeTest_FindIntervalsIntervals__int_ComputeCenter)
{
    IntervalTreeTest_FindIntervalsIntervals<int, ComputeCenter>();
}

}  // seqan

#endif
