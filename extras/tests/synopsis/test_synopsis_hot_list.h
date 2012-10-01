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
// Tests for the hot list data structures and algorithms of the synopsis
// module.
// ==========================================================================

#ifndef TEST_SYNOPSIS_TEST_SYNOPSIS_HOT_LIST_H_
#define TEST_SYNOPSIS_TEST_SYNOPSIS_HOT_LIST_H_

#include <seqan/basic.h>
#include <seqan/synopsis.h>
#include <seqan/random.h>

// ==========================================================================
// Generic Test Code
// ==========================================================================

template <typename THotList>
void testGenericHotListTopK(THotList & hotList)
{
    using namespace seqan;

    String<int> items;
    for (int i = 0; i < 10; ++i)
        appendValue(items, 0);
    for (int i = 0; i < 100; ++i)
        appendValue(items, 1);
    for (int i = 0; i < 1000; ++i)
        appendValue(items, 2);
    for (int i = 3; i < 100; ++i)
        appendValue(items, i);
    
    Rng<> rng(1);
    shuffle(items, rng);

    clear(hotList);
    for (unsigned i = 0; i < length(items); ++i)
        registerItem(hotList, items[i]);

    typedef typename Size<THotList>::Type TSize;

    String<Triple<int, TSize, TSize> > result;
    getItems(result, hotList);
    SEQAN_ASSERT_GEQ(length(result), 3u);
    SEQAN_ASSERT_EQ(result[0].i1, 2);
    SEQAN_ASSERT_EQ(result[1].i1, 1);
    SEQAN_ASSERT_EQ(result[2].i1, 0);
}

template <typename THotList>
void testGenericHotListIceberg(THotList & hotList)
{
    using namespace seqan;

    String<int> items;
    for (int i = 0; i < 100; ++i)
        appendValue(items, 0);
    for (int i = 0; i < 110; ++i)
        appendValue(items, 1);
    for (int i = 0; i < 120; ++i)
        appendValue(items, 2);
    for (int i = 3; i < 10; ++i)
        appendValue(items, i);
    
    Rng<> rng(1);
    shuffle(items, rng);

    clear(hotList);
    for (unsigned i = 0; i < length(items); ++i)
        registerItem(hotList, items[i]);

    typedef typename Size<THotList>::Type TSize;

    String<Triple<int, TSize, TSize> > result;
    getItems(result, hotList);
    SEQAN_ASSERT_EQ(length(result), 3u);
    SEQAN_ASSERT_EQ(result[0].i1, 2);
    SEQAN_ASSERT_EQ(result[1].i1, 1);
    SEQAN_ASSERT_EQ(result[2].i1, 0);
}

// ==========================================================================
// Actual Tests
// ==========================================================================

// Detailed white-box test.
SEQAN_DEFINE_TEST(test_synopsis_frequency_hot_list_frequent)
{
    using namespace seqan;

    HotList<int, Frequent> hotList(3);
    SEQAN_ASSERT_EQ(registerItem(hotList, 0), 1u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 1), 1u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 2), 0u);  // clear all
    SEQAN_ASSERT_EQ(registerItem(hotList, 0), 1u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 1), 1u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 0), 2u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 1), 2u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 2), 0u);  // decrement by  one
    SEQAN_ASSERT_EQ(registerItem(hotList, 0), 2u);
    SEQAN_ASSERT_EQ(registerItem(hotList, 1), 2u);
}

// Black-box test with larger data amount.
SEQAN_DEFINE_TEST(test_synopsis_frequency_hot_list_frequent_large)
{
    using namespace seqan;

    HotList<int, Frequent> hotList(4);
    testGenericHotListTopK(hotList);
}

// Black-box test with larger data amount.
SEQAN_DEFINE_TEST(test_synopsis_frequency_hot_list_lossy_counting_large)
{
    using namespace seqan;

    HotList<int, LossyCounting> hotList(4);
    testGenericHotListIceberg(hotList);
}

// Black-box test with larger data amount.
SEQAN_DEFINE_TEST(test_synopsis_frequency_hot_list_space_saving_large)
{
    using namespace seqan;

    HotList<int, SpaceSaving> hotList(4);
    testGenericHotListTopK(hotList);
}

#endif  // TEST_SYNOPSIS_TEST_SYNOPSIS_HOT_LIST_H_
