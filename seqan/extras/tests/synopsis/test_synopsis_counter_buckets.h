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
// Tests for the CounterBuckets data structure.
// ==========================================================================

#ifndef TEST_SYNOPSIS_TEST_SYNOPSIS_COUNTER_BUCKETS_H_
#define TEST_SYNOPSIS_TEST_SYNOPSIS_COUNTER_BUCKETS_H_

#include <list>

SEQAN_DEFINE_TEST(test_synopsis_counter_buckets_simple)
{
    using namespace seqan;

    typedef CounterBuckets<unsigned, unsigned> TCounterBuckets;
    typedef Iterator<TCounterBuckets, Buckets>::Type TBucketsIterator;
    typedef Iterator<TCounterBuckets, Entries>::Type TEntriesIterator;

    TCounterBuckets cb;
    TEntriesIterator it1 = addCounter(cb, 1);
    TEntriesIterator it2 = addCounter(cb, 2);
    TEntriesIterator it3 = addCounter(cb, 3);

    SEQAN_ASSERT_EQ(value(*it1), 1u);
    SEQAN_ASSERT_EQ(value(*it2), 1u);
    SEQAN_ASSERT_EQ(value(*it3), 1u);

    increaseCounter(cb, it2);
    increaseCounter(cb, it2);
    increaseCounter(cb, it2);
    increaseCounter(cb, it3);

    SEQAN_ASSERT_EQ(value(*it1), 1u);
    SEQAN_ASSERT_EQ(value(*it2), 4u);
    SEQAN_ASSERT_EQ(value(*it3), 2u);

    decreaseAllCounters(cb);

    SEQAN_ASSERT_EQ(value(*it1), 0u);
    SEQAN_ASSERT_EQ(value(*it2), 3u);
    SEQAN_ASSERT_EQ(value(*it3), 1u);
}

#endif  // TEST_SYNOPSIS_TEST_SYNOPSIS_COUNTER_BUCKETS_H_
