// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for parallel queue.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_
#define TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

SEQAN_DEFINE_TEST(test_parallel_queue_simple)
{
    seqan::ConcurrentQueue<int> queue;
    int x = -1;
    SEQAN_ASSERT_NOT(tryPopFront(x, queue));

    appendValue(queue, 3);
    appendValue(queue, 6);

    SEQAN_ASSERT(tryPopFront(x, queue));
    SEQAN_ASSERT_EQ(x, 3);
    SEQAN_ASSERT(tryPopFront(x, queue));
    SEQAN_ASSERT_EQ(x, 6);
    SEQAN_ASSERT_NOT(tryPopFront(x, queue));

    for (unsigned i = 0; i < 10; ++i)
    {
        seqan::ConcurrentQueue<int> queue(i);
        int x = -1;
        SEQAN_ASSERT_NOT(tryPopFront(x, queue));
        SEQAN_ASSERT(empty(queue));

        appendValue(queue, 3);
        appendValue(queue, 6);

        SEQAN_ASSERT(tryPopFront(x, queue));
        SEQAN_ASSERT_EQ(x, 3);
        SEQAN_ASSERT(tryPopFront(x, queue));
        SEQAN_ASSERT_EQ(x, 6);
        SEQAN_ASSERT_NOT(tryPopFront(x, queue));
    }
}

SEQAN_DEFINE_TEST(test_parallel_queue_resize)
{
    // in a queue of capacity 3 try all 3 states of being empty
    for (int ofs = 0; ofs < 3; ++ofs)
    {
        seqan::ConcurrentQueue<int> queue(2);

        for (int i = 0; i < ofs; ++i)
        {
            int x = -1;
            appendValue(queue, i);
            SEQAN_ASSERT(tryPopFront(x, queue));
            SEQAN_ASSERT_EQ(x, i);
        }

        SEQAN_ASSERT(empty(queue));

        // fill queue and enforce capacity upgrade
        for (int i = 0; i < 2; ++i)
        {
            appendValue(queue, 10 + i);
            SEQAN_ASSERT_EQ(length(queue), (unsigned)(i + 1));
        }
        SEQAN_ASSERT_EQ(capacity(queue), 3u);
        appendValue(queue, 12);
        SEQAN_ASSERT_GT(capacity(queue), 3u);
        SEQAN_ASSERT_EQ(length(queue), 3u);

        for (int i = 0; i < 3; ++i)
        {
            int x = -1;
            SEQAN_ASSERT(tryPopFront(x, queue));
            SEQAN_ASSERT_EQ(x, 10 + i);
            SEQAN_ASSERT_EQ(length(queue), (unsigned)(2 - i));
        }

        int x = -1;
        SEQAN_ASSERT_NOT(tryPopFront(x, queue));
        SEQAN_ASSERT(empty(queue));
    }
}

SEQAN_DEFINE_TEST(test_parallel_queue_parallel_access)
{
    typedef seqan::ConcurrentQueue<unsigned> TQueue;

    TQueue queue;
    seqan::String<unsigned> random;
    seqan::Rng<seqan::MersenneTwister> rng(0);

    unsigned chkSum = 0;

    resize(random, 1000000);
    for (unsigned i = 0; i < length(random); ++i)
    {
//        random[i] = pickRandomNumber(rng);
        random[i] = i;
        chkSum ^= random[i];
    }

    volatile unsigned chkSum2 = 0;
    size_t threadCount = omp_get_max_threads();
    seqan::Splitter<unsigned> splitter(0, length(random), (threadCount + 1) / 2);
    SEQAN_OMP_PRAGMA(parallel)
    {
        if ((omp_get_thread_num() & 1) == 0)
        {
            seqan::ScopedWriteLock<TQueue> writeLock(queue);
//            SEQAN_OMP_PRAGMA(critical(cout))
//            {
//                std::cout << "start writer thread: " << omp_get_thread_num() << std::endl;
//            }
            unsigned tid = omp_get_thread_num() / 2;
            for (unsigned i = splitter[tid]; i != splitter[tid + 1]; ++i)
                appendValue(queue, random[i]);
        }

        if (threadCount < 2 || (omp_get_thread_num() & 1) == 1)
        {
            seqan::ScopedReadLock<TQueue> readLock(queue);
            // wait for queue to be filled initialy
            while (capacity(queue) == 0)
            {}

//            SEQAN_OMP_PRAGMA(critical(cout))
//            {
//                std::cout << "start reader thread: " << omp_get_thread_num() << std::endl;
//            }

            unsigned chkSumLocal = 0, val = 0;
            while (popFront(val, queue))
                chkSumLocal ^= val;
            seqan::atomicXor(chkSum2, chkSumLocal);
        }
    }
//    std::cout << "len: " << length(queue) << std::endl;
//    std::cout << "cap: " << capacity(queue) << std::endl;
    SEQAN_ASSERT_EQ(chkSum, chkSum2);
}


#endif  // TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_
