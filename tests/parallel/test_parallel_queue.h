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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for parallel queue.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_
#define TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_

#include <random>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

#include <chrono>

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
        SEQAN_ASSERT_EQ(capacity(queue), 2u);
        appendValue(queue, 12);
        SEQAN_ASSERT_GT(capacity(queue), 2u);
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

SEQAN_DEFINE_TEST(test_parallel_queue_non_pod)
{
//    typedef std::string TValue;
    typedef seqan::CharString TValue;

    // in a queue of capacity 3 try all 3 states of being empty
    for (int ofs = 1; ofs < 10; ++ofs)
    {
        seqan::ConcurrentQueue<TValue> queue(10);

        for (int i = 0; i < ofs; ++i)
        {
            TValue x;
            TValue y = "al_";
            y[2] = '0' + i;
            appendValue(queue, TValue(y));
            SEQAN_ASSERT(tryPopFront(x, queue));
            SEQAN_ASSERT_EQ(x, y);
        }
        SEQAN_ASSERT(empty(queue));

        for (int i = 0; i < 4; ++i)
        {
            TValue y = "al_";
            y[2] = '0' + i;
            appendValue(queue, TValue(y));
        }
    }
}

template <typename TResizeTag, typename TParallelPop, typename TParallelPush>
void testMPMCQueue(size_t initialCapacity)
{
    typedef seqan::ConcurrentQueue<unsigned> TQueue;

    TQueue queue(initialCapacity);
    seqan::String<unsigned> random;
    std::mt19937 rng(0);

    unsigned chkSum = 0;

    resize(random, 100000);
    for (unsigned i = 0; i < length(random); ++i)
    {
        random[i] = rng();
//        random[i] = i;
        chkSum ^= random[i];
    }
//    std::cout <<chkSum<<std::endl;

    volatile unsigned chkSum2 = 0;
    size_t threadCount = std::thread::hardware_concurrency();

    // limit thread count as virtualbox (used by Travis) seems to have problems with thread congestion
    if (threadCount > 4)
        threadCount = 4;

    size_t writerCount = threadCount / 2;
    if (seqan::IsSameType<TParallelPush, seqan::Serial>::VALUE)
        writerCount = 1;

    if (seqan::IsSameType<TParallelPop, seqan::Serial>::VALUE)
        threadCount = writerCount + 1;

    std::cout << "threads: " << threadCount << std::endl;
    std::cout << "writers: " << writerCount << std::endl;

    SEQAN_ASSERT_GEQ(threadCount, 2u);
    seqan::Splitter<unsigned> splitter(0, length(random), writerCount);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    std::vector<std::thread> workers;
    for (size_t tid = 0; tid < threadCount; ++tid)
    {
        workers.push_back(std::thread([&,tid]()
        {
            if (tid < writerCount)
            {
                seqan::ScopedWriteLock<TQueue> writeLock(queue);
                // barrier for all writers to set up
                waitForWriters(queue, writerCount);

    //            printf("start writer #%ld\n", tid);
                for (unsigned j = splitter[tid]; j != splitter[tid + 1]; ++j)
                {
                    appendValue(queue, random[j], TResizeTag(), TParallelPush());
                }
    //            printf("stop writer #%ld %d\n", tid, splitter[tid + 1] - splitter[tid]);
            }

            if (tid >= writerCount)
            {
                seqan::ScopedReadLock<TQueue> readLock(queue);
                // barrier for all writers to set up
                waitForFirstValue(queue);

    //            printf("start reader #%lu\n",  (long unsigned)tid);
                unsigned chkSumLocal = 0, val = 0, cnt = 0;
                while (popFront(val, queue, TParallelPop()))
                {
                    chkSumLocal ^= val;
                    ++cnt;
    //                if ((cnt & 0xff) == 0)
    //                    printf("%ld ", tid);
                }
                seqan::atomicXor(chkSum2, chkSumLocal);
                printf("stop reader #%lu %d\n", (long unsigned)tid, cnt);
            }
        }));
    }

    for (auto &t : workers)
        t.join();

    std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
    double timeSpan = std::chrono::duration_cast<std::chrono::duration<double> >(stop - start).count();
    std::cout << "throughput: " << (uint64_t)(length(random) / timeSpan) << " values/s" << std::endl;


//    std::cout << "len: " << length(queue) << std::endl;
    std::cout << "cap: " << capacity(queue) << std::endl;
//    std::cout << "pushed: " << queue.pushed << std::endl;
//    std::cout << "popped: " << queue.popped << std::endl;
    SEQAN_ASSERT_EQ(chkSum, chkSum2);
}

SEQAN_DEFINE_TEST(test_parallel_queue_spsc_dynamicsize)
{
    testMPMCQueue<seqan::Generous, seqan::Serial, seqan::Serial>(0u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_spsc_fixedsize)
{
    testMPMCQueue<seqan::Limit, seqan::Serial, seqan::Serial>(30u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_spmc_dynamicsize)
{
    testMPMCQueue<seqan::Generous, seqan::Parallel, seqan::Serial>(0u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_spmc_fixedsize)
{
    testMPMCQueue<seqan::Limit, seqan::Parallel, seqan::Serial>(30u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_mpsc_dynamicsize)
{
    testMPMCQueue<seqan::Generous, seqan::Serial, seqan::Parallel>(0u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_mpsc_fixedsize)
{
    testMPMCQueue<seqan::Limit, seqan::Serial, seqan::Parallel>(30u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_mpmc_dynamicsize)
{
    testMPMCQueue<seqan::Generous, seqan::Parallel, seqan::Parallel>(0u);
}

SEQAN_DEFINE_TEST(test_parallel_queue_mpmc_fixedsize)
{
    testMPMCQueue<seqan::Limit, seqan::Parallel, seqan::Parallel>(30u);
}

#endif  // TEST_PARALLEL_TEST_PARALLEL_QUEUE_H_
