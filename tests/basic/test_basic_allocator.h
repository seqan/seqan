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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef TESTS_BASIC_TEST_BASIC_ALLOCATOR_H_
#define TESTS_BASIC_TEST_BASIC_ALLOCATOR_H_

#include <iostream>
#include <memory>
#include <vector>
#include <map>

// TODO(holtgrew): Get rid fo using namespace
using namespace seqan;

// ==========================================================================
// Helper Code
// ==========================================================================

// The following helper class is passed as the parent allocator to all tested
// allocators.  It allows one to count the number of allocations and deallocations
// and checks whether all allocated memory blocks are correctly deallocated
// when it is destructed.

struct TestAllocator
{
    mutable std::map<char *, size_t> data_allocated;
    mutable std::map<char *, size_t> data_deallocated;

    TestAllocator() {}

    ~TestAllocator()
    {
        std::map<char *, size_t>::iterator it = data_allocated.begin();
        while (it != data_allocated.end())
        {
            SEQAN_ASSERT_MSG(data_deallocated.count(it->first), "Memory block not deallocated.");
            deallocate(int(), it->first, it->second);
            ++it;
        }
    }
};

template <typename TValue, typename TSize, typename TUsage>
void allocate(TestAllocator & me,
              TValue * & data_,
              TSize count,
              Tag<TUsage> const)
{
    SEQAN_ASSERT_GT(count, static_cast<TSize>(0));
    allocate(int(), data_, count);
    me.data_allocated[(char *) data_] = count;
}

template <typename TValue, typename TSize, typename TUsage>
void deallocate(TestAllocator & me,
                TValue * data_,
                TSize count,
                Tag<TUsage> const)
{
    SEQAN_ASSERT_MSG(me.data_allocated.count((char *) data_), "memory block was not allocated");
    SEQAN_ASSERT_MSG(me.data_allocated[(char *) data_] == count, "memory block was allocated with different size");
    SEQAN_ASSERT_MSG(!me.data_deallocated.count((char *) data_), "memory block already deallocated");

    me.data_deallocated[(char *) data_] = count;
}

int countAllocs(TestAllocator & me)
{
    return me.data_allocated.size();
}

int countDeallocs(TestAllocator & me)
{
    return me.data_deallocated.size();
}

// ==========================================================================
// Tests
// ==========================================================================

SEQAN_DEFINE_TEST(test_basic_allocator_simple)
{
    int * dat1;
    int * dat2;

    Allocator<SimpleAlloc<TestAllocator> > allo1;
    allocate(allo1, dat1, 100);
    allocate(allo1, dat2, 105);
    deallocate(allo1, dat1, 100);
    allocate(allo1, dat2, 201);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(allo1)), 3);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(allo1)), 1);

    clear(allo1);

    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(allo1)), 3);
}

SEQAN_DEFINE_TEST(test_basic_allocator_pool)
{
    int * dat1;
    int * dat2;

    typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
    Allocator<SinglePool<20 * sizeof(int), TParentAlloc> > allo1;
    allocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);
    deallocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);

    SEQAN_ASSERT_EQ(dat1, dat2);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 1);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 0);

    allocate(allo1, dat1, 100);
    deallocate(allo1, dat1, 100);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 2);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 1);

    clear(allo1);

    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 2);
}

SEQAN_DEFINE_TEST(test_basic_allocator_multi_pool)
{
    int * dat1;
    int * dat2;

    typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
    Allocator<MultiPool<TParentAlloc> > allo1;
    allocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);
    deallocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);

    SEQAN_ASSERT_EQ(dat1, dat2);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 1);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 0);

    allocate(allo1, dat1, 30);
    deallocate(allo1, dat1, 30);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 2);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 0);

    clear(allo1);

    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 2);
}

SEQAN_DEFINE_TEST(test_basic_allocator_chunk_pool)
{
    int * dat1;
    int * dat2;

    typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
    Allocator<ChunkPool<4, 20, TParentAlloc> > allo1;
    allocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);
    deallocate(allo1, dat1, 20);
    allocate(allo1, dat2, 20);

    SEQAN_ASSERT_EQ(dat1, dat2);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 1);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 0);

    allocate(allo1, dat1, 30);
    deallocate(allo1, dat1, 30);

    SEQAN_ASSERT_EQ(countAllocs(parentAllocator(parentAllocator(allo1))), 2);
    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 1);

    clear(allo1);

    SEQAN_ASSERT_EQ(countDeallocs(parentAllocator(parentAllocator(allo1))), 2);
}

#endif  // #ifndef TESTS_BASIC_TEST_BASIC_ALLOCATOR_H_
