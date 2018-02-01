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
// Tests for misc simple atomic operations.
// ==========================================================================

#ifndef TEST_PARALLEL_TEST_PARALLEL_ATOMIC_MISC_H_
#define TEST_PARALLEL_TEST_PARALLEL_ATOMIC_MISC_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parallel.h>

template <typename T>
void atomicMinTestImpl(T const &)
{
    using namespace seqan;

    int const ARR_SIZE = 4 * 1024;

    // Generate pseudorandom numbers to compute minimum of.
    T expectedMin = MaxValue<T>::VALUE;
    T arr[ARR_SIZE];
    T volatile x = 3;
    for (int i = 0; i < ARR_SIZE; ++i) {
        x = 32456 * x + 9874;
        arr[i] = x;
        expectedMin = _min(expectedMin, x);
    }

    // Compute minimum in parallel.  Not efficient, but many conflicts
    // should occur.
    x = MaxValue<T>::VALUE;
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < ARR_SIZE; ++i)
        atomicMin(x, arr[i]);

    SEQAN_ASSERT_EQ(expectedMin, x);
}

template <typename T>
void atomicMaxTestImpl(T const &)
{
    using namespace seqan;

    int const ARR_SIZE = 4 * 1024;

    // Generate pseudorandom numbers to compute maximum of.
    T expectedMax = MinValue<T>::VALUE;
    T arr[ARR_SIZE];
    T x = 3;
    for (int i = 0; i < ARR_SIZE; ++i) {
        x = 32456 * x + 9874;
        arr[i] = x;
        expectedMax = _max(expectedMax, x);
    }

    // Compute maximum in parallel.  Not efficient, but many conflicts
    // should occur.
    x = MinValue<T>::VALUE;
    SEQAN_OMP_PRAGMA(parallel for schedule(static, 1))
    for (int i = 0; i < ARR_SIZE; ++i)
        atomicMax(x, arr[i]);

    SEQAN_ASSERT_EQ(expectedMax, x);
}

SEQAN_DEFINE_TEST(test_parallel_atomic_min)
{
    using namespace seqan;
    typedef unsigned short SEQAN_ushort;
    typedef unsigned long SEQAN_ulong;

    // Tests are limited to the types where MSVC allows atomic
    // Compare-And-Swap, also 64 bit CAS is not available on 32 bit Intel.
    atomicMinTestImpl(short());
    atomicMinTestImpl(SEQAN_ushort());
    atomicMinTestImpl(long());
    atomicMinTestImpl(SEQAN_ulong());
#if SEQAN_IS_64_BIT
    atomicMinTestImpl(int64_t());
    atomicMinTestImpl(uint64_t());
#endif  // #if SEQAN_IS_64_BIT
}

SEQAN_DEFINE_TEST(test_parallel_atomic_max)
{
    using namespace seqan;
    typedef unsigned short SEQAN_ushort;
    typedef unsigned long SEQAN_ulong;

    // Tests are limited to the types where MSVC allows atomic
    // Compare-And-Swap, also 64 bit CAS is not available on 32 bit Intel.
    atomicMaxTestImpl(short());
    atomicMaxTestImpl(SEQAN_ushort());
    atomicMaxTestImpl(long());
    atomicMaxTestImpl(SEQAN_ulong());
#if SEQAN_IS_64_BIT
    atomicMaxTestImpl(int64_t());
    atomicMaxTestImpl(uint64_t());
#endif  // #if SEQAN_IS_64_BIT
}

#endif  // TEST_PARALLEL_TEST_PARALLEL_ATOMIC_MISC_H_
