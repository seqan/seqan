// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Test the interface of the Seed class for specializations Simple Seed and
// Chained Seed.
// ==========================================================================

#ifndef TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_
#define TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_

#include <seqan/basic.h>  // Includes testing infrastructure.
#include <seqan/file.h>   // Required to print strings in tests.

#include <seqan/seeds.h>  // Include module under test.

#include <set>

// Test constructors of the indirect iterator.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_constructors)
{
    using namespace seqan;

    SEQAN_ASSERT_FAIL("Write me!");

    // Default constructor.
    {
    }
    // Construct from wrapped iterator.
    {
    }
    // Copy constructor.
    {
    }
}

// Test the metafunctions.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_metafunctions)
{
    using namespace seqan;

    SEQAN_ASSERT_FAIL("Write me!");

    // Iterator
    {
    }

    // Const Iterator
    {
    }
}

// Test the common iterator functions.
SEQAN_DEFINE_TEST(test_seeds_basic_iter_indirect_basic_functions)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

#endif  // #ifndef TEST_SEEDS_TEST_BASIC_ITER_INDIRECT_H_
