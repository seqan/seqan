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
// Tests for fundamental conversion code.
// ==========================================================================

#ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONVERSION_H_
#define SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONVERSION_H_

// Test default implementation of Convert<> metafunction.
SEQAN_DEFINE_TEST(test_basic_fundamental_convert_metafunction)
{
    using namespace seqan;

    SEQAN_ASSERT((+SameType_<typename Convert<int, unsigned>::Type, int>::VALUE));
    SEQAN_ASSERT((+SameType_<typename Convert<double, unsigned>::Type, double>::VALUE));
}

// Test implementation of Convert<> function.
SEQAN_DEFINE_TEST(test_basic_fundamental_convert_function)
{
    using namespace seqan;

    // Compiler would warn if result type incorrect.
    SEQAN_ASSERT_EQ(convert<unsigned>(3), 3u);
    SEQAN_ASSERT_EQ(convert<int>(3u), 3);
}

#endif  // #ifndef SEQAN_TESTS_BASIC_TEST_BASIC_FUNDAMENTAL_CONVERSION_H_
