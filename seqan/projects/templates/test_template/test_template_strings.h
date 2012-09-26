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
// Author: @@Your Name@@ <@@Your Email@@>
// ==========================================================================
// @@By its name, this would test the header template/strings.h.@@
// ==========================================================================

// @@Lines and sections with two at-chars are comments that should either
// be removed or replaced by tests defined by you.@@

// @@Replace symbol accordingly to your file name!@@
#ifndef TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_
#define TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_

#include <seqan/basic.h>     // @@For the testing infrastructure.@@
#include <seqan/sequence.h>  // @@Replace with your header under test.@@


// @@Simple test on string operations and demonstration of assertions.
// Replace this with your own test.@@
SEQAN_DEFINE_TEST(test_template_strings_example1) {
    using namespace seqan;

    // Define some constant test data for comparison...
    CharString const kStr1 = "test 1";
    CharString const kStr2 = "test 2";

    // Append to a string and make equality assertion on the result.
    CharString myStr = "test ";
    append(myStr, "1");
    SEQAN_ASSERT_EQ(kStr1, myStr);

    // Demonstration of other assertions.
    SEQAN_ASSERT_GT(kStr2, myStr);
    SEQAN_ASSERT_GEQ(kStr2, myStr);
    SEQAN_ASSERT_LT(myStr, kStr2);
    SEQAN_ASSERT_LEQ(myStr, kStr2);
}

// @@Replace symbol accordingly to your file name!@@
#endif  // TEST_TEMPLATE_TEST_TEMPLATE_STRINGS_H_
