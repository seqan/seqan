// ==========================================================================
// %(TITLE)s
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
// Author: %(AUTHOR)s
// ==========================================================================

#ifndef %(HEADER_GUARD)s
#define %(HEADER_GUARD)s

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>

// A test for strings.
SEQAN_DEFINE_TEST(test_%(NAME)s_strings_example1)
{
    using namespace seqan;

    // Define some constant test data for comparison...
    CharString const STRING1 = "test 1";
    CharString const STRING2 = "test 2";

    // Append to a string and make equality assertion on the result.
    CharString myStr = "test ";
    append(myStr, "1");
    SEQAN_ASSERT_EQ(STRING1, myStr);

    // Demonstration of other assertions.
    SEQAN_ASSERT_GT(STRING2, myStr);
    SEQAN_ASSERT_GEQ(STRING2, myStr);
    SEQAN_ASSERT_LT(myStr, STRING2);
    SEQAN_ASSERT_LEQ(STRING2, STRING2);
}

#endif  // %(HEADER_GUARD)s
