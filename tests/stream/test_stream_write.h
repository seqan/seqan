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
//         David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Tests for the File Stream.
// ==========================================================================

#ifndef TEST_STREAM_TEST_STREAM_WRITE_H_
#define TEST_STREAM_TEST_STREAM_WRITE_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "test_stream_generic.h"

using namespace seqan;

// ==========================================================================
// Types
// ==========================================================================

template <typename TTarget_>
class WriteTest : public Test
{
public:
    typedef TTarget_ TTarget;
};

// --------------------------------------------------------------------------
// FileStream Specs
// --------------------------------------------------------------------------

typedef
    TagList<CharString
    >
    WriteTargets;

SEQAN_TYPED_TEST_CASE(WriteTest, WriteTargets);

// --------------------------------------------------------------------------
// FileStream Tests
// --------------------------------------------------------------------------

// Simple example of writing a std::string.
SEQAN_TYPED_TEST(WriteTest, StdString)
{
    typename TestFixture::TTarget t1, t2;

    std::string s1 = "This is a string!\nWith two lines.";
    write(t1, s1);
    SEQAN_ASSERT_EQ(t1, s1);

    std::string s2 = "";
    write(t2, s2);
    SEQAN_ASSERT_EQ(t2, s2);

}


#endif  // TEST_STREAM_TEST_STREAM_WRITE_H_
