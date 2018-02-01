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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#include "test_index_helpers.h"

using namespace seqan;

// ========================================================================== 
// Test Classes
// ========================================================================== 

SEQAN_TYPED_TEST_CASE(IndexTest, FMIndexTypes);

// ==========================================================================
// Index Tests
// ========================================================================== 

// --------------------------------------------------------------------------
// Test indexCreate()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(IndexTest, IndexCreate)
{
    SEQAN_ASSERT(indexCreate(this->index));

    // Construct without any text should't work.
    clear(this->index);
    SEQAN_ASSERT_NOT(indexCreate(this->index));

// TODO(esiragusa): test using fibre tags.
//    SEQAN_ASSERT(indexCreate(index, FibreSA()));
//    clear(index);
//    SEQAN_ASSERT(indexCreate(index, FibreLF()));
//    clear(index);
//    SEQAN_ASSERT(indexCreate(index, FibreSALF()));

// TODO(esiragusa): use setHost() to set the text.
//    setHost(index, text);
//    SEQAN_ASSERT(indexCreate(index));
}

// --------------------------------------------------------------------------
// Test length()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(IndexTest, Length)
{
    SEQAN_ASSERT_EQ(length(this->index), lengthSum(this->text));
    clear(this->index);
    SEQAN_ASSERT_EQ(length(this->index), 0u);
}

// --------------------------------------------------------------------------
// Test clear() and empty()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(IndexTest, ClearEmpty)
{
    SEQAN_ASSERT_NOT(empty(this->index));
    indexCreate(this->index);
    SEQAN_ASSERT_NOT(empty(this->index));
    clear(this->index);
    SEQAN_ASSERT(empty(this->index));
}

// --------------------------------------------------------------------------
// Test open() and save()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(IndexTest, OpenSave)
{
    indexCreate(this->index);

    const char * fileName = SEQAN_TEMP_FILENAME();
    SEQAN_ASSERT(save(this->index, fileName));
    clear(this->index);

    SEQAN_ASSERT(open(this->index, fileName));
    SEQAN_ASSERT_NOT(empty(this->index));
    SEQAN_ASSERT_EQ(length(this->index), lengthSum(this->text));
}

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
