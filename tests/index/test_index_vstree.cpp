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

#include <seqan/basic.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// ==========================================================================
// Types
// ========================================================================== 

// ========================================================================== 
// Test Classes
// ========================================================================== 

// --------------------------------------------------------------------------
// Class IndexIteratorTest
// --------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
class IndexIteratorTest : public IndexTest<TIndex>
{
public:
    typedef IndexTest<TIndex>                       TBase;
    typedef typename Iterator<TIndex, TSpec>::Type  TIter;

    // Disable base class index creation since iterators uses lazy construction.
    void setUp()
    {
        createText(this->text, typename TBase::TValue());
    }
};

// --------------------------------------------------------------------------
// Class TopDownIndexIteratorTest
// --------------------------------------------------------------------------

template <typename TIndex>
class TopDownIndexIteratorTest : public IndexIteratorTest<TIndex, TopDown<> > {};

SEQAN_TYPED_TEST_CASE(TopDownIndexIteratorTest, TrieIndexTypes);

// ==========================================================================
// TopDown Iterator Tests
// ==========================================================================

// --------------------------------------------------------------------------
// Test Constructor
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(TopDownIndexIteratorTest, Constructor)
{
    typename TestFixture::TIter it(this->index);

    SEQAN_ASSERT(isRoot(it));
    // NOTE(esiragusa): An empty index must yield isLeaf(it) == true
    SEQAN_ASSERT_NOT(isLeaf(it));

    // NOTE(esiragusa): Actual behavior is:
    SEQAN_ASSERT_EQ(countOccurrences(it), lengthSum(this->text) + countSequences(this->index));
    // NOTE(esiragusa): Correct behavior should be:
//    SEQAN_ASSERT_EQ(countOccurrences(it), lengthSum(this->text));
}

// --------------------------------------------------------------------------
// Test goDown()
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(TopDownIndexIteratorTest, GoDown)
{
    typename TestFixture::TIter it(this->index);

    SEQAN_ASSERT(goDown(it));
    SEQAN_ASSERT_EQ(repLength(it), 1u);

    goRoot(it);

    SEQAN_ASSERT(goDown(it, front(concat(this->text))));
    SEQAN_ASSERT_EQ(parentEdgeLabel(it), front(concat(this->text)));
    SEQAN_ASSERT_EQ(repLength(it), 1u);

    goRoot(it);

    // NOTE(esiragusa): Actual behavior is:
//    SEQAN_ASSERT(goDown(it, prefix(concat(this->text), 5u)));
//    SEQAN_ASSERT_EQ(repLength(it), 5u);
}

// --------------------------------------------------------------------------
// Test goRight()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test isRoot()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test isLeaf()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test range()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test countOccurrences()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test getOccurrences()
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// Test representative()
// --------------------------------------------------------------------------

// ==========================================================================
// ParentLinks TopDown Iterator Tests
// ==========================================================================

// --------------------------------------------------------------------------
// Test goUp()
// --------------------------------------------------------------------------

// ========================================================================== 
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
