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

#include <vector>
#include <seqan/basic.h>
#include <seqan/index.h>

#include "test_index_helpers.h"

using namespace seqan;

// ========================================================================== 
// Test Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class IndexFinderTest
// --------------------------------------------------------------------------

template <typename TIndex_>
class IndexFinderTest : public IndexTest<TIndex_>
{
public:
    typedef TIndex_             TIndex;
    typedef IndexTest<TIndex>   TBase;
    typedef Finder<TIndex>      TFinder;

    TFinder finder;

    IndexFinderTest() :
        finder(TBase::index)
    {}
};

SEQAN_TYPED_TEST_CASE(IndexFinderTest, UnidirectionalIndexTypes);

// --------------------------------------------------------------------------
// Test Finder
// --------------------------------------------------------------------------

SEQAN_TYPED_TEST(IndexFinderTest, FindFirstChar)
{
    SEQAN_ASSERT(find(this->finder, prefix(concat(this->text), 1u)));
    SEQAN_ASSERT(find(this->finder, "A"));
}

SEQAN_TYPED_TEST(IndexFinderTest, DefaultFinder)
{
    typedef typename TestFixture::TIndex TIndex;
    Finder<TIndex> finder;
    find(finder, "needle");
}

SEQAN_TYPED_TEST(IndexFinderTest, StdString)
{
    typedef Index<std::string, IndexSa<> > TIndex;
    std::string str;
    TIndex index(str);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

SEQAN_TYPED_TEST(IndexFinderTest, StdVector)
{
    typedef Index<std::vector<char>, IndexSa<> > TIndex;
    std::vector<char> vec;
    TIndex index(vec);
    Iterator<TIndex, TopDown<> >::Type iter(index);
}

// ==========================================================================
// Functions
// ========================================================================== 

int main(int argc, char const ** argv)
{
    TestSystem::init(argc, argv);
    return TestSystem::runAll();
}
