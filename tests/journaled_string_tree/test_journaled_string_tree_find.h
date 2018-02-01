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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Tests for the journaled string tree find module.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_FIND_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_FIND_H_

#include <sstream>
#include <vector>

#include <seqan/find.h>

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>
#include <seqan/stream.h>

#include "test_journaled_string_tree_mock.h"

using namespace seqan;

template <typename TTraverser>
struct TestHelperHitCollector_
{
    TTraverser* travPtr;
    StringSet<String<unsigned> > positions;

    TestHelperHitCollector_(TTraverser & trav) : travPtr(&trav)
    {
        resize(positions, length(JstMockGenerator::_createComplexTestSet()), Exact());
    }

    inline void
    operator()()
    {
        auto pos = position(*travPtr);
        for (auto p : pos)
        {
            appendValue(positions[p.i1], p.i2);
        }
    }
};

template <typename TNeedle, typename TSpec, typename TSize>
inline void _testFindJst(Pattern<TNeedle, TSpec> & p,
                         TSize contextLength,
                         TSize windowLength)
{
    typedef JournaledStringTree<DnaString> TJst;

    auto jst = JstMockGenerator::_createComplexJst<TJst>();

    typename Traverser<TJst>::Type trav(jst, contextLength, windowLength);

    JstExtension<Pattern<TNeedle, TSpec> > ext(p);

    TestHelperHitCollector_<decltype(trav)> delegate(trav);

#if defined(JST_FIND_DEBUG)
    clear(__testSet);
    resize(__testSet, length(jst));
#endif

    find(trav, ext, delegate);

    auto _compareTestSetM = JstMockGenerator::_createComplexTestSet();
    for (unsigned pos = 0; pos < length(_compareTestSetM); ++pos)
    {
        Finder<DnaString> finder(_compareTestSetM[pos]);
        unsigned counter = 0;
        while (find(finder, p))
        {
            if (IsSameType<TSpec, Horspool>::VALUE)
                SEQAN_ASSERT_EQ(delegate.positions[pos][counter], position(finder) + length(needle(p)) - 1);
            else
                SEQAN_ASSERT_EQ(delegate.positions[pos][counter], finder.data_endPos - 1);
            ++counter;
        }
        SEQAN_ASSERT(counter == length(delegate.positions[pos]));
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_find_horspool)
{
    for (auto ndl : JstMockGenerator::generateNeedles())
    {
        Pattern<DnaString, Horspool> p(ndl);
        _testFindJst(p, length(needle(p)), length(needle(p)));
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_find_shiftand)
{
    for (auto ndl : JstMockGenerator::generateNeedles())
    {
        Pattern<DnaString, ShiftAnd> p(ndl);
        _testFindJst(p, length(needle(p)), length(needle(p)));
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_find_shiftor)
{
    for (auto ndl : JstMockGenerator::generateNeedles())
    {
        Pattern<DnaString, ShiftOr> p(ndl);
        _testFindJst(p, length(needle(p)), length(needle(p)));
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_find_myers)
{
    for (auto ndl : JstMockGenerator::generateNeedles())
    {
        Pattern<DnaString, MyersUkkonen> p(ndl, -2);
        _testFindJst(p, length(needle(p)) + 2, length(needle(p)) + 2);
    }

#if defined(JST_FIND_DEBUG)
    _printTestSet();

    for(unsigned i = 0; i < length(__testSet); ++i)
        if (__testSet[i] != JstMockGenerator::_createComplexTestSet()[i])
            std::cout << i << std::endl;
#endif
}

#endif // TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_FIND_H_
