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
// Tests for the journaled string tree.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSER_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSER_H_

#include <sstream>

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>
#include <seqan/stream.h>

#include "test_journaled_string_tree_mock.h"

using namespace seqan;

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_constructor)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    {  // Default ctor.
        TTraverser traverser;
        SEQAN_ASSERT(traverser._contPtr == nullptr);
        SEQAN_ASSERT(traverser._branchLength == 1);
        SEQAN_ASSERT(traverser._contextSize == 1);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    {  // Constructor with jst
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser traverser(jst);
        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._branchLength == 1u);
        SEQAN_ASSERT(traverser._contextSize == 1u);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    {  // Constructor with jst, contextSize
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser traverser(jst, 2);
        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._branchLength == 2u);
        SEQAN_ASSERT(traverser._contextSize == 2u);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    {  // Constructor with jst, contextSize, branchSize
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser traverser(jst, 2, 5);
        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._branchLength == 5u);
        SEQAN_ASSERT(traverser._contextSize == 2u);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    {  // Copy C'tor
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser traverser(jst, 2, 5);
        TTraverser travCopy(traverser);

        SEQAN_ASSERT(travCopy._contPtr == &jst);
        SEQAN_ASSERT(travCopy._branchLength == 5u);
        SEQAN_ASSERT(travCopy._contextSize == 2u);
        SEQAN_ASSERT(travCopy._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*travCopy._stackPtr) == true);

        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._branchLength == 5u);
        SEQAN_ASSERT(traverser._contextSize == 2u);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    { // Move C'tor
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser travMove(TTraverser(jst, 2, 5));

        SEQAN_ASSERT(travMove._contPtr == &jst);
        SEQAN_ASSERT(travMove._branchLength == 5u);
        SEQAN_ASSERT(travMove._contextSize == 2u);
        SEQAN_ASSERT(travMove._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*travMove._stackPtr) == true);
    }

    { // Copy Assignment
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser traverser(jst, 2, 5);
        TTraverser travCopy;
        travCopy = traverser;

        SEQAN_ASSERT(travCopy._contPtr == &jst);
        SEQAN_ASSERT(travCopy._branchLength == 5u);
        SEQAN_ASSERT(travCopy._contextSize == 2u);
        SEQAN_ASSERT(travCopy._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*travCopy._stackPtr) == true);

        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._branchLength == 5u);
        SEQAN_ASSERT(traverser._contextSize == 2u);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    { // Move Assignment
        TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

        TTraverser travMove;
        travMove = TTraverser(jst, 2, 5);

        SEQAN_ASSERT(travMove._contPtr == &jst);
        SEQAN_ASSERT(travMove._branchLength == 5u);
        SEQAN_ASSERT(travMove._contextSize == 2u);
        SEQAN_ASSERT(travMove._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*travMove._stackPtr) == true);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_traverser)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef Traverser<TJst>::Type TTraverser;
    // Test Traverser Metafunction!

    bool res = IsSameType<TTraverser, TraverserImpl<TJst, JstTraversalSpec<void> > >::VALUE;
    SEQAN_ASSERT(res);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_init)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();
    auto observer = makeObserverList();
    {
        TTraverser test(jst, 1);
        init(test, observer);

        SEQAN_ASSERT(test._contPtr == &jst);
        SEQAN_ASSERT(test._branchLength == 1u);
        SEQAN_ASSERT(test._contextSize == 1u);
        SEQAN_ASSERT(length(*test._stackPtr) == 1u);
        SEQAN_ASSERT_EQ(*(back(*test._stackPtr).curEdgeIt), 'A');
    }

    {
        TTraverser test(jst, 5);
        init(test, observer);

        SEQAN_ASSERT(test._contPtr == &jst);
        SEQAN_ASSERT(test._branchLength == 5u);
        SEQAN_ASSERT(test._contextSize == 5u);
        SEQAN_ASSERT(length(*test._stackPtr) == 2u);
        SEQAN_ASSERT(*(back(*test._stackPtr).curEdgeIt) == 'G');
    }

    {
        TTraverser test(jst, 10);
        init(test, observer);

        SEQAN_ASSERT(test._contPtr == &jst);
        SEQAN_ASSERT(test._branchLength == 10u);
        SEQAN_ASSERT(test._contextSize == 10u);
        SEQAN_ASSERT_EQ(length(*test._stackPtr), 4u);
        SEQAN_ASSERT(*(back(*test._stackPtr).curEdgeIt) == 'G');
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_set_container)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

    TTraverser test;
    SEQAN_ASSERT_EQ(test._contPtr, nullptr);

    setContainer(test, jst);
    SEQAN_ASSERT_EQ(test._contPtr, &jst);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_container)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

    TTraverser test(jst);
    SEQAN_ASSERT_EQ(&container(test), &jst);

    TTraverser const testC = test;
    SEQAN_ASSERT_EQ(&container(test), &jst);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_context_iterator)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();
    TTraverser test(jst);
    auto observer = makeObserverList();
    init(test, observer);

    SEQAN_ASSERT(contextIterator(test) == begin(impl::member(jst, JstSourceMember()), Standard()));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_advance)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    auto observer = makeObserverList();

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();
    TTraverser test(jst, 1, 10);
    init(test, observer);
    SEQAN_ASSERT(*contextIterator(test) == 'A');
    advance(test, 1);
    SEQAN_ASSERT(*contextIterator(test) == 'C');
    advance(test, 6);
    SEQAN_ASSERT(*contextIterator(test) == 'T');
    advance(test, 2);
    SEQAN_ASSERT(*contextIterator(test) == 'C');
    advance(test, 7);
    SEQAN_ASSERT(*contextIterator(test) == 'A');
    advance(test, 10);
    SEQAN_ASSERT(*contextIterator(test) == 'G');
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_at_end)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();
    TTraverser test(jst, 1);

    SEQAN_ASSERT_EQ(atEnd(test), false);
    for (unsigned i = 0; i < 31; ++i)
    {
        advance(test, 1);
    }
    SEQAN_ASSERT_EQ(atEnd(test), true);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_is_base)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();
    TTraverser test(jst, 1);

    SEQAN_ASSERT(!isBase(test));
    for (unsigned i = 0; i < 31; ++i)
    {
        advance(test, 1);
    }
    SEQAN_ASSERT(isBase(test));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_context_size)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

    TTraverser test(jst);

    SEQAN_ASSERT_EQ(contextSize(test), 1u);

    setContextSize(test, 5);
    SEQAN_ASSERT_EQ(contextSize(test), 5u);
    SEQAN_ASSERT_EQ(test._needInitialization, true);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_branch_size)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = JstMockGenerator::_createSimpleJst<TJst>();

    TTraverser test(jst);

    SEQAN_ASSERT_EQ(branchSize(test), 1u);

    setBranchSize(test, 5);
    SEQAN_ASSERT_EQ(branchSize(test), 5u);
    SEQAN_ASSERT_EQ(test._needInitialization, true);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_basic_traversal)
{
    typedef JournaledStringTree<DnaString> TJst;
    auto jst = JstMockGenerator::_createComplexJst<TJst>();

    StringSet<DnaString> testSeqs;
    resize(testSeqs, length(jst));

    std::stringstream strStream;

    typename Traverser<TJst>::Type sub(jst, 1, 1);
    auto observer = makeObserverList();
    init(sub, observer, SelectFirstProxy());

    while (!atEnd(sub))
    {
        auto const & pos = position(sub);
        for (auto it = begin(pos), itEnd = end(pos); it != itEnd; ++it)
        {
            auto const & p = *it;
            appendValue(testSeqs[p.i1], *(impl::activeNode(sub).curEdgeIt));
        }

        advance(sub, 1, SelectFirstProxy());

    }

    for (unsigned i = 0; i < length(impl::buffer(sub)._journaledSet); ++i)
    {
        SEQAN_ASSERT_MSG(testSeqs[i] == impl::buffer(sub)._journaledSet[i], "Sequence %d did not match", i);
    }
}

#endif // TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSER_H_
