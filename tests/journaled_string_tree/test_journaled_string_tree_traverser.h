// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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

#include "test_config_reader.h"
#include "test_journaled_string_tree.h"

using namespace seqan;

template <typename TJst>
inline TJst _createSimpleJst()
{
    typename Host<typename Source<TJst>::Type>::Type seq = "AGATCGAGCGAGCTAGCGACTCAG";
    TJst jst(seq, 10, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    insertNode(jst, 1, 3, ids, DeltaTypeDel());
    insertNode(jst, 8, "CGTA", ids, DeltaTypeIns());
    insertNode(jst, 10, 'C', ids, DeltaTypeSnp());
    insertNode(jst, 15, 2, ids, DeltaTypeDel());
    insertNode(jst, 20, 'A', ids, DeltaTypeSnp());

    SEQAN_ASSERT(create(jst));
    return jst;
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_constructor)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    {  // Default ctor.
        TTraverser traverser;
        SEQAN_ASSERT(traverser._contPtr == nullptr);
        SEQAN_ASSERT(traverser._contextSize == 1);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(empty(*traverser._stackPtr) == true);
    }

    {  // Constructor with jst
        TJst jst = _createSimpleJst<TJst>();

        TTraverser traverser(jst);
        SEQAN_ASSERT(traverser._contPtr == &jst);
        SEQAN_ASSERT(traverser._contextSize == 1);
        SEQAN_ASSERT(traverser._stackPtr.get() != nullptr);
        SEQAN_ASSERT(length(*traverser._stackPtr) == 1u);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_traverser)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef TraverserImpl<TJst, JstTraversalSpec<> > TTraverser;

    TJst jst = _createSimpleJst<TJst>();


    { // Without observable
        typename Traverser<TJst>::Type subject = traverser(jst);
        SEQAN_ASSERT(subject._contPtr == &jst);
        SEQAN_ASSERT(subject._contextSize == 1);
        SEQAN_ASSERT(length(*subject._stackPtr) == 1);
    }
    { // With observable
        // TODO(rrahn): Write me!
    }

}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_traversal)
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/tests/journaled_string_tree/testConfig.txt");

    TestConfigFileIn_ configIn(toCString(path));

    TestConfigHeader_<DnaString> header;
    readHeader(header, configIn);
    DnaString seq = header.ref;

    JournaledStringTree<DnaString, TestJstPosConfig_> jst(seq, 10, length(context(configIn)));
    _createJst(jst, configIn);

    SEQAN_ASSERT(create(jst));

    StringSet<DnaString> testSeqs;
    resize(testSeqs, dimension(jst));

    std::stringstream strStream;

    auto sub = traverser(jst);

//    unsigned lastId = length(*sub._stackPtr);
//    String<bool, Packed<> > lastCov;
    while (!atEnd(sub))
    {
        //stepSize = presentContext(alg);
//        if (lastId != length(*sub._stackPtr))
//        {
//            std::cout << "Node ID: " << lastId << " Label: " << strStream.str() << " Coverage: " << _printCoverage(lastCov) << std::endl;
//            strStream.str("");
//        }

        unsigned count = 0;
        for (auto it = begin(back(*sub._stackPtr).coverage); it != end(back(*sub._stackPtr).coverage); ++it, ++count)
            if (*it)
                appendValue(testSeqs[count], *(back(*sub._stackPtr).curEdgeIt));

//        std::cout << "Current Sequences: " << std::endl;
//        count = 0;
//        for (auto seq : testSeqs)
//        {
//            std::cout << "seq ";
//            std::cout.fill('0');
//            std::cout.width(2);
//            std::cout << count << ": ";
//            std::cout << seq << '\n';
//            ++count;
//        }

//        lastId = length(*sub._stackPtr);
//        strStream << *(back(*sub._stackPtr).curEdgeIt);
//        lastCov = back(*sub._stackPtr).coverage;
        goNext(sub);

    }

    for (unsigned i = 0; i < length(jst._buffer._journaledSet); ++i)
        SEQAN_ASSERT(testSeqs[i] == jst._buffer._journaledSet[i]);
}


SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_at_end)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_string_context)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_go_next)
{
    SEQAN_ASSERT_FAIL("Write me!");
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_traverser_container)
{
    SEQAN_ASSERT_FAIL("Write me!");
}
#endif // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_TRAVERSER_H_
