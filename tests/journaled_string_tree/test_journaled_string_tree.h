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

#ifndef EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
#define EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_journaled_string_tree_constructor)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef Size<TJst>::Type TSize;

    {  // 1st Constructor
        TJst jst(100);

        SEQAN_ASSERT_EQ(jst._depth, static_cast<TSize>(100u));
        SEQAN_ASSERT(empty(jst._deltaMap));
        SEQAN_ASSERT(empty(jst._baseSeq));
    }

    {  // 2nd Constructor
        DnaString hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";

        TJst jst(101, hostSeq);
        SEQAN_ASSERT_EQ(jst._depth, static_cast<TSize>(101));
        SEQAN_ASSERT(empty(jst._deltaMap));
        SEQAN_ASSERT_EQ(jst._baseSeq, hostSeq);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_depth)
{
    typedef Size<JournaledStringTree<Dna5String> >::Type TSize;

    JournaledStringTree<Dna5String> jst(10);
    SEQAN_ASSERT_EQ(depth(jst), static_cast<TSize>(10));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_base_seq)
{
    JournaledStringTree<Dna5String> jst(10);
    SEQAN_ASSERT_EQ(length(baseSeq(jst)), 0u);

    JournaledStringTree<Dna5String> jst2(10, "CGTATAGGANNAGAT");
    SEQAN_ASSERT_EQ(baseSeq(jst2), "CGTATAGGANNAGAT");
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_add_node)
{
    JournaledStringTree<Dna5String> jst(10, "CGTATAGGANNAGAT");
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);

    SEQAN_ASSERT(addNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(addNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(addNode(jst, 5, "CGTA", ids, DeltaTypeIns()));
    SEQAN_ASSERT_NOT(addNode(jst, 5, "CGTA", ids, DeltaTypeIns()));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_erase_node)
{
    JournaledStringTree<Dna5String> jst(10, "CGTATAGGANNAGAT");
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);

    SEQAN_ASSERT(addNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(addNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(addNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    SEQAN_ASSERT(eraseNode(jst, 5, DeltaTypeIns()));
    SEQAN_ASSERT(eraseNode(jst, 5, DeltaTypeSnp()));
    SEQAN_ASSERT_NOT(eraseNode(jst, 5, DeltaTypeSnp()));
    SEQAN_ASSERT_NOT(eraseNode(jst, 4, DeltaTypeIns()));
    SEQAN_ASSERT(eraseNode(jst, 2, DeltaTypeDel()));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_clear)
{
    JournaledStringTree<Dna5String> jst(10, "CGTATAGGANNAGAT");
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);

    SEQAN_ASSERT(addNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(addNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(addNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    clear(jst);
    SEQAN_ASSERT(jst._depth == 0);
    SEQAN_ASSERT(empty(jst._baseSeq));
    SEQAN_ASSERT(empty(jst._deltaMap));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_set_base_seq)
{
    JournaledStringTree<Dna5String> jst(10, "CGTATAGGANNAGAT");
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);

    SEQAN_ASSERT(addNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(addNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(addNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    setBaseSeq(jst, "GTTATAGAGGNAGTAGAAAAANN");
    SEQAN_ASSERT(jst._depth == 10u);
    SEQAN_ASSERT_EQ(jst._baseSeq, "GTTATAGAGGNAGTAGAAAAANN");
    SEQAN_ASSERT(empty(jst._deltaMap));
}

#endif // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
