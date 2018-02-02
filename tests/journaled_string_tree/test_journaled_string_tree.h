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

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_journaled_string_tree_constructor)
{
    typedef JournaledStringTree<DnaString>          TJst;
    typedef JournaledStringTree<DnaString const>    TJst2;
    typedef Size<TJst>::Type                        TSize;

    DnaString hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    DnaString const hostSeqConst = hostSeq;

    {  // Custom c'tor with only the set size.
        TJst jst(10);

        SEQAN_ASSERT(jst._dimension == 10u);
        SEQAN_ASSERT(empty(jst._map));
        SEQAN_ASSERT(empty(jst._source));
    }

    {  // Custom c'tor with base seq and set size.

        TJst jst1(hostSeq, 100);

        SEQAN_ASSERT_EQ(jst1._dimension, static_cast<TSize>(100));
        SEQAN_ASSERT(empty(jst1._map));
        SEQAN_ASSERT(!empty(jst1._source));
        SEQAN_ASSERT(jst1._source == hostSeq);
        SEQAN_ASSERT(&host(jst1._source) == &hostSeq);

        TJst2 jst2(hostSeqConst, 100);
        SEQAN_ASSERT(jst2._source == hostSeqConst);

        TJst jst3("ACGTAGAACTTGA", 100);
        SEQAN_ASSERT_EQ(jst3._dimension, static_cast<TSize>(100));
        SEQAN_ASSERT(empty(jst3._map));
        SEQAN_ASSERT(!empty(jst3._source));
        SEQAN_ASSERT(jst3._source == "ACGTAGAACTTGA");
    }

    {  // Copy c'tor
    }

    {  // Move c'tor
    }
}

// TODO(rrahn): Add test for assignment.

SEQAN_DEFINE_TEST(test_journaled_string_tree_length)
{
    typedef Size<JournaledStringTree<Dna5String> >::Type TSize;
    JournaledStringTree<Dna5String> jst(100);
    SEQAN_ASSERT_EQ(length(jst), static_cast<TSize>(100));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_host)
{
    Dna5String seq = "CGTATAGGANNAGAT";

    {
        JournaledStringTree<Dna5String> jst(100);
        SEQAN_ASSERT(empty(host(jst)));
    }

    {
        JournaledStringTree<Dna5String> jst(seq, 100);
        SEQAN_ASSERT(&host(jst) == &seq);
        JournaledStringTree<Dna5String> const jstConst = jst;
        SEQAN_ASSERT(&host(jstConst) == &seq);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_set_host)
{
    Dna5String seq = "CGTATAGGANNAGAT";

    JournaledStringTree<Dna5String> jst(100);
    SEQAN_ASSERT(empty(host(jst)));
    setHost(jst, seq);
    SEQAN_ASSERT(&host(jst) == &seq);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_insert)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    String<bool, Packed<> > cov;
    resize(cov, 100, false);
    cov[0] = true;
    cov[3] = true;
    cov[9] = true;
    cov[99] = true;

    insert(jst, 5, 'C', ids, DeltaTypeSnp());
    insert(jst, 2, 3, ids, DeltaTypeDel());
    insert(jst, 5, "CGTA", ids, DeltaTypeIns());

    auto it = begin(impl::member(jst, JstDeltaMapMember()), Standard());
    typedef typename RemoveReference<decltype(*it)>::Type TEntry;
    typedef typename RemoveReference<decltype(getDeltaRecord(*it))>::Type TRecord;

    SEQAN_ASSERT_EQ(*it, TEntry(2, TRecord(DELTA_TYPE_DEL, 0), cov, DeltaEndType::IS_LEFT));
    ++it;
    SEQAN_ASSERT_EQ(*it, TEntry(4, TRecord(DELTA_TYPE_DEL, 0), cov, DeltaEndType::IS_RIGHT));
    ++it;
    SEQAN_ASSERT_EQ(*it, TEntry(5, TRecord(DELTA_TYPE_SNP, 0), cov, DeltaEndType::IS_BOTH));
    ++it;
    SEQAN_ASSERT_EQ(*it, TEntry(5, TRecord(DELTA_TYPE_INS, 0), cov, DeltaEndType::IS_BOTH));

}

SEQAN_DEFINE_TEST(test_journaled_string_tree_erase)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    insert(jst, 5, 'C', ids, DeltaTypeSnp());
    insert(jst, 2, 3, ids, DeltaTypeDel());
    insert(jst, 5, "CGTA", ids, DeltaTypeIns());

    SEQAN_ASSERT_EQ(erase(jst, 5, DeltaTypeIns()), 1u);
    SEQAN_ASSERT_EQ(erase(jst, 5, DeltaTypeSnp()), 1u);
    SEQAN_ASSERT_EQ(erase(jst, 5, DeltaTypeSnp()), 0u);
    SEQAN_ASSERT_EQ(erase(jst, 4, DeltaTypeIns()), 0u);
    SEQAN_ASSERT_EQ(erase(jst, 2, DeltaTypeDel()), 2u);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_clear)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 99);

    insert(jst, 5, 'C', ids, DeltaTypeSnp());
    insert(jst, 2, 3, ids, DeltaTypeDel());
    insert(jst, 5, "CGTA", ids, DeltaTypeIns());

    clear(jst);
    SEQAN_ASSERT(jst._dimension == 0);
    SEQAN_ASSERT(empty(jst._source));
    SEQAN_ASSERT(empty(jst._map));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_resize)
{
    Dna5String seq = "CGTATAGGANNAGAT";

    {  // Resize on empty jst
        JournaledStringTree<Dna5String> jst(seq, 100);
        resize(jst, 50);
        SEQAN_ASSERT_EQ(length(jst), 50u);

        resize(jst, 150);
        SEQAN_ASSERT_EQ(length(jst), 150u);
    }

    {
        JournaledStringTree<Dna5String> jst(seq, 100);

        String<unsigned> ids;
        appendValue(ids, 0);
        appendValue(ids, 3);
        appendValue(ids, 99);

        insert(jst, 5, 'C', ids, DeltaTypeSnp());
        insert(jst, 2, 3, ids, DeltaTypeDel());
        insert(jst, 5, "CGTA", ids, DeltaTypeIns());

        resize(jst, 150);
        auto& map = impl::member(jst, JstDeltaMapMember());
        for (auto& entry : map)
        {
            SEQAN_ASSERT_EQ(length(getDeltaCoverage(entry)), 150u);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[0], true);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[1], false);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[99], true);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[100], false);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[125], false);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[149], false);
        }

        resize(jst, 50);
        for (auto& entry : map)
        {
            SEQAN_ASSERT_EQ(length(getDeltaCoverage(entry)), 50u);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[0], true);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[1], false);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[3], true);
            SEQAN_ASSERT_EQ(getDeltaCoverage(entry)[49], false);
        }
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_empty)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);

    SEQAN_ASSERT(empty(jst));
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    insert(jst, 5, 'C', ids, DeltaTypeSnp());

    SEQAN_ASSERT(!empty(jst));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_size)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);

    SEQAN_ASSERT(size(jst) == 0u);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    insert(jst, 5, 'C', ids, DeltaTypeSnp());
    insert(jst, 8, 3, ids, DeltaTypeDel());

    SEQAN_ASSERT_EQ(size(jst), 3u);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_max_size)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 100);

    typedef typename Member<JournaledStringTree<Dna5String>, JstDeltaMapMember>::Type TDeltaMap;
    typedef typename Size<TDeltaMap>::Type TSize;

    SEQAN_ASSERT_EQ(maxSize(jst), std::numeric_limits<TSize>::max());
}

#endif // TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
