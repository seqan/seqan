// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
#include <seqan/sequence.h>
#include <seqan/journaled_string_tree.h>

#include "test_journaled_string_tree_traversal.h"

#include "test_delta_map.h"

using namespace seqan;

template <typename TJst, typename TExpectedType>
void _testJstContainerMF(TExpectedType &)
{
    typedef typename Container<TJst>::Type TContainer;
    bool res = IsSameType<TContainer, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

template <typename TJst, typename TExpectedType>
void _testJstGetStringSetMF(TExpectedType &)
{
    typedef typename GetStringSet<TJst>::Type TStringSet;

    bool res = IsSameType<TStringSet, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

template <typename TJst, typename TExpectedType>
void _testJstHostMF(TExpectedType &)
{
    typedef typename Host<TJst>::Type THost;

    bool res = IsSameType<THost, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_container_mf)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;

    TDeltaMap map;
    _testJstContainerMF<JournaledStringTree<TDeltaMap> >(map);

    const TDeltaMap constMap;
    _testJstContainerMF<const JournaledStringTree<TDeltaMap> >(constMap);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_get_string_tree_mf)
{
    StringSet<String<Dna, Journaled<Alloc<>, SortedArray> >, Owner<JournaledSet> > set;
    _testJstGetStringSetMF<JournaledStringTree<DeltaMap<unsigned, Dna> > >(set);

    const StringSet<String<Dna5, Journaled<Alloc<>, SortedArray> >, Owner<JournaledSet> > constSet;
    _testJstGetStringSetMF<const JournaledStringTree<DeltaMap<unsigned, Dna5> > >(constSet);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_host_mf)
{
    String<Dna, Alloc<> > host;
    _testJstHostMF<JournaledStringTree<DeltaMap<unsigned, Dna> > >(host);

    const String<Dna5, Alloc<> > constHost;
    _testJstHostMF<const JournaledStringTree<DeltaMap<unsigned, Dna5> > >(constHost);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_constructor)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;
    typedef Size<TJst>::Type TSize;

    {  // Default Constructor
        JournaledStringTree<TDeltaMap> jst;

        SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
        SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
        SEQAN_ASSERT_EQ(jst._activeBlock, 0u);
        SEQAN_ASSERT_EQ(jst._emptyJournal, true);
        SEQAN_ASSERT_EQ(length(jst._journalSet), 0u);
        SEQAN_ASSERT_EQ(length(jst._blockVPOffset), 0u);
        SEQAN_ASSERT_EQ(length(jst._activeBlockVPOffset), 0u);
    }

    {  // Parameterized Constructor
        String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);

        JournaledStringTree<TDeltaMap> jst(hostSeq, deltaMap);

        SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
        SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
        SEQAN_ASSERT_EQ(jst._activeBlock, 0u);
        SEQAN_ASSERT_EQ(jst._emptyJournal, true);
        SEQAN_ASSERT_EQ(&value(jst._container), &deltaMap);
        SEQAN_ASSERT_EQ(host(jst._journalSet), hostSeq);
        SEQAN_ASSERT_EQ(length(jst._journalSet), 8u);
        SEQAN_ASSERT_EQ(length(jst._blockVPOffset), 8u);
        SEQAN_ASSERT_EQ(length(jst._activeBlockVPOffset), 8u);
        for (unsigned i = 0; i < length(jst._journalSet); ++i)
        {
            SEQAN_ASSERT_EQ(host(jst._journalSet[i]), hostSeq);
            SEQAN_ASSERT_EQ(jst._blockVPOffset[i], 0);
            SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[i], 0);
        }
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_init)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;
    typedef Size<TJst>::Type TSize;

    JournaledStringTree<TDeltaMap> jst;

    SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
    SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 0u);
    SEQAN_ASSERT_EQ(jst._emptyJournal, true);
    SEQAN_ASSERT_EQ(length(jst._journalSet), 0u);
    SEQAN_ASSERT_EQ(length(jst._blockVPOffset), 0u);
    SEQAN_ASSERT_EQ(length(jst._activeBlockVPOffset), 0u);

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

     init(jst, hostSeq, deltaMap);

    SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
    SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 0u);
    SEQAN_ASSERT_EQ(jst._emptyJournal, true);
    SEQAN_ASSERT_EQ(&value(jst._container), &deltaMap);
    SEQAN_ASSERT_EQ(host(jst._journalSet), hostSeq);
    SEQAN_ASSERT_EQ(length(jst._journalSet), 8u);
    SEQAN_ASSERT_EQ(length(jst._blockVPOffset), 8u);
    SEQAN_ASSERT_EQ(length(jst._activeBlockVPOffset), 8u);
    for (unsigned i = 0; i < length(jst._journalSet); ++i)
    {
        SEQAN_ASSERT_EQ(host(jst._journalSet[i]), hostSeq);
        SEQAN_ASSERT_EQ(jst._blockVPOffset[i], 0);
        SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[i], 0);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_container)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";

    {  // Empty container
        TJst jst;
        SEQAN_ASSERT_EQ(empty(container(jst)), true);
    }

    {  // Set container
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);
        TJst jst(hostSeq, deltaMap);
        SEQAN_ASSERT_EQ(&container(jst), &deltaMap);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_string_set)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";

    {  // Empty container
        TJst jst;
        SEQAN_ASSERT_EQ(empty(stringSet(jst)), true);
    }

    {  // Set container
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);
        TJst jst(hostSeq, deltaMap);
        SEQAN_ASSERT_EQ(host(stringSet(jst)), hostSeq);
        SEQAN_ASSERT_EQ(length(stringSet(jst)), 8u);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_host)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";

    {  // Empty container
        TJst jst;
        SEQAN_ASSERT_EQ(empty(host(jst)), true);
    }

    {  // Set container
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);
        TJst jst(hostSeq, deltaMap);
        SEQAN_ASSERT_EQ(host(jst), hostSeq);
    }

    {  // Set container
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);
        TJst jst(hostSeq, deltaMap);
        SEQAN_ASSERT_EQ(host(jst), hostSeq);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_set_block_size)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;
    typedef typename Size<TJst>::Type TSize;

    JournaledStringTree<TDeltaMap> jst;
    SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);

    setBlockSize(jst, 3);
    SEQAN_ASSERT_EQ(jst._blockSize, 3u);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_block_size)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;
    typedef typename Size<TJst>::Type TSize;

    JournaledStringTree<TDeltaMap> jst;
    SEQAN_ASSERT_EQ(getBlockSize(jst), MaxValue<TSize>::VALUE);

    setBlockSize(jst, 3);
    SEQAN_ASSERT_EQ(getBlockSize(jst), 3u);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_full_journal_required)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";

    {  // Empty container
        TJst jst;
        SEQAN_ASSERT_EQ(fullJournalRequired(jst), true);
    }

    {  // Set container
        TDeltaMap deltaMap;
        _testDetlaMapfill(deltaMap);
        TJst jst(hostSeq, deltaMap);
        SEQAN_ASSERT_EQ(fullJournalRequired(jst), true);

        setBlockSize(jst, 3);
        SEQAN_ASSERT_EQ(fullJournalRequired(jst), false);
    }
}

template <typename TSet, typename TJst>
void
_createJournalString(TSet & set,
                    TJst const & jst,
                    unsigned beginPos,
                    unsigned endPos)
{
    typedef typename Container<TJst const>::Type TDeltaMap;
    typedef typename Iterator<TDeltaMap>::Type TIter;

    for (TIter it = begin(container(jst), Standard()) + endPos; it != begin(container(jst), Standard()) + beginPos; --it)
    {
        for (unsigned i = 0; i < length(set); ++i)
        {
            if (deltaCoverage(it - 1)[i] == true)
            {
                switch (deltaType(it - 1))
                {
                    case DELTA_TYPE_SNP:
                        replace(set[i], deltaPosition(it - 1), deltaPosition(it - 1) + 1, deltaValue(it - 1, DeltaTypeSnp()));
                        break;
                    case DELTA_TYPE_DEL:
                        erase(set[i], deltaPosition(it - 1), deltaPosition(it - 1) + deltaValue(it - 1, DeltaTypeDel()));
                        break;
                    case DELTA_TYPE_INS:
                        insert(set[i], deltaPosition(it - 1), deltaValue(it - 1, DeltaTypeIns()));
                        break;
                    case DELTA_TYPE_SV:
                        erase(set[i], deltaPosition(it - 1), deltaPosition(it - 1) + deltaValue(it - 1, DeltaTypeSV()).i1);
                        insert(set[i], deltaPosition(it - 1), deltaValue(it - 1, DeltaTypeSV()).i2);
                        break;
                }
            }
        }
    }
}

template <typename TJst>
void _testJournaledStringTreeJournalNextBlock(TJst & jst)
{
    typedef typename RemoveConst<typename Host<TJst>::Type >::Type TString;

    StringSet<TString> set;
    resize(set, 4, Exact());
    set[0] = host(jst);
    set[1] = host(jst);
    set[2] = host(jst);
    set[3] = host(jst);

    _createJournalString(set, jst, 0, length(container(jst)));
    
//    std::cerr << "[DEBUG] Output container entries\n";
//    for (unsigned i = 0; i < length(container(jst)._entries); ++i)
//        std::cerr << "Entry " << i << ": " << container(jst)._entries[i];
//    std::cerr << std::endl;

    {  // Journal full block
        create(jst, 5, Serial());

        TString test0 = stringSet(jst)[0];
        TString test1 = stringSet(jst)[1];
        TString test2 = stringSet(jst)[2];
        TString test3 = stringSet(jst)[3];

        SEQAN_ASSERT_EQ(set[0], test0);
        SEQAN_ASSERT_EQ(set[1], test1);
        SEQAN_ASSERT_EQ(set[2], test2);
        SEQAN_ASSERT_EQ(set[3], test3);
    }

    {  // Journal in blocks
        setBlockSize(jst, 3);
        reinit(jst);
        create(jst, 5, Serial());

//                 01234567890123  4567890123
        set[0] = "CACGTGGATCTGTAAAATGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                 012345678901234567890123
        set[1] = "CACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                012345678901234567890123
        set[2] = "ACGTGGATCTATACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                 012345678901    47890123
        set[3] = "CACGTGGATCTGTAAAATCGGACGGACTTGACGGGAGTACGAGCATCGACT";

        SEQAN_ASSERT_EQ(stringSet(jst)[0], set[0]);
        SEQAN_ASSERT_EQ(stringSet(jst)[1], set[1]);
        SEQAN_ASSERT_EQ(stringSet(jst)[2], set[2]);
        SEQAN_ASSERT_EQ(stringSet(jst)[3], set[3]);

        SEQAN_ASSERT_EQ(createNext(jst, 5, Serial()), true);

//                  012345678901234567890123
//        set[X] = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                012345678901234567890123
        set[0] = "ACGTGGATCTGTACTGACGGTAGGACTTGACGGGAGTACGAGCATCGACT";
//                0123456789012347890    123
        set[1] = "ACGTGGATCTGTACTCGGCACGTCGGACTTGACGGGAGTACGAGCATCGACT";
//                0123456789012347895
        set[2] = "ACGTGGATCTGTACTCGGCTTGACGGGAGTACGAGCATCGACT";
//                0123456789012347890    123
        set[3] = "ACGTGGATCTGTACTCGGCACGTCGGACTTGACGGGAGTACGAGCATCGACT";

        SEQAN_ASSERT_EQ(stringSet(jst)[0], set[0]);
        SEQAN_ASSERT_EQ(stringSet(jst)[1], set[1]);
        SEQAN_ASSERT_EQ(stringSet(jst)[2], set[2]);
        SEQAN_ASSERT_EQ(stringSet(jst)[3], set[3]);

        SEQAN_ASSERT_EQ(createNext(jst, 5, Serial()), true);

        //                012345678901234567890123
        set[0] = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                012345678901234567890    123
        set[1] = "ACGTGGATCTGTACTGACGGAACGTCGGACTTGACGGGAGTACGAGCATCGACT";
//                012345678901234567890123
        set[2] = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
//                012345678901234567890    123
        set[3] = "ACGTGGATCTGTACTGACGGAACGTCGGACTTGACGGGAGTACGAGCATCGACT";

        SEQAN_ASSERT_EQ(stringSet(jst)[0], set[0]);
        SEQAN_ASSERT_EQ(stringSet(jst)[1], set[1]);
        SEQAN_ASSERT_EQ(stringSet(jst)[2], set[2]);
        SEQAN_ASSERT_EQ(stringSet(jst)[3], set[3]);

        SEQAN_ASSERT_EQ(createNext(jst, 5, Serial()), false);
    }

}

SEQAN_DEFINE_TEST(test_journaled_string_tree_journal_next_block)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);
    setCoverageSize(deltaMap, 4);
    getDeltaCoverage(deltaMap._entries[0])[0] = true;
    getDeltaCoverage(deltaMap._entries[0])[1] = true;
    getDeltaCoverage(deltaMap._entries[0])[2] = false;
    getDeltaCoverage(deltaMap._entries[0])[3] = true;

    getDeltaCoverage(deltaMap._entries[1])[0] = false;
    getDeltaCoverage(deltaMap._entries[1])[1] = false;
    getDeltaCoverage(deltaMap._entries[1])[2] = true;
    getDeltaCoverage(deltaMap._entries[1])[3] = false;

    getDeltaCoverage(deltaMap._entries[2])[0] = true;
    getDeltaCoverage(deltaMap._entries[2])[1] = false;
    getDeltaCoverage(deltaMap._entries[2])[2] = false;
    getDeltaCoverage(deltaMap._entries[2])[3] = true;

    getDeltaCoverage(deltaMap._entries[3])[0] = false;
    getDeltaCoverage(deltaMap._entries[3])[1] = true;
    getDeltaCoverage(deltaMap._entries[3])[2] = true;
    getDeltaCoverage(deltaMap._entries[3])[3] = true;

    getDeltaCoverage(deltaMap._entries[4])[0] = true;
    getDeltaCoverage(deltaMap._entries[4])[1] = false;
    getDeltaCoverage(deltaMap._entries[4])[2] = false;
    getDeltaCoverage(deltaMap._entries[4])[3] = false;

    getDeltaCoverage(deltaMap._entries[5])[0] = false;
    getDeltaCoverage(deltaMap._entries[5])[1] = false;
    getDeltaCoverage(deltaMap._entries[5])[2] = true;
    getDeltaCoverage(deltaMap._entries[5])[3] = false;

    getDeltaCoverage(deltaMap._entries[6])[0] = false;
    getDeltaCoverage(deltaMap._entries[6])[1] = true;
    getDeltaCoverage(deltaMap._entries[6])[2] = false;
    getDeltaCoverage(deltaMap._entries[6])[3] = true;

    getDeltaCoverage(deltaMap._entries[7])[0] = false;
    getDeltaCoverage(deltaMap._entries[7])[1] = true;
    getDeltaCoverage(deltaMap._entries[7])[2] = false;
    getDeltaCoverage(deltaMap._entries[7])[3] = true;

    TJst jst(hostSeq, deltaMap);
    _testJournaledStringTreeJournalNextBlock(jst);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_reinit)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;
    typedef Size<TJst>::Type TSize;
    typedef GetStringSet<TJst>::Type TJournalSet;

    JournaledStringTree<TDeltaMap> jst;
    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    init(jst, hostSeq, deltaMap);
    create(jst, 5, Serial());

    TJournalSet journalSet = stringSet(jst);

    SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
    SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 1u);
    SEQAN_ASSERT_EQ(jst._emptyJournal, false);

    reinit(jst);

    SEQAN_ASSERT_EQ(jst._blockSize, MaxValue<TSize>::VALUE);
    SEQAN_ASSERT_EQ(jst._numBlocks, 1u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 0u);
    SEQAN_ASSERT_EQ(jst._emptyJournal, false);
    for (unsigned i = 0; i < length(jst._journalSet); ++i)
        SEQAN_ASSERT_EQ(jst._journalSet[i], journalSet[i]);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_local_to_global_pos)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    TJst jst;
    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);
    init(jst, hostSeq, deltaMap);

    setCoverageSize(deltaMap, 4);
    getDeltaCoverage(deltaMap._entries[0])[0] = true;
    getDeltaCoverage(deltaMap._entries[0])[1] = true;
    getDeltaCoverage(deltaMap._entries[0])[2] = false;
    getDeltaCoverage(deltaMap._entries[0])[3] = true;

    getDeltaCoverage(deltaMap._entries[1])[0] = false;
    getDeltaCoverage(deltaMap._entries[1])[1] = false;
    getDeltaCoverage(deltaMap._entries[1])[2] = true;
    getDeltaCoverage(deltaMap._entries[1])[3] = false;

    getDeltaCoverage(deltaMap._entries[2])[0] = true;
    getDeltaCoverage(deltaMap._entries[2])[1] = false;
    getDeltaCoverage(deltaMap._entries[2])[2] = false;
    getDeltaCoverage(deltaMap._entries[2])[3] = true;

    getDeltaCoverage(deltaMap._entries[3])[0] = false;
    getDeltaCoverage(deltaMap._entries[3])[1] = true;
    getDeltaCoverage(deltaMap._entries[3])[2] = true;
    getDeltaCoverage(deltaMap._entries[3])[3] = true;

    getDeltaCoverage(deltaMap._entries[4])[0] = true;
    getDeltaCoverage(deltaMap._entries[4])[1] = false;
    getDeltaCoverage(deltaMap._entries[4])[2] = false;
    getDeltaCoverage(deltaMap._entries[4])[3] = false;

    getDeltaCoverage(deltaMap._entries[5])[0] = false;
    getDeltaCoverage(deltaMap._entries[5])[1] = false;
    getDeltaCoverage(deltaMap._entries[5])[2] = true;
    getDeltaCoverage(deltaMap._entries[5])[3] = false;

    getDeltaCoverage(deltaMap._entries[6])[0] = false;
    getDeltaCoverage(deltaMap._entries[6])[1] = true;
    getDeltaCoverage(deltaMap._entries[6])[2] = false;
    getDeltaCoverage(deltaMap._entries[6])[3] = true;

    getDeltaCoverage(deltaMap._entries[7])[0] = false;
    getDeltaCoverage(deltaMap._entries[7])[1] = true;
    getDeltaCoverage(deltaMap._entries[7])[2] = false;
    getDeltaCoverage(deltaMap._entries[7])[3] = true;

    setBlockSize(jst, 3);
    create(jst, 5, Serial());

    SEQAN_ASSERT_EQ(jst._blockSize, 3u);
    SEQAN_ASSERT_EQ(jst._numBlocks, 3u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 1u);

    SEQAN_ASSERT_EQ(jst._blockVPOffset[0], 0);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[0], 3);
    SEQAN_ASSERT_EQ(localToGlobalPos(0, 0, jst), 0u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[1], 0);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[1], 1);
    SEQAN_ASSERT_EQ(localToGlobalPos(0, 1, jst), 0u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[2], 0);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[2], 0);
    SEQAN_ASSERT_EQ(localToGlobalPos(0, 2, jst), 0u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[3], 0);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[3], 3);
    SEQAN_ASSERT_EQ(localToGlobalPos(0, 3, jst), 0u);

    createNext(jst, 5, Serial());

    SEQAN_ASSERT_EQ(jst._blockSize, 3u);
    SEQAN_ASSERT_EQ(jst._numBlocks, 3u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 2u);

    SEQAN_ASSERT_EQ(jst._blockVPOffset[0], 3);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[0], 0);
    SEQAN_ASSERT_EQ(localToGlobalPos(2, 0, jst), 5u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[1], 1);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[1], -2);
    SEQAN_ASSERT_EQ(localToGlobalPos(2, 1, jst), 3u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[2], 0);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[2], -7);
    SEQAN_ASSERT_EQ(localToGlobalPos(2, 2, jst), 2u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[3], 3);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[3], -2);
    SEQAN_ASSERT_EQ(localToGlobalPos(2, 3, jst), 5u);

    createNext(jst, 5, Serial());

    SEQAN_ASSERT_EQ(jst._blockSize, 3u);
    SEQAN_ASSERT_EQ(jst._numBlocks, 3u);
    SEQAN_ASSERT_EQ(jst._activeBlock, 3u);

    SEQAN_ASSERT_EQ(jst._blockVPOffset[0], 3);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[0], 0);
    SEQAN_ASSERT_EQ(localToGlobalPos(20, 0, jst), 23u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[1], -1);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[1], 4);
    SEQAN_ASSERT_EQ(localToGlobalPos(20, 1, jst), 19u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[2], -7);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[2], 0);
    SEQAN_ASSERT_EQ(localToGlobalPos(20, 2, jst), 13u);
    SEQAN_ASSERT_EQ(jst._blockVPOffset[3], 1);
    SEQAN_ASSERT_EQ(jst._activeBlockVPOffset[3], 4);
    SEQAN_ASSERT_EQ(localToGlobalPos(20, 3, jst), 21u);
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_save_open)
{
    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    String<Dna> hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);
    setCoverageSize(deltaMap, 4);
    getDeltaCoverage(deltaMap._entries[0])[0] = true;
    getDeltaCoverage(deltaMap._entries[0])[1] = true;
    getDeltaCoverage(deltaMap._entries[0])[2] = false;
    getDeltaCoverage(deltaMap._entries[0])[3] = true;

    getDeltaCoverage(deltaMap._entries[1])[0] = false;
    getDeltaCoverage(deltaMap._entries[1])[1] = false;
    getDeltaCoverage(deltaMap._entries[1])[2] = true;
    getDeltaCoverage(deltaMap._entries[1])[3] = false;

    getDeltaCoverage(deltaMap._entries[2])[0] = true;
    getDeltaCoverage(deltaMap._entries[2])[1] = false;
    getDeltaCoverage(deltaMap._entries[2])[2] = false;
    getDeltaCoverage(deltaMap._entries[2])[3] = true;

    getDeltaCoverage(deltaMap._entries[3])[0] = false;
    getDeltaCoverage(deltaMap._entries[3])[1] = true;
    getDeltaCoverage(deltaMap._entries[3])[2] = true;
    getDeltaCoverage(deltaMap._entries[3])[3] = true;

    getDeltaCoverage(deltaMap._entries[4])[0] = true;
    getDeltaCoverage(deltaMap._entries[4])[1] = false;
    getDeltaCoverage(deltaMap._entries[4])[2] = false;
    getDeltaCoverage(deltaMap._entries[4])[3] = false;

    getDeltaCoverage(deltaMap._entries[5])[0] = false;
    getDeltaCoverage(deltaMap._entries[5])[1] = false;
    getDeltaCoverage(deltaMap._entries[5])[2] = true;
    getDeltaCoverage(deltaMap._entries[5])[3] = false;

    getDeltaCoverage(deltaMap._entries[6])[0] = false;
    getDeltaCoverage(deltaMap._entries[6])[1] = true;
    getDeltaCoverage(deltaMap._entries[6])[2] = false;
    getDeltaCoverage(deltaMap._entries[6])[3] = true;

    getDeltaCoverage(deltaMap._entries[7])[0] = false;
    getDeltaCoverage(deltaMap._entries[7])[1] = true;
    getDeltaCoverage(deltaMap._entries[7])[2] = false;
    getDeltaCoverage(deltaMap._entries[7])[3] = true;

    TJst jst(hostSeq, deltaMap);

    CharString testDir = SEQAN_PARENT_PATH(__FILE__);

    CharString wrongPath = "wrongPath";
    append(wrongPath, testDir);
    append(wrongPath, "unknown.gdf");

    CharString goldRefFile = testDir;
    append(goldRefFile, "gold_ref.fa");

    CharString tmpTestDir = SEQAN_PARENT_PATH(SEQAN_TEMP_FILENAME());
    CharString tmpGdfFile = tmpTestDir;
    append(tmpGdfFile, "testGdf.gdf");

    CharString tmpRefFile = tmpTestDir;
    append(tmpRefFile, "testRef.fa");

    GdfHeader header;
    header.referenceMode = GdfIOMode::SAVE_REFERENCE_MODE_ENABLED;
#ifdef SEQAN_EXCEPTIONS  // Only test if exceptions are enabled.

    SEQAN_TRY  // Try save with invalid name store.
    {
        save(jst, header, tmpGdfFile);
    }
    SEQAN_CATCH(GdfIOException e)
    {
        SEQAN_ASSERT(isEqual(e.what(), "Gdf_IO_Exception: (Too few sequence names - Needed 4 but 0 were provided!)"));
    }

    appendValue(header.nameStore, "Seq1");
    appendValue(header.nameStore, "Seq2");
    appendValue(header.nameStore, "Seq3");
    appendValue(header.nameStore, "Seq4");

    header.referenceFilename = "unknown";
    SEQAN_TRY  // Try save with unknown reference file.
    {
        save(jst, header, tmpGdfFile);
    }
    SEQAN_CATCH(GdfIOException e)
    {
        std::stringstream errMessage;
        errMessage << "Gdf_IO_Exception: (Cannot open file: unknown!)";
        SEQAN_ASSERT_EQ(e.what(), errMessage.str());
    }

    SEQAN_TRY  // Try save with wrong path.
    {
        save(jst, header, wrongPath);
    }
    SEQAN_CATCH(GdfIOException e)
    {
        std::stringstream errMessage;
        errMessage << "Gdf_IO_Exception: (Cannot open file: "<< wrongPath << "!)";
        SEQAN_ASSERT_EQ(e.what(), errMessage.str());
    }

    header.referenceId = "reference";
    header.referenceFilename = tmpRefFile;
    SEQAN_TRY  // Try save with everything set correctly.
    {
        save(jst, header, tmpGdfFile);
    }
    SEQAN_CATCH(Exception e)
    {
        SEQAN_FAIL("Error while save!");
    }
#else  // SEQAN_EXCEPTIONS
    save(jst, header, tmpGdfFile);
#endif  // SEQAN_EXCEPTIONS

    // Compare the reference sequences.
    SEQAN_ASSERT(_compareTextFiles(toCString(tmpRefFile), toCString(goldRefFile)));

    CharString refId;
    CharString refFilename;
    String<CharString> seqIds;

#ifdef SEQAN_EXCEPTIONS

    {  // Try load from wrong file.
        SEQAN_TRY
        {
            JournaledStringTree<TDeltaMap> jstLoaded;
            GdfHeader headerLoaded;
            open(jstLoaded, headerLoaded, wrongPath);
        }
        SEQAN_CATCH(GdfIOException e)
        {
            std::stringstream errMessage;
            errMessage << "Gdf_IO_Exception: (Unknown file <" << wrongPath << ">!)";
            SEQAN_ASSERT_EQ(e.what(), errMessage.str());
        }
    }
    {  // Try load from file with unknown reference file.
        SEQAN_TRY
        {
            JournaledStringTree<TDeltaMap> jstLoaded;
            GdfHeader headerLoaded;
            CharString unknownRefExample = testDir;
            append(unknownRefExample, "unknownRefExample.gdf");
            open(jstLoaded, headerLoaded, unknownRefExample);
        }
        SEQAN_CATCH(GdfIOException e)
        {
            SEQAN_ASSERT(isEqual(e.what(), "Gdf_IO_Exception: (Unknown file </wrong/path/to/reference.fa>!)"));
        }
    }
    {  // Try load from file with wrong crc or refId.
        SEQAN_TRY
        {
            CharString wrongRefFile = testDir;
            append(wrongRefFile, "wrong_ref.fa");

            CharString wrongCrcExample = tmpTestDir;
            append(wrongCrcExample, "wrongCrcExample.gdf");


            GdfHeader testHeader;
            testHeader.referenceId = "reference";
            testHeader.referenceFilename = wrongRefFile;
            testHeader.nameStore = header.nameStore;
            testHeader.referenceMode = GdfIOMode::SAVE_REFERENCE_MODE_DISABLED;

            SEQAN_TRY
            {
                save(jst, testHeader, wrongCrcExample);
            }
            SEQAN_CATCH(Exception e)
            {
                SEQAN_ASSERT_FAIL("Error during save!");
            }

            JournaledStringTree<TDeltaMap> jstLoaded;
            GdfHeader headerLoaded;
            open(jstLoaded, headerLoaded, wrongCrcExample);
        }
        SEQAN_CATCH(GdfIOException e)
        {
#ifdef __SSE4_2__
            SEQAN_ASSERT(isEqual(e.what(), "Gdf_IO_Exception: (The id of the reference is \'another ref\' but should be \'reference\' and the crc is \'0\' but should be \'0\'!)"));
#else  // __SSE4_2__
            SEQAN_ASSERT(isEqual(e.what(), "Gdf_IO_Exception: (The id of the reference is \'another ref\' but should be \'reference\' and the crc is \'0\' but should be \'0\'!)"));
#endif // __SSE4_2__
        }
    }

    JournaledStringTree<TDeltaMap> jstLoaded;
    {  // Try load correct file.
        SEQAN_TRY
        {
            GdfHeader headerLoaded;
            open(jstLoaded, headerLoaded, tmpGdfFile);
        }
        SEQAN_CATCH(Exception e)
        {
            std::cerr << e.what() << std::endl;
            SEQAN_FAIL("Error while open!");
        }
    }
#else
    GdfHeader headerLoaded;
    open(jstLoaded, headerLoaded, tmpGdfFile);
#endif

    _testJournaledStringTreeJournalNextBlock(jstLoaded);
}

// ----------------------------------------------------------------------------
// Test Traversal Concept.
// ----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_journaled_string_tree_jst_traversal_concept)
{
    using namespace seqan;

    typedef DeltaMap<unsigned, Dna> TDeltaMap;
    typedef JournaledStringTree<TDeltaMap> TJst;

    typedef DummyCaller_<TJst> TDummy;

    bool res =  Is<JstTraversalConcept<TDummy> >::VALUE;
    SEQAN_ASSERT_EQ(res, true);
}

#endif // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
