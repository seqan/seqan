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
// Test for the delta map.
// ==========================================================================

#ifndef EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
#define EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

struct TestDeltaMapConfig
{
    typedef unsigned    TDeltaPos;
    typedef Dna         TSnpValue;
    typedef DnaString   TInsValue;
    typedef unsigned    TDelValue;
    typedef Pair<unsigned, DnaString> TSVValue;
};

template <typename TEntries, typename TStore>
inline void createMock(TEntries & entries, TStore & store)
{
    typedef typename DeltaValue<TStore, DeltaTypeSV>::Type TSV;
    typedef typename Value<TEntries>::Type TEntry;
    typedef typename DeltaRecord<TEntry>::Type TDeltaRecord;
    typedef typename DeltaCoverage<TEntry>::Type TCoverage;

    TCoverage cov1;
    resize(cov1, 10);
    for (unsigned i = 0; i < length(cov1); ++i)
        cov1[i] = (!(i % 3)) ? true : false;

    TCoverage cov2 = cov1;
    for (unsigned i = 0; i < length(cov2); ++i)
        cov2[i] = (!(i % 5)) ? true : false;

    TCoverage cov3 = cov1;
    for (unsigned i = 0; i < length(cov3); ++i)
        cov3[i] = (!(i % 2)) ? true : false;

    appendValue(entries, TEntry(0, TDeltaRecord(DELTA_TYPE_SNP, 1), cov1));
    appendValue(entries, TEntry(1, TDeltaRecord(DELTA_TYPE_DEL, 2), cov2));
    appendValue(entries, TEntry(1, TDeltaRecord(DELTA_TYPE_SV, 0), cov1));
    appendValue(entries, TEntry(2, TDeltaRecord(DELTA_TYPE_SNP, 0), cov2));
    appendValue(entries, TEntry(4, TDeltaRecord(DELTA_TYPE_INS, 0), cov3));
    appendValue(entries, TEntry(5, TDeltaRecord(DELTA_TYPE_DEL, 0), cov3));
    appendValue(entries, TEntry(20, TDeltaRecord(DELTA_TYPE_DEL, 1), cov2));

    appendValue(store._snpData, 'C');
    appendValue(store._snpData, 'A');
    appendValue(store._insData, "ACGT");
    appendValue(store._delData, 1);
    appendValue(store._delData, 2);
    appendValue(store._delData, 3);
    appendValue(store._svData, TSV(2, "TGAT"));
}

SEQAN_DEFINE_TEST(test_delta_map_insert)
{

    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Member<TDeltaMap, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Member<TDeltaMap, DeltaMapStoreMember>::Type TDeltaStore;

    TEntries testEntries;
    TDeltaStore store;
    createMock(testEntries, store);

    TCoverage cov1;
    resize(cov1, 10);
    for (unsigned i = 0; i < length(cov1); ++i)
        cov1[i] = (!(i % 3)) ? true : false;

    TCoverage cov2 = cov1;
    for (unsigned i = 0; i < length(cov2); ++i)
        cov2[i] = (!(i % 5)) ? true : false;

    TCoverage cov3 = cov1;
    for (unsigned i = 0; i < length(cov3); ++i)
        cov3[i] = (!(i % 2)) ? true : false;

    TDeltaMap deltaMap;
    SEQAN_ASSERT(insert(deltaMap,  2, store._snpData[0], cov2, DeltaTypeSnp()));
    SEQAN_ASSERT(insert(deltaMap,  5, store._delData[0], cov3, DeltaTypeDel()));
    SEQAN_ASSERT(insert(deltaMap,  1,  store._svData[0], cov1, DeltaTypeSV()));
    SEQAN_ASSERT(insert(deltaMap,  0, store._snpData[1], cov1, DeltaTypeSnp()));
    SEQAN_ASSERT(insert(deltaMap,  4, store._insData[0], cov3, DeltaTypeIns()));
    SEQAN_ASSERT(insert(deltaMap, 20, store._delData[1], cov2, DeltaTypeDel()));
    SEQAN_ASSERT(insert(deltaMap,  1, store._delData[2], cov2, DeltaTypeDel()));
    SEQAN_ASSERT_NOT(insert(deltaMap,  1, store._delData[2], cov2, DeltaTypeDel()));

    for (unsigned i = 0; i < length(testEntries); ++i)
    {
        SEQAN_ASSERT_EQ(deltaMap._entries[i], testEntries[i]);
        unsigned recordPos = getDeltaRecord(deltaMap._entries[i]).i2;
        switch (static_cast<DeltaType>(getDeltaRecord(deltaMap._entries[i]).i1))
        {
            case DELTA_TYPE_SNP: SEQAN_ASSERT_EQ(deltaValue(deltaMap._deltaStore, recordPos, DeltaTypeSnp()),
                                                 deltaValue(store, recordPos, DeltaTypeSnp()));
                break;
            case DELTA_TYPE_INS: SEQAN_ASSERT_EQ(deltaValue(deltaMap._deltaStore, recordPos, DeltaTypeIns()),
                                                 deltaValue(store, recordPos, DeltaTypeIns()));
                break;
            case DELTA_TYPE_DEL: SEQAN_ASSERT_EQ(deltaValue(deltaMap._deltaStore, recordPos, DeltaTypeDel()),
                                                 deltaValue(store, recordPos, DeltaTypeDel()));
                break;
            case DELTA_TYPE_SV: SEQAN_ASSERT_EQ(deltaValue(deltaMap._deltaStore, recordPos, DeltaTypeSV()),
                                                deltaValue(store, recordPos, DeltaTypeSV()));
                break;
        }
    }
}

SEQAN_DEFINE_TEST(test_delta_map_erase)
{
    DeltaMap<TestDeltaMapConfig> deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    SEQAN_ASSERT(erase(deltaMap,  2, DeltaTypeSnp()));
    SEQAN_ASSERT(erase(deltaMap,  5, DeltaTypeDel()));
    SEQAN_ASSERT(erase(deltaMap,  1, DeltaTypeSV()));
    SEQAN_ASSERT(erase(deltaMap,  0, DeltaTypeSnp()));
    SEQAN_ASSERT(erase(deltaMap,  4, DeltaTypeIns()));
    SEQAN_ASSERT(erase(deltaMap, 20, DeltaTypeDel()));
    SEQAN_ASSERT(erase(deltaMap,  1, DeltaTypeDel()));
    SEQAN_ASSERT_EQ(length(deltaMap._entries), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._snpData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._delData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._insData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._svData), 0u);
    SEQAN_ASSERT_NOT(erase(deltaMap,  1, DeltaTypeDel()));
}

SEQAN_DEFINE_TEST(test_delta_map_find)
{
    DeltaMap<TestDeltaMapConfig> deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    SEQAN_ASSERT_EQ(value(find(deltaMap,  2, DeltaTypeSnp())), deltaMap._entries[3]);
    SEQAN_ASSERT_EQ(value(find(deltaMap,  1, DeltaTypeSV())), deltaMap._entries[2]);
    SEQAN_ASSERT_EQ(value(find(deltaMap,  20, DeltaTypeDel())), deltaMap._entries[6]);
    SEQAN_ASSERT(find(deltaMap,  1, DeltaTypeIns()) == end(deltaMap, Standard()));
    SEQAN_ASSERT(find(deltaMap,  6, DeltaTypeSnp()) == end(deltaMap, Standard()));
}

SEQAN_DEFINE_TEST(test_delta_map_size)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(size(deltaMap), 0u);
    createMock(deltaMap._entries, deltaMap._deltaStore);
    SEQAN_ASSERT_EQ(size(deltaMap), 7u);
}

SEQAN_DEFINE_TEST(test_delta_map_empty)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(empty(deltaMap), true);
    createMock(deltaMap._entries, deltaMap._deltaStore);
    SEQAN_ASSERT_EQ(empty(deltaMap), false);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_copy_constructor)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator itMap = begin(deltaMap, Standard());

    TIterator copyItMap = itMap;
    TConstIterator constCopyItMap = copyItMap;
    TConstIterator constCopyItMap2 = constCopyItMap;

    SEQAN_ASSERT_EQ(copyItMap._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(copyItMap._mapIter, itMap._mapIter);

    SEQAN_ASSERT_EQ(constCopyItMap._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(constCopyItMap._mapIter, itMap._mapIter);

    SEQAN_ASSERT_EQ(constCopyItMap2._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(constCopyItMap2._mapIter, itMap._mapIter);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_assign)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator itMap = begin(deltaMap, Standard());

    TIterator copyItMap;
    TConstIterator constCopyItMap;
    TConstIterator constCopyItMap2;

    copyItMap = itMap;

    constCopyItMap = copyItMap;
    constCopyItMap2 = constCopyItMap;

    SEQAN_ASSERT_EQ(copyItMap._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(copyItMap._mapIter, itMap._mapIter);

    SEQAN_ASSERT_EQ(constCopyItMap._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(constCopyItMap._mapIter, itMap._mapIter);

    SEQAN_ASSERT_EQ(constCopyItMap2._mapPtr, itMap._mapPtr);
    SEQAN_ASSERT_EQ(constCopyItMap2._mapIter, itMap._mapIter);
}

template <typename TMap>
inline void
_testDeltaMapIterator(TMap & deltaMap)
{
    typedef typename Iterator<TMap, Standard>::Type TIterator;

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it), 0u);

    unsigned counter = 0;
    for (; it != end(deltaMap, Standard()); ++it, ++counter)
        SEQAN_ASSERT_EQ(value(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, size(deltaMap));

    for (; it != begin(deltaMap, Standard()); --it, --counter)
        SEQAN_ASSERT_EQ(value(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

    for (; !(it == end(deltaMap, Standard())); it++, ++counter)
        SEQAN_ASSERT_EQ(value(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, size(deltaMap));

    for (; !(it == begin(deltaMap, Standard())); it--, --counter)
        SEQAN_ASSERT_EQ(value(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

}

SEQAN_DEFINE_TEST(test_delta_map_iterator)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);
    const TDeltaMap deltaMap2 = deltaMap;

    _testDeltaMapIterator(deltaMap);
    _testDeltaMapIterator(deltaMap2);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_value)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(value(it), deltaMap._entries[0]);
    SEQAN_ASSERT_EQ(value(++it), deltaMap._entries[1]);
    SEQAN_ASSERT_EQ(value(++it), deltaMap._entries[2]);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(value(it2++), deltaMap._entries[0]);
    SEQAN_ASSERT_EQ(value(it2++), deltaMap._entries[1]);
    SEQAN_ASSERT_EQ(value(it2), deltaMap._entries[2]);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_type)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaType(it),   DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_INS);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_INS);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_position)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 2u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 4u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 2u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 4u);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_value)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
    typedef DeltaValue<TDeltaMap, DeltaTypeSV>::Type TSV;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeSnp()), TSnp('A'));
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeDel()), TDel(3));
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeSV()), TSV(2, "TGAT"));
    ++it;
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeIns()), TIns("ACGT"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeSnp()), TSnp('A'));
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeDel()), TDel(3));
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeSV()), TSV(2, "TGAT"));
    ++it2;
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeIns()), TIns("ACGT"));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_coverage)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaCoverage<TDeltaMap>::Type TCoverage;

    TDeltaMap deltaMap;
    createMock(deltaMap._entries, deltaMap._deltaStore);

    TCoverage cov1;
    resize(cov1, 10);
    for (unsigned i = 0; i < length(cov1); ++i)
        cov1[i] = (!(i % 3)) ? true : false;

    TCoverage cov2 = cov1;
    for (unsigned i = 0; i < length(cov2); ++i)
        cov2[i] = (!(i % 5)) ? true : false;

    TCoverage cov3 = cov1;
    for (unsigned i = 0; i < length(cov3); ++i)
        cov3[i] = (!(i % 2)) ? true : false;

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaCoverage(it++), cov1);
    SEQAN_ASSERT_EQ(deltaCoverage(it++), cov2);
    SEQAN_ASSERT_EQ(deltaCoverage(it++), cov1);
    SEQAN_ASSERT_EQ(deltaCoverage(it++), cov2);
    SEQAN_ASSERT_EQ(deltaCoverage(it++), cov3);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaCoverage(it2++), cov1);
    SEQAN_ASSERT_EQ(deltaCoverage(it2++), cov2);
    SEQAN_ASSERT_EQ(deltaCoverage(it2++), cov1);
    SEQAN_ASSERT_EQ(deltaCoverage(it2++), cov2);
    SEQAN_ASSERT_EQ(deltaCoverage(it2++), cov3);
}


#endif  // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
