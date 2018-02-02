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
// Test for the delta map.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_

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

template <typename TCoverage, typename TSize>
inline void
_getCoverage(TCoverage & cov, TSize size, TSize remainder)
{
    resize(cov, size);
    for (unsigned i = 0; i < length(cov); ++i)
        cov[i] = (!(i % remainder)) ? true : false;
}

template <typename TConfig, typename TSpec>
inline void createMock(DeltaMap<TConfig, TSpec> & map)
{
    typedef DeltaMap<TConfig, TSpec> TDeltMap;
    typedef typename Member<TDeltMap, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Member<TDeltMap, DeltaMapStoreMember>::Type TStore;
    typedef typename DeltaValue<TStore, DeltaTypeSV>::Type TSV;
    typedef typename Value<TEntries>::Type TEntry;
    typedef typename DeltaCoverage<TEntry>::Type TCoverage;

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);
    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    TCoverage cov3;
    _getCoverage(cov3, 10, 2);

    insert(map, 5, 1, cov3, DeltaTypeDel());
    insert(map, 1, TSV(2, "TGAT"), cov1, DeltaTypeSV());
    insert(map, 4, "ACGT", cov3, DeltaTypeIns());
    insert(map, 2, 'C', cov2, DeltaTypeSnp());
    insert(map, 20, 2, cov2, DeltaTypeDel());
    insert(map, 0, 'A', cov1, DeltaTypeSnp());
    insert(map, 1, 3, cov2, DeltaTypeDel());
}

SEQAN_DEFINE_TEST(test_delta_map_insert)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Member<TDeltaMap, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Value<TEntries>::Type TEntry;
    typedef DeltaRecord<TEntry>::Type TRecord;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);
    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    TCoverage cov3;
    _getCoverage(cov3, 10, 2);

    insert(deltaMap,  2, deltaMap._deltaStore._snpData[0], cov2, DeltaTypeSnp());
    insert(deltaMap,  1,  deltaMap._deltaStore._svData[0], cov1, DeltaTypeSV());
    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov2, DeltaTypeDel());
    insert(deltaMap,  3, 1, cov1, DeltaTypeDel());
    insert(deltaMap,  3, 4, cov1, DeltaTypeDel());

    SEQAN_ASSERT_EQ(deltaMap._entries[0], TEntry(0, TRecord(DELTA_TYPE_SNP, 1), cov1, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[1], TEntry(1, TRecord(DELTA_TYPE_DEL, 2), cov2, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[2], TEntry(1, TRecord(DELTA_TYPE_DEL, 3), cov2, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[3], TEntry(1, TRecord(DELTA_TYPE_SV, 0), cov1, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[4], TEntry(1, TRecord(DELTA_TYPE_SV, 1), cov1, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[5], TEntry(2, TRecord(DELTA_TYPE_SNP, 0), cov2, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[6], TEntry(2, TRecord(DELTA_TYPE_SNP, 2), cov2, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[7], TEntry(2, TRecord(DELTA_TYPE_SV, 0), cov1, DeltaEndType::IS_RIGHT));
    SEQAN_ASSERT_EQ(deltaMap._entries[8], TEntry(2, TRecord(DELTA_TYPE_SV, 1), cov1, DeltaEndType::IS_RIGHT));
    SEQAN_ASSERT_EQ(deltaMap._entries[9], TEntry(3, TRecord(DELTA_TYPE_DEL, 2), cov2, DeltaEndType::IS_RIGHT));
    SEQAN_ASSERT_EQ(deltaMap._entries[10], TEntry(3, TRecord(DELTA_TYPE_DEL, 3), cov2, DeltaEndType::IS_RIGHT));
    SEQAN_ASSERT_EQ(deltaMap._entries[11], TEntry(3, TRecord(DELTA_TYPE_DEL, 4), cov1, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[12], TEntry(3, TRecord(DELTA_TYPE_DEL, 5), cov1, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[13], TEntry(4, TRecord(DELTA_TYPE_INS, 0), cov3, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[14], TEntry(5, TRecord(DELTA_TYPE_DEL, 0), cov3, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(deltaMap._entries[15], TEntry(6, TRecord(DELTA_TYPE_DEL, 5), cov1, DeltaEndType::IS_RIGHT));
    SEQAN_ASSERT_EQ(deltaMap._entries[16], TEntry(20, TRecord(DELTA_TYPE_DEL, 1), cov2, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT_EQ(deltaMap._entries[17], TEntry(21, TRecord(DELTA_TYPE_DEL, 1), cov2, DeltaEndType::IS_RIGHT));
}

SEQAN_DEFINE_TEST(test_delta_map_erase)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov2, DeltaTypeDel());

    SEQAN_ASSERT_EQ(erase(deltaMap,  2, DeltaTypeSnp()), 1u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  5, DeltaTypeDel()), 1u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  1, DeltaTypeSV()), 2u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  0, DeltaTypeSnp()), 1u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  4, DeltaTypeIns()), 1u);
    SEQAN_ASSERT_EQ(erase(deltaMap, 20, DeltaTypeDel()), 2u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  1, DeltaTypeDel()), 4u);
    SEQAN_ASSERT_EQ(length(deltaMap._entries), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._snpData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._delData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._insData), 0u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaStore._svData), 0u);
    SEQAN_ASSERT_EQ(erase(deltaMap,  1, DeltaTypeDel()), 0u);
}

SEQAN_DEFINE_TEST(test_delta_map_lower_bound)
{

    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Member<TDeltaMap, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Value<TEntries>::Type TEntry;
    typedef DeltaRecord<TEntry>::Type TRecord;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);
    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    TCoverage cov3;
    _getCoverage(cov3, 10, 2);

    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov1, DeltaTypeDel());
    SEQAN_ASSERT_EQ(*lowerBound(deltaMap, 7, DeltaTypeSnp()), TEntry(20, TRecord(DELTA_TYPE_DEL, 1), cov2, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT(lowerBound(deltaMap, 22, DeltaTypeSnp()) == end(deltaMap, Standard()));
    SEQAN_ASSERT_EQ(*lowerBound(deltaMap, 4, DeltaTypeIns()), TEntry(4, TRecord(DELTA_TYPE_INS, 0), cov3, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(*lowerBound(deltaMap, 1, DeltaTypeDel()), TEntry(1, TRecord(DELTA_TYPE_DEL, 2), cov2, DeltaEndType::IS_LEFT));
}

SEQAN_DEFINE_TEST(test_delta_map_upper_bound)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    typedef typename Member<TDeltaMap, DeltaMapEntriesMember>::Type TEntries;
    typedef typename Value<TEntries>::Type TEntry;
    typedef DeltaRecord<TEntry>::Type TRecord;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);
    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    TCoverage cov3;
    _getCoverage(cov3, 10, 2);

    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov1, DeltaTypeDel());
    SEQAN_ASSERT_EQ(*upperBound(deltaMap, 7, DeltaTypeSnp()), TEntry(20, TRecord(DELTA_TYPE_DEL, 1), cov2, DeltaEndType::IS_LEFT));
    SEQAN_ASSERT(upperBound(deltaMap, 22, DeltaTypeSnp()) == end(deltaMap, Standard()));
    SEQAN_ASSERT_EQ(*upperBound(deltaMap, 4, DeltaTypeIns()), TEntry(5, TRecord(DELTA_TYPE_DEL, 0), cov3, DeltaEndType::IS_BOTH));
    SEQAN_ASSERT_EQ(*upperBound(deltaMap, 1, DeltaTypeDel()), TEntry(1, TRecord(DELTA_TYPE_SV, 0), cov1, DeltaEndType::IS_LEFT));

}

SEQAN_DEFINE_TEST(test_delta_map_count)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);

    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov1, DeltaTypeDel());

    SEQAN_ASSERT_EQ(count(deltaMap, 7, DeltaTypeSnp()), 0u);
    SEQAN_ASSERT_EQ(count(deltaMap, 4, DeltaTypeIns()), 1u);
    SEQAN_ASSERT_EQ(count(deltaMap, 1, DeltaTypeDel()), 2u);
}

SEQAN_DEFINE_TEST(test_delta_map_equal_range)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);

    insert(deltaMap,  1, deltaMap._deltaStore._delData[2], cov1, DeltaTypeDel());

    auto range = equalRange(deltaMap, 7, DeltaTypeSnp());
    SEQAN_ASSERT(range.i1 == range.i2);

    range = equalRange(deltaMap, 22, DeltaTypeSnp());
    SEQAN_ASSERT(range.i1 == end(deltaMap, Standard()));
    SEQAN_ASSERT(range.i2 == end(deltaMap, Standard()));

    range = equalRange(deltaMap, 4, DeltaTypeIns());
    SEQAN_ASSERT_EQ(*range.i1, *lowerBound(deltaMap, 4, DeltaTypeIns()));
    SEQAN_ASSERT_EQ(*range.i2, *upperBound(deltaMap, 4, DeltaTypeIns()));

    range = equalRange(deltaMap, 1, DeltaTypeDel());
    SEQAN_ASSERT_EQ(*range.i1, *lowerBound(deltaMap, 1, DeltaTypeDel()));
    SEQAN_ASSERT_EQ(*range.i2, *upperBound(deltaMap, 1, DeltaTypeDel()));
}

SEQAN_DEFINE_TEST(test_delta_map_find)
{
    DeltaMap<TestDeltaMapConfig> deltaMap;
    createMock(deltaMap);

    SEQAN_ASSERT_EQ(*find(deltaMap,  2, DeltaTypeSnp()), deltaMap._entries[3]);
    SEQAN_ASSERT_EQ(*find(deltaMap,  1, DeltaTypeSV()), deltaMap._entries[2]);
    SEQAN_ASSERT_EQ(*find(deltaMap,  20, DeltaTypeDel()), deltaMap._entries[8]);
    SEQAN_ASSERT(find(deltaMap,  1, DeltaTypeIns()) == end(deltaMap, Standard()));
    SEQAN_ASSERT(find(deltaMap,  6, DeltaTypeSnp()) == end(deltaMap, Standard()));
}

SEQAN_DEFINE_TEST(test_delta_map_size)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(size(deltaMap), 0u);
    createMock(deltaMap);
    SEQAN_ASSERT_EQ(size(deltaMap), 10u);
}

SEQAN_DEFINE_TEST(test_delta_map_empty)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(empty(deltaMap), true);
    createMock(deltaMap);
    SEQAN_ASSERT_EQ(empty(deltaMap), false);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_copy_constructor)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap);

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
    createMock(deltaMap);

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
    SEQAN_ASSERT_EQ(getDeltaPosition(*it), 0u);

    unsigned counter = 0;
    for (; it != end(deltaMap, Standard()); ++it, ++counter)
        SEQAN_ASSERT_EQ(*(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, size(deltaMap));

    for (; it != begin(deltaMap, Standard()); --it, --counter)
        SEQAN_ASSERT_EQ(*(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

    for (; !(it == end(deltaMap, Standard())); it++, ++counter)
        SEQAN_ASSERT_EQ(*(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, size(deltaMap));

    for (; !(it == begin(deltaMap, Standard())); it--, --counter)
        SEQAN_ASSERT_EQ(*(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

}

SEQAN_DEFINE_TEST(test_delta_map_iterator)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;

    TDeltaMap deltaMap;
    createMock(deltaMap);
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
    createMock(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(*(it), deltaMap._entries[0]);
    SEQAN_ASSERT_EQ(*(++it), deltaMap._entries[1]);
    SEQAN_ASSERT_EQ(*(++it), deltaMap._entries[2]);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(*(it2++), deltaMap._entries[0]);
    SEQAN_ASSERT_EQ(*(it2++), deltaMap._entries[1]);
    SEQAN_ASSERT_EQ(*(it2), deltaMap._entries[2]);
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
    createMock(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeSnp()), TSnp('A'));
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeDel()), TDel(3));
    SEQAN_ASSERT_EQ(deltaValue(it++, DeltaTypeSV()), TSV(2, "TGAT"));
    it+=3;
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeIns()), TIns("ACGT"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeSnp()), TSnp('A'));
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeDel()), TDel(3));
    SEQAN_ASSERT_EQ(deltaValue(it2++, DeltaTypeSV()), TSV(2, "TGAT"));
    it2+=3;
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeIns()), TIns("ACGT"));
}

SEQAN_DEFINE_TEST(test_delta_map_entry_delta_type)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(getDeltaType(*it),   DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaType(*++it), DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(getDeltaType(*++it), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(getDeltaType(*++it), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaType(*++it), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(getDeltaType(*(it+2)), DELTA_TYPE_INS);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(getDeltaType(*it2++), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaType(*it2++), DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(getDeltaType(*it2++), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(getDeltaType(*it2++), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaType(*it2++), DELTA_TYPE_SV);
}

SEQAN_DEFINE_TEST(test_delta_map_entry_delta_position)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(getDeltaPosition(*it), 0u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 1u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 1u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 2u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 2u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 3u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*++it), 4u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 0u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 1u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 1u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 2u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 2u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 3u);
    SEQAN_ASSERT_EQ(getDeltaPosition(*it2++), 4u);
}

SEQAN_DEFINE_TEST(test_delta_map_entry_delta_coverage)
{
    typedef DeltaMap<TestDeltaMapConfig, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaCoverage<TDeltaMap>::Type TCoverage;

    TDeltaMap deltaMap;
    createMock(deltaMap);

    TCoverage cov1;
    _getCoverage(cov1, 10, 3);
    TCoverage cov2;
    _getCoverage(cov2, 10, 5);
    TCoverage cov3;
    _getCoverage(cov3, 10, 2);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it++), cov2);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it++), cov2);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*(it+2)), cov3);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it2++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it2++), cov2);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it2++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it2++), cov2);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*it2++), cov1);
    SEQAN_ASSERT_EQ(getDeltaCoverage(*(it2+2)), cov3);
}

#endif  // TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
