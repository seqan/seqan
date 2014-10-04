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
// Test for the delta map.
// ==========================================================================

#ifndef EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
#define EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

template <typename TDeltaMap>
inline void
_generateSnp(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TSnp snp0 = 'C';
    TSnp snp1 = 'A';

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    cov0[0] = true;
    insert(deltaMap, 10, snp1, cov0, DeltaTypeSnp());
    cov1[1] = true;
    insert(deltaMap, 20, snp0, cov1, DeltaTypeSnp());
}

template <typename TDeltaMap>
inline void
_generateDel(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TDel del0 = 2;
    TDel del1 = 5;
    cov0[2] = true;
    insert(deltaMap, 15, del0, cov0, DeltaTypeDel());
    cov1[3] = true;
    insert(deltaMap, 20, del1, cov1, DeltaTypeDel());
}

template <typename TDeltaMap>
inline void
_generateIns(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TIns ins0 = "ACGT";
    TIns ins1 = "C";
    cov0[4] = true;
    TIns ins = "ACGT";
    insert(deltaMap, 21, ins0, cov0, DeltaTypeIns());
    cov1[5] = true;
    insert(deltaMap, 0, ins1, cov1, DeltaTypeIns());
}

template <typename TDeltaMap>
inline void
_generateIndel(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TInDel;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TInDel indel0(2, "AAAA");
    TInDel indel1(3, "GTA");
    cov0[6] = true;
    insert(deltaMap, 12, indel0, cov0, DeltaTypeSV());
    cov1[7] = true;
    insert(deltaMap, 19, indel1, cov1, DeltaTypeSV());
}

template <typename TDeltaMap>
inline void
_testDetlaMapfill(TDeltaMap & deltaMap)
{
    _generateSnp(deltaMap);
    _generateDel(deltaMap);
    _generateIns(deltaMap);
    _generateIndel(deltaMap);
}

SEQAN_DEFINE_TEST(test_delta_map_insert)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaTypeSV>::Type TInDel;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[0]), 0u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[1]), 10u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[2]), 12u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[3]), 15u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[4]), 19u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[5]), 20u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[6]), 20u);
    SEQAN_ASSERT_EQ(getDeltaPosition(deltaMap._entries[7]), 21u);

    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[0])[5], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[1])[0], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[2])[6], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[3])[2], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[4])[7], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[5])[3], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[6])[1], true);
    SEQAN_ASSERT_EQ(getDeltaCoverage(deltaMap._entries[7])[4], true);

    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[0]).i2, 1u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[1]).i2, 0u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[2]).i2, 0u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[3]).i2, 0u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[4]).i2, 1u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[5]).i2, 1u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[6]).i2, 1u);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[7]).i2, 0u);

    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[0]).i1, DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[1]).i1, DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[2]).i1, DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[3]).i1, DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[4]).i1, DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[5]).i1, DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[6]).i1, DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(getDeltaRecord(deltaMap._entries[7]).i1, DELTA_TYPE_INS);

    SEQAN_ASSERT_EQ(deltaMap._deltaStore._insData[1], TIns("C"));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._snpData[0], TSnp('A'));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._svData[0], TInDel(2, "AAAA"));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._delData[0], TDel(2));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._svData[1], TInDel(3, "GTA"));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._delData[1], TDel(5));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._snpData[1], TSnp('C'));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._insData[0], TIns("ACGT"));
}

SEQAN_DEFINE_TEST(test_delta_map_length)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(length(deltaMap), 0u);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(length(deltaMap), 2u);
    _generateDel(deltaMap);
    SEQAN_ASSERT_EQ(length(deltaMap), 4u);
    _generateIns(deltaMap);
    SEQAN_ASSERT_EQ(length(deltaMap), 6u);
    _generateIndel(deltaMap);
    SEQAN_ASSERT_EQ(length(deltaMap), 8u);
}

SEQAN_DEFINE_TEST(test_delta_map_empty)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(empty(deltaMap), true);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(empty(deltaMap), false);
}

SEQAN_DEFINE_TEST(test_delta_map_coverage_size)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 0u);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 8u);
    _generateDel(deltaMap);
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 8u);
}

SEQAN_DEFINE_TEST(test_delta_map_set_coverage_size)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 0u);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 8u);

    setCoverageSize(deltaMap, 1u);
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 1u);
    SEQAN_ASSERT_EQ(length(getDeltaCoverage(deltaMap._entries[0])), 1u);
    SEQAN_ASSERT_EQ(length(getDeltaCoverage(deltaMap._entries[1])), 1u);

    setCoverageSize(deltaMap, 10u);
    SEQAN_ASSERT_EQ(getCoverageSize(deltaMap), 10u);
    SEQAN_ASSERT_EQ(length(getDeltaCoverage(deltaMap._entries[0])), 10u);
    SEQAN_ASSERT_EQ(length(getDeltaCoverage(deltaMap._entries[1])), 10u);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_copy_constructor)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

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
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

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

template <typename TContainer>
inline void
_testDeltaMapIterator(TContainer & deltaMap)
{
    typedef typename Iterator<TContainer, Standard>::Type TIterator;

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it), 0u);

    unsigned counter = 0;
    for (; it != end(deltaMap, Standard()); ++it, ++counter)
        SEQAN_ASSERT_EQ(value(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; it != begin(deltaMap, Standard()); --it, --counter)
        SEQAN_ASSERT_EQ(value(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

    for (; !(it == end(deltaMap, Standard())); it++, ++counter)
        SEQAN_ASSERT_EQ(value(it), deltaMap._entries[counter]);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; !(it == begin(deltaMap, Standard())); it--, --counter)
        SEQAN_ASSERT_EQ(value(it - 1), deltaMap._entries[counter - 1]);
    SEQAN_ASSERT_EQ(counter, 0u);

}

SEQAN_DEFINE_TEST(test_delta_map_iterator)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);
    const TDeltaMap deltaMap2 = deltaMap;

    _testDeltaMapIterator(deltaMap);
    _testDeltaMapIterator(deltaMap2);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_value)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

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
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaType(it), DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(deltaType(++it), DELTA_TYPE_DEL);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_SV);
    SEQAN_ASSERT_EQ(deltaType(it2++), DELTA_TYPE_DEL);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_position)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 10u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 12u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 15u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 10u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 12u);
    SEQAN_ASSERT_EQ(deltaPosition(it2), 15u);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_value)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaTypeSnp>::Type TSnp;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeSnp()), TSnp('A'));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeSnp()), TSnp('A'));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_ins)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaTypeIns>::Type TIns;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeIns()), TIns("C"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeIns()), TIns("C"));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_del)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaTypeDel>::Type TDel;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;++it;++it;
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeDel()), TDel(2));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;++it2;++it2;
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeDel()), TDel(2));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_sv)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaTypeSV>::Type TSV;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;++it;
    SEQAN_ASSERT_EQ(deltaValue(it, DeltaTypeSV()), TSV(2, "AAAA"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;++it2;
    SEQAN_ASSERT_EQ(deltaValue(it2, DeltaTypeSV()), TSV(2, "AAAA"));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_coverage)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaCoverage<TDeltaMap>::Type TDeltaCoverage;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TDeltaCoverage bitVec;
    resize(bitVec, 8, false);
    bitVec[5] = true;
    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaCoverage(it), bitVec);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaCoverage(it2), bitVec);
}


#endif  // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_DELTA_MAP_H_
