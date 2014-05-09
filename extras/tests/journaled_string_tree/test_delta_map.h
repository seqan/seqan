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

template <typename TDeltaMap, typename TExpectedType>
void _testDeltaMapValueMF(TExpectedType &)
{
    typedef typename Value<TDeltaMap>::Type TValue;
    bool res = IsSameType<TValue, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

template <typename TDeltaMap, typename TExpectedType>
void _testDeltaMapReferenceMF(TExpectedType &)
{
    typedef typename Reference<TDeltaMap>::Type TReference;
    typedef TExpectedType & TestType;
    bool res = IsSameType<TReference, TestType>::VALUE;
    SEQAN_ASSERT(res);
}

template <typename TDeltaMap, typename DeltaType::TValue TYPE, typename TExpectedType>
void _testDeltaMapDeltaValueMF(TExpectedType &)
{
    typedef typename DeltaValue<TDeltaMap, TYPE>::Type TDeltaValue;
    bool res = IsSameType<TDeltaValue, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

template <typename TDeltaMap, typename TExpectedType>
void _testDeltaMapDeltaCoverageMF(TExpectedType &)
{
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;
    bool res = IsSameType<TCoverage, TExpectedType>::VALUE;
    SEQAN_ASSERT(res);
}

SEQAN_DEFINE_TEST(test_delta_map_value_mf)
{
    short var1 = 0;
    _testDeltaMapValueMF<DeltaMap<short, char, Default> >(var1);
    const int var2 = 0;
    _testDeltaMapValueMF<const DeltaMap<int, Dna, Default> >(var2);
    char var3 = 'a';
    _testDeltaMapValueMF<DeltaMap<char, char, Default> >(var3);
    const Dna5 var4 = 'A';
    _testDeltaMapValueMF<const DeltaMap<Dna5, AminoAcid, Default> >(var4);
}

SEQAN_DEFINE_TEST(test_delta_map_reference_mf)
{
    short var1 = 0;
    _testDeltaMapReferenceMF<DeltaMap<short, char, Default> >(var1);
    const int var2 = 0;
    _testDeltaMapReferenceMF<const DeltaMap<int, Dna, Default> >(var2);
    char var3 = 'a';
    _testDeltaMapReferenceMF<DeltaMap<char, char, Default> >(var3);
    const Dna5 var4 = 'A';
    _testDeltaMapReferenceMF<const DeltaMap<Dna5, AminoAcid, Default> >(var4);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_value_mf)
{
    char snp1;
    _testDeltaMapDeltaValueMF<DeltaMap<short, char, Default>, DeltaType::DELTA_TYPE_SNP>(snp1);
    const Dna snp2;
    _testDeltaMapDeltaValueMF<const DeltaMap<short, Dna, Default>, DeltaType::DELTA_TYPE_SNP>(snp2);

    String<char> ins1;
    _testDeltaMapDeltaValueMF<DeltaMap<short, char, Default>, DeltaType::DELTA_TYPE_INS>(ins1);
    const String<Dna> ins2;
    _testDeltaMapDeltaValueMF<const DeltaMap<short, Dna, Default>, DeltaType::DELTA_TYPE_INS>(ins2);

    typedef Size<DeltaMap<short, char, Default> >::Type TSize_;
    TSize_ del1;
    _testDeltaMapDeltaValueMF<DeltaMap<short, char, Default>, DeltaType::DELTA_TYPE_DEL>(del1);
    const TSize_ del2 = 0;
    _testDeltaMapDeltaValueMF<const DeltaMap<short, Dna, Default>, DeltaType::DELTA_TYPE_DEL>(del2);

    Pair<TSize_, String<char> > indel1;
    _testDeltaMapDeltaValueMF<DeltaMap<short, char, Default>, DeltaType::DELTA_TYPE_INDEL>(indel1);
    const Pair<TSize_, String<Dna> > indel2;
    _testDeltaMapDeltaValueMF<const DeltaMap<short, Dna, Default>, DeltaType::DELTA_TYPE_INDEL>(indel2);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_coverage_mf)
{
    typedef String<bool, Packed<> > TBitVec;
    TBitVec bitVec;
    _testDeltaMapDeltaCoverageMF<DeltaMap<short, char, Default> >(bitVec);
    const TBitVec bitVecConst;
    _testDeltaMapDeltaCoverageMF<const DeltaMap<short, Dna, Default> >(bitVecConst);
}

template <typename TDeltaMap>
inline void
_generateSnp(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TSnp;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TSnp snp0 = 'C';
    TSnp snp1 = 'A';

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    cov0[0] = true;
    insert(deltaMap, 10, snp1, cov0);
    cov1[1] = true;
    insert(deltaMap, 20, snp0, cov1);
}

template <typename TDeltaMap>
inline void
_generateDel(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TDel del0 = 2;
    TDel del1 = 5;
    cov0[2] = true;
    insert(deltaMap, 15, del0, cov0);
    cov1[3] = true;
    insert(deltaMap, 20, del1, cov1);
}

template <typename TDeltaMap>
inline void
_generateIns(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TIns;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TIns ins0 = "ACGT";
    TIns ins1 = "C";
    cov0[4] = true;
    TIns ins = "ACGT";
    insert(deltaMap, 21, ins0, cov0);
    cov1[5] = true;
    insert(deltaMap, 0, ins1, cov1);
}

template <typename TDeltaMap>
inline void
_generateIndel(TDeltaMap & deltaMap)
{
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TInDel;
    typedef typename DeltaCoverage<TDeltaMap>::Type TCoverage;

    TCoverage cov0;
    resize(cov0, 8, false);
    TCoverage cov1 = cov0;

    TInDel indel0(2, "AAAA");
    TInDel indel1(3, "GTA");
    cov0[6] = true;
    insert(deltaMap, 12, indel0, cov0);
    cov1[7] = true;
    insert(deltaMap, 19, indel1, cov1);
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
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TSnp;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TIns;
    typedef typename DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TInDel;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    SEQAN_ASSERT_EQ(deltaMap._keys[0], 0u);
    SEQAN_ASSERT_EQ(deltaMap._keys[1], 10u);
    SEQAN_ASSERT_EQ(deltaMap._keys[2], 12u);
    SEQAN_ASSERT_EQ(deltaMap._keys[3], 15u);
    SEQAN_ASSERT_EQ(deltaMap._keys[4], 19u);
    SEQAN_ASSERT_EQ(deltaMap._keys[5], 20u);
    SEQAN_ASSERT_EQ(deltaMap._keys[6], 20u);
    SEQAN_ASSERT_EQ(deltaMap._keys[7], 21u);

    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[0][5], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[1][0], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[2][6], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[3][2], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[4][7], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[5][3], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[6][1], true);
    SEQAN_ASSERT_EQ(deltaMap._deltaCoverageStore._coverageData[7][4], true);

    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[0]), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[1]), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[2]), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[3]), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[4]), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[5]), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[6]), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(deltaMap._deltaStore._varDataMap[7]), 0u);

    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[0]), DeltaType::DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[1]), DeltaType::DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[2]), DeltaType::DELTA_TYPE_INDEL);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[3]), DeltaType::DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[4]), DeltaType::DELTA_TYPE_INDEL);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[5]), DeltaType::DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[6]), DeltaType::DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(deltaMap._deltaStore._varDataMap[7]), DeltaType::DELTA_TYPE_INS);

    SEQAN_ASSERT_EQ(deltaMap._deltaStore._insData[1], TIns("C"));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._snpData[0], TSnp('A'));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._indelData[0], TInDel(2, "AAAA"));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._delData[0], TDel(2));
    SEQAN_ASSERT_EQ(deltaMap._deltaStore._indelData[1], TInDel(3, "GTA"));
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
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 0u);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 8u);
    _generateDel(deltaMap);
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 8u);
}

SEQAN_DEFINE_TEST(test_delta_map_set_coverage_size)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;

    TDeltaMap deltaMap;
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 0u);
    _generateSnp(deltaMap);
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 8u);

    setCoverageSize(deltaMap, 1u);
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 1u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaCoverageStore._coverageData[0]), 1u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaCoverageStore._coverageData[1]), 1u);

    setCoverageSize(deltaMap, 10u);
    SEQAN_ASSERT_EQ(coverageSize(deltaMap), 10u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaCoverageStore._coverageData[0]), 10u);
    SEQAN_ASSERT_EQ(length(deltaMap._deltaCoverageStore._coverageData[1]), 10u);
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

    SEQAN_ASSERT_EQ(*copyItMap._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*copyItMap._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*copyItMap._deltaCoverageIter, *itMap._deltaCoverageIter);

    SEQAN_ASSERT_EQ(*constCopyItMap._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*constCopyItMap._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*constCopyItMap._deltaCoverageIter, *itMap._deltaCoverageIter);

    SEQAN_ASSERT_EQ(*constCopyItMap2._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*constCopyItMap2._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*constCopyItMap2._deltaCoverageIter, *itMap._deltaCoverageIter);
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

    SEQAN_ASSERT_EQ(*copyItMap._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*copyItMap._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*copyItMap._deltaCoverageIter, *itMap._deltaCoverageIter);

    SEQAN_ASSERT_EQ(*constCopyItMap._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*constCopyItMap._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*constCopyItMap._deltaCoverageIter, *itMap._deltaCoverageIter);

    SEQAN_ASSERT_EQ(*constCopyItMap2._mapIter, *itMap._mapIter);
    SEQAN_ASSERT_EQ(*constCopyItMap2._deltaStoreIter, *itMap._deltaStoreIter);
    SEQAN_ASSERT_EQ(*constCopyItMap2._deltaCoverageIter, *itMap._deltaCoverageIter);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(*it, 0u);

    unsigned counter = 0;
    for (; it != end(deltaMap, Standard()); ++it, ++counter);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; it != begin(deltaMap, Standard()); --it, --counter);
    SEQAN_ASSERT_EQ(counter, 0u);

    for (; !(it == end(deltaMap, Standard())); it++, ++counter);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; !(it == begin(deltaMap, Standard())); it--, --counter);
    SEQAN_ASSERT_EQ(counter, 0u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(*it2, 0u);

    for (; it2 != end(deltaMap2, Standard()); ++it2, ++counter);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; it2 != begin(deltaMap2, Standard()); --it2, --counter);
    SEQAN_ASSERT_EQ(counter, 0u);

    for (; !(it2 == end(deltaMap2, Standard())); it2++, ++counter);
    SEQAN_ASSERT_EQ(counter, length(deltaMap));

    for (; !(it2 == begin(deltaMap2, Standard())); it2--, --counter);
    SEQAN_ASSERT_EQ(counter, 0u);
}


SEQAN_DEFINE_TEST(test_delta_map_iterator_value)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(value(it), 0u);
    SEQAN_ASSERT_EQ(value(++it), 10u);
    SEQAN_ASSERT_EQ(value(++it), 12u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(value(it2++), 0u);
    SEQAN_ASSERT_EQ(value(it2++), 10u);
    SEQAN_ASSERT_EQ(value(it2), 12u);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_type)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaType(it), DeltaType::DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(deltaType(++it), DeltaType::DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(++it), DeltaType::DELTA_TYPE_INDEL);
    SEQAN_ASSERT_EQ(deltaType(++it), DeltaType::DELTA_TYPE_DEL);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaType(it2++), DeltaType::DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(deltaType(it2++), DeltaType::DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(deltaType(it2++), DeltaType::DELTA_TYPE_INDEL);
    SEQAN_ASSERT_EQ(deltaType(it2++), DeltaType::DELTA_TYPE_DEL);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_position)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(++it), 0u);

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 1u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(it2++), 0u);
    SEQAN_ASSERT_EQ(deltaPosition(it2), 0u);
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_snp)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_SNP>::Type TSnp;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;
    SEQAN_ASSERT_EQ(deltaSnp(it), TSnp('A'));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;
    SEQAN_ASSERT_EQ(deltaSnp(it2), TSnp('A'));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_ins)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INS>::Type TIns;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    SEQAN_ASSERT_EQ(deltaIns(it), TIns("C"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    SEQAN_ASSERT_EQ(deltaIns(it2), TIns("C"));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_del)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_DEL>::Type TDel;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;++it;++it;
    SEQAN_ASSERT_EQ(deltaDel(it), TDel(2));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;++it2;++it2;
    SEQAN_ASSERT_EQ(deltaDel(it2), TDel(2));
}

SEQAN_DEFINE_TEST(test_delta_map_iterator_delta_indel)
{
    typedef DeltaMap<unsigned, char, Default> TDeltaMap;
    typedef Iterator<TDeltaMap, Standard>::Type TIterator;
    typedef Iterator<TDeltaMap const, Standard>::Type TConstIterator;
    typedef DeltaValue<TDeltaMap, DeltaType::DELTA_TYPE_INDEL>::Type TInDel;

    TDeltaMap deltaMap;
    _testDetlaMapfill(deltaMap);

    TIterator it = begin(deltaMap, Standard());
    ++it;++it;
    SEQAN_ASSERT_EQ(deltaIndel(it), TInDel(2, "AAAA"));

    const TDeltaMap deltaMap2 = deltaMap;
    TConstIterator it2 = begin(deltaMap2, Standard());
    ++it2;++it2;
    SEQAN_ASSERT_EQ(deltaIndel(it2), TInDel(2, "AAAA"));
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
