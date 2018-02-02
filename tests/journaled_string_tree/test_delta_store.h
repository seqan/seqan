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
// Tests for the delta store.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_DELTA_STORE_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_DELTA_STORE_H_

#include <seqan/basic.h>
#include <seqan/journaled_string_tree.h>

using namespace seqan;

template <typename TSnp, typename TDel>
impl::DeltaStore<TSnp, TDel> createMock()
{
    impl::DeltaStore<TSnp, TDel> store;

    resize(store._snpData, 3);
    store._snpData[0] = 'A';
    store._snpData[1] = 'G';
    store._snpData[2] = 'C';

    resize(store._insData, 2);
    store._insData[0] = "AAAAT";
    store._insData[1] = "CG";

    resize(store._delData, 1);
    store._delData[0] = 3;

    resize(store._svData, 2);
    store._svData[0].i1 = 3; store._svData[0].i2 = "CGTAA";
    store._svData[1].i1 = 4; store._svData[1].i2 = "AA";

    return store;
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_is_delta_type)
{
    SEQAN_ASSERT(isDeltaType(DELTA_TYPE_SNP, DeltaTypeSnp()));
    SEQAN_ASSERT_NOT(isDeltaType(DELTA_TYPE_SNP, DeltaTypeDel()));
    SEQAN_ASSERT_NOT(isDeltaType(DELTA_TYPE_SNP, DeltaTypeIns()));
    SEQAN_ASSERT_NOT(isDeltaType(DELTA_TYPE_SNP, DeltaTypeSV()));
    SEQAN_ASSERT(isDeltaType(DELTA_TYPE_DEL, DeltaTypeDel()));
    SEQAN_ASSERT(isDeltaType(DELTA_TYPE_INS, DeltaTypeIns()));
    SEQAN_ASSERT(isDeltaType(DELTA_TYPE_SV, DeltaTypeSV()));
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_select_delta_type)
{
    SEQAN_ASSERT_EQ(selectDeltaType(DeltaTypeSnp()), DELTA_TYPE_SNP);
    SEQAN_ASSERT_EQ(selectDeltaType(DeltaTypeDel()), DELTA_TYPE_DEL);
    SEQAN_ASSERT_EQ(selectDeltaType(DeltaTypeIns()), DELTA_TYPE_INS);
    SEQAN_ASSERT_EQ(selectDeltaType(DeltaTypeSV()), DELTA_TYPE_SV);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_get_delta_store)
{
    impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();

    SEQAN_ASSERT(getDeltaStore(store, DeltaTypeSnp()) == store._snpData);
    SEQAN_ASSERT(getDeltaStore(store, DeltaTypeDel()) == store._delData);
    SEQAN_ASSERT(getDeltaStore(store, DeltaTypeIns()) == store._insData);
    SEQAN_ASSERT(getDeltaStore(store, DeltaTypeSV())  == store._svData);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_add_delta_value)
{
    typedef Pair<unsigned, DnaString> TPair;
    { // Add into empty store.
        impl::DeltaStore<Dna, unsigned> store;
        SEQAN_ASSERT_EQ(addDeltaValue(store, 'C', DeltaTypeSnp()), 0u);
        SEQAN_ASSERT_EQ(store._snpData[0], 'C');
        SEQAN_ASSERT_EQ(addDeltaValue(store, "ATG", DeltaTypeIns()), 0u);
        SEQAN_ASSERT_EQ(store._insData[0], "ATG");
        SEQAN_ASSERT_EQ(addDeltaValue(store, 2, DeltaTypeDel()), 0u);
        SEQAN_ASSERT_EQ(store._delData[0], 2u);
        SEQAN_ASSERT_EQ(addDeltaValue(store, TPair(3, "CGT"), DeltaTypeSV()), 0u);
        SEQAN_ASSERT_EQ(store._svData[0].i1, 3u);
        SEQAN_ASSERT_EQ(store._svData[0].i2, "CGT");
    }

    { // Add into filled store.
        impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();
        SEQAN_ASSERT_EQ(addDeltaValue(store, 'C', DeltaTypeSnp()), 3u);
        SEQAN_ASSERT_EQ(store._snpData[3], 'C');
        SEQAN_ASSERT_EQ(addDeltaValue(store, "CG", DeltaTypeIns()), 2u);
        SEQAN_ASSERT_EQ(store._insData[2], "CG");
        SEQAN_ASSERT_EQ(addDeltaValue(store, 2, DeltaTypeDel()), 1u);
        SEQAN_ASSERT_EQ(store._delData[1], 2u);
        SEQAN_ASSERT_EQ(addDeltaValue(store, TPair(3, "CGT"), DeltaTypeSV()), 2u);
        SEQAN_ASSERT_EQ(store._svData[2].i1, 3u);
        SEQAN_ASSERT_EQ(store._svData[2].i2, "CGT");
    }
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_erase_delta_value)
{   
    { // Erase from empty store.
        impl::DeltaStore<Dna, unsigned> store;
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 0, DeltaTypeSnp()), 0u);
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 3, DeltaTypeDel()), 0u);
    }

    { // Erase from filled store.
        impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 1, DeltaTypeSnp()), 2u);
        SEQAN_ASSERT_EQ(store._snpData[1], 'C');
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 0, DeltaTypeIns()), 1u);
        SEQAN_ASSERT_EQ(store._insData[0], "CG");
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 0, DeltaTypeDel()), 0u);
        SEQAN_ASSERT(empty(store._delData));
        SEQAN_ASSERT_EQ(eraseDeltaValue(store, 1, DeltaTypeSV()), 1u);
        SEQAN_ASSERT_EQ(store._svData[0].i1, 3u);
        SEQAN_ASSERT_EQ(store._svData[0].i2, "CGTAA");
    }
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_delta_value)
{
    typedef Pair<unsigned, DnaString> TPair;

    impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();

    SEQAN_ASSERT_EQ(deltaValue(store, 2, DeltaTypeSnp()), 'C');
    SEQAN_ASSERT_EQ(deltaValue(store, 1, DeltaTypeIns()), "CG");
    SEQAN_ASSERT_EQ(deltaValue(store, 0, DeltaTypeDel()), 3u);
    SEQAN_ASSERT_EQ(deltaValue(store, 1, DeltaTypeSV()), TPair(4, "AA"));
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_clear)
{
    { // Clear empty store.
        impl::DeltaStore<Dna, unsigned> store;
        clear(store);
        SEQAN_ASSERT_EQ(length(store._snpData), 0u);
        SEQAN_ASSERT_EQ(length(store._insData), 0u);
        SEQAN_ASSERT_EQ(length(store._delData), 0u);
        SEQAN_ASSERT_EQ(length(store._svData), 0u);
    }

    { // Clear filled store.
        impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();
        clear(store);
        SEQAN_ASSERT_EQ(length(store._snpData), 0u);
        SEQAN_ASSERT_EQ(length(store._insData), 0u);
        SEQAN_ASSERT_EQ(length(store._delData), 0u);
        SEQAN_ASSERT_EQ(length(store._svData), 0u);
    }
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_deletion_size)
{
    impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();

    SEQAN_ASSERT_EQ(deletionSize(store, 2, DeltaTypeSnp()), 1u);
    SEQAN_ASSERT_EQ(deletionSize(store, 1, DeltaTypeIns()), 0u);
    SEQAN_ASSERT_EQ(deletionSize(store, 0, DeltaTypeDel()), 3u);
    SEQAN_ASSERT_EQ(deletionSize(store, 1, DeltaTypeSV()), 4u);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_insertion_size)
{
    impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();

    SEQAN_ASSERT_EQ(insertionSize(store, 2, DeltaTypeSnp()), 1u);
    SEQAN_ASSERT_EQ(insertionSize(store, 1, DeltaTypeIns()), 2u);
    SEQAN_ASSERT_EQ(insertionSize(store, 0, DeltaTypeDel()), 0u);
    SEQAN_ASSERT_EQ(insertionSize(store, 1, DeltaTypeSV()), 2u);
}

SEQAN_DEFINE_TEST(test_delta_map_delta_store_net_size)
{
    impl::DeltaStore<Dna, unsigned> store = createMock<Dna, unsigned>();

    SEQAN_ASSERT_EQ(netSize(store, 2, DeltaTypeSnp()), 0);
    SEQAN_ASSERT_EQ(netSize(store, 1, DeltaTypeIns()), 2);
    SEQAN_ASSERT_EQ(netSize(store, 0, DeltaTypeDel()), -3);
    SEQAN_ASSERT_EQ(netSize(store, 1, DeltaTypeSV()), -2);
}

#endif // TESTS_JOURNALED_STRING_TREE_TEST_DELTA_STORE_H_
