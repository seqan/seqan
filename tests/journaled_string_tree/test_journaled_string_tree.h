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

#include "test_config_reader.h"

using namespace seqan;

struct TestDataConfig
{
    
    String<unsigned>            deltaPositions;
    String<DeltaType>           deltaTypes;
    String<String<unsigned> >   ids;
};

SEQAN_DEFINE_TEST(test_journaled_string_tree_constructor)
{
    typedef JournaledStringTree<DnaString> TJst;
    typedef Size<TJst>::Type TSize;
    typedef Host<TJst>::Type TDeltaMap;
    typedef DeltaCoverage<TDeltaMap>::Type TCoverage;

    DnaString hostSeq = "ACGTGGATCTGTACTGACGGACGGACTTGACGGGAGTACGAGCATCGACT";
    DnaString const hostSeqConst = hostSeq;
    {  // 1st Constructor.

        TJst jst1(hostSeq, 10, 100);

        SEQAN_ASSERT_EQ(jst1._dimension, static_cast<TSize>(100));
        SEQAN_ASSERT_EQ(jst1._branchLength, static_cast<TSize>(10));
        SEQAN_ASSERT(!empty(jst1._mapHolder));
        SEQAN_ASSERT(!dependent(jst1._mapHolder));
        SEQAN_ASSERT(!empty(jst1._source));
        SEQAN_ASSERT_EQ(host(jst1._source), hostSeq);

        TJst jst2(hostSeqConst, 10, 100);
        SEQAN_ASSERT_EQ(host(jst2._source), hostSeqConst);
    }

    {  // 2nd Constructor
        TDeltaMap map;
        TCoverage cov;
        resize(cov, 10, Exact());
        insert(map, 2, 'C', cov, DeltaTypeSnp());

        TJst jst(hostSeq, 11, map);
        SEQAN_ASSERT_EQ(jst._dimension, static_cast<TSize>(10));
        SEQAN_ASSERT_EQ(jst._branchLength, static_cast<TSize>(11));
        SEQAN_ASSERT(!empty(jst._mapHolder));
        SEQAN_ASSERT(dependent(jst._mapHolder));
        SEQAN_ASSERT_EQ(host(jst._source), hostSeq);
    }
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_dimension)
{
    typedef Size<JournaledStringTree<Dna5String> >::Type TSize;
    Dna5String seq = "ANGATAGAC";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);
    SEQAN_ASSERT_EQ(dimension(jst), static_cast<TSize>(100));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_insert_node)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    SEQAN_ASSERT(insertNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(insertNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(insertNode(jst, 5, "CGTA", ids, DeltaTypeIns()));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_erase_node)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 9);
    appendValue(ids, 99);

    SEQAN_ASSERT(insertNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(insertNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(insertNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    SEQAN_ASSERT(eraseNode(jst, 5, DeltaTypeIns()));
    SEQAN_ASSERT(eraseNode(jst, 5, DeltaTypeSnp()));
    SEQAN_ASSERT_NOT(eraseNode(jst, 5, DeltaTypeSnp()));
    SEQAN_ASSERT_NOT(eraseNode(jst, 4, DeltaTypeIns()));
    SEQAN_ASSERT(eraseNode(jst, 2, DeltaTypeDel()));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_clear)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);
    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 99);

    SEQAN_ASSERT(insertNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(insertNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(insertNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    clear(jst);
    SEQAN_ASSERT(jst._dimension == 0);
    SEQAN_ASSERT(jst._branchLength == 0);
    SEQAN_ASSERT(empty(jst._source));
    SEQAN_ASSERT(empty(jst._mapHolder));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_container)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);

    SEQAN_ASSERT(empty(host(jst)));

    String<unsigned> ids;
    appendValue(ids, 0);
    appendValue(ids, 3);
    appendValue(ids, 99);

    SEQAN_ASSERT(insertNode(jst, 5, 'C', ids, DeltaTypeSnp()));
    SEQAN_ASSERT(insertNode(jst, 2, 3, ids, DeltaTypeDel()));
    SEQAN_ASSERT(insertNode(jst, 5, "CGTA", ids, DeltaTypeIns()));

    SEQAN_ASSERT(!empty(host(jst)));
}


SEQAN_DEFINE_TEST(test_journaled_string_tree_branch_length)
{
    typedef Size<JournaledStringTree<Dna5String> >::Type TSize;
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);

    SEQAN_ASSERT_EQ(branchLength(jst), static_cast<TSize>(10));
}

SEQAN_DEFINE_TEST(test_journaled_string_tree_set_branch_length)
{
    typedef Size<JournaledStringTree<Dna5String> >::Type TSize;
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String> jst(seq, 10, 100);

    SEQAN_ASSERT_EQ(branchLength(jst), static_cast<TSize>(10));
    setBranchLength(jst, 15);
    SEQAN_ASSERT_EQ(branchLength(jst), static_cast<TSize>(15));
}

struct TestJstPosConfig_
{
    typedef __uint16 TDeltaPos;
    typedef Dna      TSnpValue;
    typedef __uint16 TDelValue;
    typedef String<Dna> TInsValue;
    typedef Pair<TDelValue, TInsValue> TSVValue;
};

SEQAN_DEFINE_TEST(test_journaled_string_tree_max_size)
{
    Dna5String seq = "CGTATAGGANNAGAT";
    JournaledStringTree<Dna5String, TestJstPosConfig_> jst(seq, 10, 100);

    SEQAN_ASSERT_EQ(maxSize(jst), MaxValue<__uint16>::VALUE);
}


template <typename TString, typename TRef, typename TDeltaMap>
inline void
_createJournaledStrings(StringSet<TString> & set,
                        TRef const & source,
                        TDeltaMap const & map,
                        unsigned const dim)
{

    for (unsigned seqId = 0; seqId < dim; ++seqId)
    {
        TString seq = source;
        auto it = end(map, Standard());
        auto itBegin = begin(map, Standard());

        while (it != itBegin)
        {
            --it;
            if (!getDeltaCoverage(*it)[seqId])
                continue;
            switch (getDeltaType(*it))
            {
                case DELTA_TYPE_SNP:
                {
                    erase(seq, getDeltaPosition(*it));
                    insertValue(seq, getDeltaPosition(*it), deltaValue(it, DeltaTypeSnp()));
                    break;
                }
                case DELTA_TYPE_DEL:
                {
                    if((*it).deltaTypeEnd != DeltaEndType::IS_RIGHT)
                        erase(seq, getDeltaPosition(*it), getDeltaPosition(*it)+ deltaValue(it, DeltaTypeDel()));
                    break;
                }
                case DELTA_TYPE_INS:
                {
                    insert(seq, getDeltaPosition(*it), deltaValue(it, DeltaTypeIns()));
                    break;
                }
                case DELTA_TYPE_SV:
                {
                    if((*it).deltaTypeEnd != DeltaEndType::IS_RIGHT)
                    {
                        erase(seq, getDeltaPosition(*it), getDeltaPosition(*it) + deltaValue(it, DeltaTypeSV()).i1);
                        insert(seq, getDeltaPosition(*it), deltaValue(it, DeltaTypeSV()).i2);
                    }
                    break;
                }
            }
        }
        appendValue(set, seq);
    }
}

template <typename TJst, typename TInFile>
inline void
_createJst(TJst & jst, TInFile & configIn)
{
    while (!atEnd(configIn))
    {
        TestConfigRecord_<Dna, unsigned> record;
        readRecord(record, configIn);

        switch (record.deltaType)
        {
            case DELTA_TYPE_SNP:
            {
                SEQAN_ASSERT(insertNode(jst, record.pos, record.snp, record.coverage, DeltaTypeSnp()));
                break;
            }
            case DELTA_TYPE_DEL:
            {
                SEQAN_ASSERT(insertNode(jst, record.pos, record.del, record.coverage, DeltaTypeDel()));
                break;
            }
            case DELTA_TYPE_INS:
            {
                SEQAN_ASSERT(insertNode(jst, record.pos, record.ins, record.coverage, DeltaTypeIns()));
                break;
            }
            case DELTA_TYPE_SV:
            {
                SEQAN_ASSERT(insertNode(jst, record.pos, record.stv, record.coverage, DeltaTypeSV()));
                break;
            }
        }
    }
}

// Define a typed test for the creation of different scenarios.
SEQAN_DEFINE_TEST(test_journaled_string_tree_create)
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/tests/journaled_string_tree/testConfig.txt");

    TestConfigFileIn_ configIn(toCString(path));

    TestConfigHeader_<DnaString> header;
    readHeader(header, configIn);
    DnaString seq = header.ref;

    JournaledStringTree<DnaString, TestJstPosConfig_> jst(seq, 10, length(context(configIn)));
    _createJst(jst, configIn);

    StringSet<DnaString> set;
    _createJournaledStrings(set, seq, host(jst), length(context(configIn)));

    SEQAN_ASSERT(create(jst));
    for (unsigned i = 0; i < length(set); ++i)
        SEQAN_ASSERT(jst._buffer._journaledSet[i] == set[i]);
}

#endif // EXTRAS_TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_H_
