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
// Mock generator for the journaled string tree.
// ==========================================================================

#ifndef TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_MOCK_H_
#define TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_MOCK_H_

#include <vector>
#include <seqan/sequence.h>

using namespace seqan;

struct TestJstPosConfig_
{
    typedef uint16_t TDeltaPos;
    typedef Dna      TSnpValue;
    typedef uint16_t TDelValue;
    typedef String<Dna> TInsValue;
    typedef Pair<TDelValue, TInsValue> TSVValue;
};

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

class JstMockGenerator
{
public:
    template <typename TJst>
    static TJst _createSimpleJst()
    {
        typename Host<typename Member<TJst, JstSourceMember>::Type>::Type const seq = "AGATCGAGCGAGCTAGCGACTCAG";
        TJst jst(seq, 100);
        String<unsigned> ids;
        appendValue(ids, 0);
        appendValue(ids, 3);
        appendValue(ids, 9);
        appendValue(ids, 99);

        insert(jst, 1, 3, ids, DeltaTypeDel());
        insert(jst, 8, "CGTA", ids, DeltaTypeIns());
        insert(jst, 10, 'C', ids, DeltaTypeSnp());
        insert(jst, 15, 2, ids, DeltaTypeDel());
        insert(jst, 20, 'A', ids, DeltaTypeSnp());

        return jst;
    }

    template <typename TJst>
    static TJst _createComplexJst()
    {
        //    POS	TYPE	S0	S1	S2	S3	S4	S5	S6	S7	S8	S9	MOD
        //     0	SNP		0	1	0	1	0	1	1	1	0	0	C
        //     0	INS		1	0	1	0	0	0	0	0	0	0	TTG
        //     0	DEL		0	0	0	0	1	0	0	0	0	1	2
        //     1	SNP		0	1	1	1	0	0	0	1	0	0	G
        //     5	STV		1	0	0	0	1	0	0	1	0	0	3,CTC
        //     6	INS		0	1	1	0	0	0	1	0	0	0	GGTAGACAACG
        //     8	DEL		1	0	0	0	0	0	1	1	0	1	10
        //     8	DEL		0	1	0	0	1	0	0	0	0	0	5
        //    10	SNP		0	0	1	0	0	1	0	0	0	0	G
        //    13	DEL		0	1	0	0	0	0	0	0	0	0	5
        //    13	DEL		0	0	0	0	1	1	0	0	0	0	8
        //    20	SNP		1	1	0	0	0	0	1	0	0	0	T
        //    40	INS		1	0	1	0	0	0	1	0	1	0	TGAC
        //    41	INS		1	0	1	0	0	0	0	1	0	1	GTAGA
        //    41	STV		0	1	0	0	0	1	1	0	0	0	5,A
        //    58	DEL		1	0	0	0	0	0	0	0	1	1	10
        //    60	SNP		0	1	1	1	1	1	1	1	0	0	G
        //    61	STV		0	0	1	0	0	1	1	1	0	0	2,TA
        //    63	DEL		0	0	0	1	0	1	0	1	0	0	3
        //    66	INS		0	0	0	0	1	0	1	1	0	0	AGGTA
        //    90	DEL		1	0	0	0	0	0	0	0	0	0	9
        //    90	DEL		0	1	0	0	0	0	0	0	0	0	10
        //    99	SNP		0	0	1	0	0	0	0	0	0	0	T
        //    100	INS		0	0	0	1	0	0	0	0	0	0	CCT

        typename Host<typename Member<TJst, JstSourceMember>::Type>::Type const seq = "ACGTGGGATTACGGACGAGCGGACTGATTTAGGAGGACAGTACGGGACTACGGACGTACTACGAGCGACTATCGACGAGAGGAGATCAGGCATTAAGCCC";
          // 01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
          // 0         1         2         3         4         5         6         7         8         9         10
        TJst jst(std::move(seq), 10);

        insert(jst, 0, 'C',                                 std::vector<unsigned>{1, 3, 5, 6, 7},       DeltaTypeSnp());
        insert(jst, 0, "TTG",                               std::vector<unsigned>{0, 2},                DeltaTypeIns());
        insert(jst, 0, 2,                                   std::vector<unsigned>{4, 9},                DeltaTypeDel());
        insert(jst, 1, 'G',                                 std::vector<unsigned>{1, 2, 3, 7},          DeltaTypeSnp());
        insert(jst, 5, Pair<unsigned, DnaString>(3, "CTC"), std::vector<unsigned>{0, 4, 7},             DeltaTypeSV());
        insert(jst, 6, "GGTAGACAACG",                       std::vector<unsigned>{1, 2, 6},             DeltaTypeIns());
        insert(jst, 8, 10,                                  std::vector<unsigned>{0, 6, 7, 9},          DeltaTypeDel());
        insert(jst, 8, 5,                                   std::vector<unsigned>{1, 4},                DeltaTypeDel());
        insert(jst, 10, 'G',                                std::vector<unsigned>{2, 5},                DeltaTypeSnp());
        insert(jst, 13, 5,                                  std::vector<unsigned>{1},                   DeltaTypeDel());
        insert(jst, 13, 8,                                  std::vector<unsigned>{4, 5},                DeltaTypeDel());
        insert(jst, 20, 'T',                                std::vector<unsigned>{0, 1, 6},             DeltaTypeSnp());
        insert(jst, 40, "TGAC",                             std::vector<unsigned>{0, 2, 6, 8},          DeltaTypeIns());
        insert(jst, 41, "GTAGA",                            std::vector<unsigned>{0, 2, 7, 9},          DeltaTypeIns());
        insert(jst, 41, Pair<unsigned, DnaString>(5, "A"),  std::vector<unsigned>{1, 5, 6},             DeltaTypeSV());
        insert(jst, 58, 10,                                 std::vector<unsigned>{0, 8, 9},             DeltaTypeDel());
        insert(jst, 60, 'G',                                std::vector<unsigned>{1, 2, 3, 4, 5, 6, 7}, DeltaTypeSnp());
        insert(jst, 61, Pair<unsigned, DnaString>(2, "TA"), std::vector<unsigned>{2, 5, 6, 7},          DeltaTypeSV());
        insert(jst, 63, 3,                                  std::vector<unsigned>{3, 5, 7},             DeltaTypeDel());
        insert(jst, 66, "AGGTA",                            std::vector<unsigned>{4, 6, 7},             DeltaTypeIns());
        insert(jst, 90, 9,                                  std::vector<unsigned>{0},                   DeltaTypeDel());
        insert(jst, 90, 10,                                 std::vector<unsigned>{1},                   DeltaTypeDel());
        insert(jst, 99, 'T',                                std::vector<unsigned>{2},                   DeltaTypeSnp());
        insert(jst, 100, "CCT",                             std::vector<unsigned>{3},                   DeltaTypeIns());
        
        return jst;
    }

    static StringSet<DnaString> _createComplexTestSet()
    {
        StringSet<DnaString> set;
        appendValue(set, "TTGACGTGCTCGCTGACTGATTTAGGAGGACAGTGACTGTAGAACGGGACTACGGACGTACTATCGACGAGAGGAGATCAGGC");
        appendValue(set, "CGGTGGGGTAGACAACGGAGCTGACTGATTTAGGAGGACAGTAACTACGGACGTACTGCGAGCGACTATCGACGAGAGGAGATCAGG");
        appendValue(set, "TTGAGGTGGGGTAGACAACGGATTGCGGACGAGCGGACTGATTTAGGAGGACAGTGACTGTAGAACGGGACTACGGACGTACTGTAAGCGACTATCGACGAGAGGAGATCAGGCATTAAGCCT");
        appendValue(set, "CGGTGGGATTACGGACGAGCGGACTGATTTAGGAGGACAGTACGGGACTACGGACGTACTGCGGACTATCGACGAGAGGAGATCAGGCATTAAGCCCCCT");
        appendValue(set, "GTGCTCGACTGATTTAGGAGGACAGTACGGGACTACGGACGTACTGCGAGCAGGTAGACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        appendValue(set, "CCGTGGGATTGCGGACTGATTTAGGAGGACAGTAACTACGGACGTACTGTAGACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        appendValue(set, "CCGTGGGGTAGACAACGGAGCTGACTGATTTAGGAGGACAGTGACTAACTACGGACGTACTGTAAGCAGGTAGACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        appendValue(set, "CGGTGCTCGCGGACTGATTTAGGAGGACAGTGTAGAACGGGACTACGGACGTACTGTAAGGTAGACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        appendValue(set, "ACGTGGGATTACGGACGAGCGGACTGATTTAGGAGGACAGTGACTACGGGACTACGGACGTACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        appendValue(set, "GTGGGAGCGGACTGATTTAGGAGGACAGTGTAGAACGGGACTACGGACGTACTATCGACGAGAGGAGATCAGGCATTAAGCCC");
        return set;
    }

    static StringSet<DnaString> generateNeedles()
    {
        std::vector<unsigned> ndlLengths = {5, 12, 23, 32, 66};
        StringSet<DnaString> ndlSet;

        auto set = _createComplexTestSet();
        for (unsigned id = 0; id < length(set); ++id)
        {
            for (unsigned pos = 0; pos < length(set[id]) - 5; ++pos)
            {
                auto ndlSize = ndlLengths[pos % length(ndlLengths)];
                ndlSize = _min(ndlSize, length(set[id]) - pos);
                appendValue(ndlSet, infix(set[id], pos, pos + ndlSize));
            }
        }
        return ndlSet;
    }
};

#endif // TESTS_JOURNALED_STRING_TREE_TEST_JOURNALED_STRING_TREE_MOCK_H_
