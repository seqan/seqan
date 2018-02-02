// ==========================================================================
//                                  parse_lm
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_PARSE_LM_TEST_PARSE_LM_H_
#define TESTS_PARSE_LM_TEST_PARSE_LM_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/parse_lm.h>

SEQAN_DEFINE_TEST(test_parse_lm_local_match_constructor)
{
    using namespace seqan;

    typedef LocalMatch<unsigned, unsigned> TLocalMatch;

    // Default constructor.
    {
        TLocalMatch localMatch;
        unsigned const maxU = std::numeric_limits<unsigned>::max();

        SEQAN_ASSERT_EQ(maxU, localMatch.id);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectId);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectBeginPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.subjectEndPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryId);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryBeginPos);
        SEQAN_ASSERT_EQ(maxU, localMatch.queryEndPos);
    }
    // Constructor with options.
    {
        TLocalMatch localMatch(0, 1, 2, 3, 4, 5, 6);

        SEQAN_ASSERT_EQ(0u, localMatch.id);
        SEQAN_ASSERT_EQ(1u, localMatch.subjectId);
        SEQAN_ASSERT_EQ(2u, localMatch.subjectBeginPos);
        SEQAN_ASSERT_EQ(3u, localMatch.subjectEndPos);
        SEQAN_ASSERT_EQ(4u, localMatch.queryId);
        SEQAN_ASSERT_EQ(5u, localMatch.queryBeginPos);
        SEQAN_ASSERT_EQ(6u, localMatch.queryEndPos);
    }
}

SEQAN_DEFINE_TEST(test_parse_lm_local_match_store_constructor)
{
    using namespace seqan;

    // Default constructor.
    LocalMatchStore<> store;
}

SEQAN_DEFINE_TEST(test_parse_lm_local_match_store_append_local_match)
{
    using namespace seqan;

    LocalMatchStore<> store;

    // Append with sequence names.
    {
        appendLocalMatch(store, "seq0", 1, 2, "seq1", 3, 4);

        SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
        SEQAN_ASSERT_EQ(store.sequenceNameStore[0], CharString("seq0"));
        SEQAN_ASSERT_EQ(store.sequenceNameStore[1], CharString("seq1"));

        SEQAN_ASSERT_EQ(length(store.matchStore), 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).id, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectId, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectBeginPos, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectEndPos, 2u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryId, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryBeginPos, 3u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryEndPos, 4u);
    }

    // Append with sequence ids.
    {
        appendLocalMatch(store, 0, 5, 6, 1, 7, 8);

        SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
        SEQAN_ASSERT_EQ(store.sequenceNameStore[0], CharString("seq0"));
        SEQAN_ASSERT_EQ(store.sequenceNameStore[1], CharString("seq1"));

        SEQAN_ASSERT_EQ(length(store.matchStore), 2u);
        SEQAN_ASSERT_EQ(back(store.matchStore).id, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectId, 0u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectBeginPos, 5u);
        SEQAN_ASSERT_EQ(back(store.matchStore).subjectEndPos, 6u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryId, 1u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryBeginPos, 7u);
        SEQAN_ASSERT_EQ(back(store.matchStore).queryEndPos, 8u);
    }
}

SEQAN_DEFINE_TEST(test_parse_lm_parse_lastz_general)
{
    using namespace seqan;

    // Write out lastz result file to temporary file.
    String<char> content =
            "#score\tname1\tstrand1\tsize1\tzstart1\tend1\tname2\tstrand2\tsize2\tzstart2\tend2\tidentity\tidPct\tcoverage\tcovPct\n"
            "8823\tNC_001401.2\t+\t4679\t0\t125\tNC_001729.1\t+\t4726\t4601\t4726\t107/125\t85.6%%\t125/4679\t2.7%%\n"
            "303631\tNC_001401.2\t-\t4679\t0\t4679\tNC_001729.1\t+\t4726\t0\t4726\t3876/4656\t83.2%%\t4679/4679\t100.0%%\n";
    // Read file.
    LocalMatchStore<> store;
    DirectionIterator<String<char>, Input>::Type iter = begin(content);

    readRecord(store, iter, LastzGeneral());
    readRecord(store, iter, LastzGeneral());
    SEQAN_TEST_EXCEPTION(ParseError,
                         seqan::readRecord(store, iter, LastzGeneral()));


    SEQAN_ASSERT_EQ(length(store.matchStore), 2u);
    SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
    SEQAN_ASSERT_EQ(store.sequenceNameStore[0], "NC_001401.2");
    SEQAN_ASSERT_EQ(store.sequenceNameStore[1], "NC_001729.1");

    typedef typename LocalMatchStore<>::TLocalMatch TLocalMatch;

    TLocalMatch const match0(0, 0, 0, 125, 1, 4601, 4726);
    SEQAN_ASSERT(match0 == store.matchStore[0]);
    TLocalMatch const match1(1, 0, 4679, 0, 1, 0, 4726);
    SEQAN_ASSERT(match1 == store.matchStore[1]);
}

SEQAN_DEFINE_TEST(test_parse_lm_parse_blastn_tabular)
{
    using namespace seqan;

    // Write out lastz result file to temporary file.
    String<char> content =
            "# BLASTN 2.2.25+\n"
            "# Query: NC_001729.1 Adeno-associated virus - 3, complete genome\n"
            "# Subject: NC_001401.2 Adeno-associated virus - 2, complete genome\n"
            "# Fields: query id, subject id, %% identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n"
            "# 1 hits found\n"
            "NC_001729.1\tNC_001401.2\t83.13\t4617\t580\t182\t1\t4525\t1\t4510\t0.0\t4026\n"
            "# BLAST processed 1 queries\n";


    // Read file.
    LocalMatchStore<> store;
    DirectionIterator<String<char>, Input>::Type iter = begin(content);
    readRecord(store, iter, BlastnTabular());
    SEQAN_TEST_EXCEPTION(ParseError,
                         seqan::readRecord(store, iter, BlastnTabular()));

    SEQAN_ASSERT_EQ(length(store.matchStore), 1u);
    SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
    SEQAN_ASSERT_EQ(store.sequenceNameStore[0], "NC_001401.2");
    SEQAN_ASSERT_EQ(store.sequenceNameStore[1], "NC_001729.1");

    typedef typename LocalMatchStore<>::TLocalMatch TLocalMatch;

    TLocalMatch const match0(0, 0, 0, 4510, 1, 0, 4525);
    SEQAN_ASSERT(match0 == store.matchStore[0]);
}

SEQAN_DEFINE_TEST(test_parse_lm_parse_stellar_gff)
{
    using namespace seqan;

    // Write out lastz result file to temporary file.
    String<char> content =
            "NC_001401.2\tStellar\teps-matches\t790\t1012\t95.0673\t+\t.\tNC_001729.1;seq2Range=787,1009;eValue=4.76915e-92;cigar=223M;mutations=11C,20G,44G,50C,78A,93G,155A,167G,176G,182G,191C\n"
            "NC_001401.2\tStellar\teps-matches\t1371\t1484\t95.614\t-\t.\tNC_001729.1;seq2Range=1368,1481;eValue=5.10729e-45;cigar=111M1I1M;mutations=15T,54G,57C,90A,99G\n";

    // Read file.
    LocalMatchStore<> store;
    DirectionIterator<String<char>, Input>::Type iter = begin(content);
    readRecord(store, iter, StellarGff());
    readRecord(store, iter, StellarGff());
    SEQAN_TEST_EXCEPTION(ParseError,
                         seqan::readRecord(store, iter, StellarGff()));

    SEQAN_ASSERT_EQ(length(store.matchStore), 2u);
    SEQAN_ASSERT_EQ(length(store.sequenceNameStore), 2u);
    SEQAN_ASSERT_EQ(store.sequenceNameStore[0], "NC_001401.2");
    SEQAN_ASSERT_EQ(store.sequenceNameStore[1], "NC_001729.1");

    typedef typename LocalMatchStore<>::TLocalMatch TLocalMatch;

    TLocalMatch const match0(0, 0, 789, 1012, 1, 786, 1009);
    SEQAN_ASSERT(match0 == store.matchStore[0]);
    TLocalMatch const match1(1, 0, 1484, 1370, 1, 1367, 1481);
    SEQAN_ASSERT(match1 == store.matchStore[1]);

    SEQAN_ASSERT_EQ(length(store.cigarStore), 2u);
    SEQAN_ASSERT_EQ(length(store.cigarStore[0]), 1u);
    SEQAN_ASSERT_EQ(store.cigarStore[0][0].operation, 'M');
    SEQAN_ASSERT_EQ(store.cigarStore[0][0].count, 223u);
    SEQAN_ASSERT_EQ(length(store.cigarStore[1]), 3u);
    SEQAN_ASSERT_EQ(store.cigarStore[1][0].operation, 'M');
    SEQAN_ASSERT_EQ(store.cigarStore[1][0].count, 111u);
    SEQAN_ASSERT_EQ(store.cigarStore[1][1].operation, 'I');
    SEQAN_ASSERT_EQ(store.cigarStore[1][1].count, 1u);
    SEQAN_ASSERT_EQ(store.cigarStore[1][2].operation, 'M');
    SEQAN_ASSERT_EQ(store.cigarStore[1][2].count, 1u);
}

#endif  // TESTS_PARSE_LM_TEST_PARSE_LM_H_
