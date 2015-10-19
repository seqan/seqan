// ==========================================================================
//                                SeqAn-Flexbar
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
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================
// This file provides tests for the general processing functions of
// seqan-flexbar which is based in the implementation of the original
// flexbar program in [1].
// [1] Dodt, M.; Roehr, J.T.; Ahmed, R.; Dieterich, C.  FLEXBAR—Flexible
// Barcode and Adapter Processing for Next-Generation Sequencing Platforms.
// Biology 2012, 1, 895-905.
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "demultiplex.h"
#include "general_processing.h"
#include "read.h"

using namespace seqan;

SEQAN_DEFINE_TEST(removeShortSeqs_test)
{
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(5);

    reads[0].id = "1";
    reads[0].seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    reads[1].id = "2";
    reads[1].seq = "AAAAAAAAAAAAAAAAAAAA";
    reads[2].id = "3";
    reads[2].seq = "AAAAA";
    reads[3].id = "4";
    reads[3].seq = "AAA";
    reads[4].id = "5";
    reads[4].seq = "A";

    std::vector<ReadPairedEnd<seqan::Dna5QString>> readsPairedEnd(5);
    for (unsigned int i = 0;i < 5;++i)
    {
        readsPairedEnd[i].id = readsPairedEnd[i].idRev = reads[i].id;
        readsPairedEnd[i].seq = readsPairedEnd[i].seqRev = reads[i].seq;
    }

    unsigned removedReads = 0;

    removedReads += removeShortSeqs(reads, 2);
    SEQAN_ASSERT_EQ(length(reads), 4u);

    removedReads += removeShortSeqs(reads, 4);
    SEQAN_ASSERT_EQ(length(reads), 3u);

    removedReads += removeShortSeqs(reads, 6);
    SEQAN_ASSERT_EQ(length(reads), 2u);

    removedReads += removeShortSeqs(reads, 24);
    SEQAN_ASSERT_EQ(length(reads), 1u);

    removedReads += removeShortSeqs(reads, 50);
    SEQAN_ASSERT_EQ(length(reads), 0u);

    //Part for paired End
    removedReads = 0;

    removedReads += removeShortSeqs(readsPairedEnd, 2);
    SEQAN_ASSERT_EQ(length(readsPairedEnd), 4u);
    SEQAN_ASSERT_EQ(removedReads, 1u);

    removedReads += removeShortSeqs(readsPairedEnd, 4);
    SEQAN_ASSERT_EQ(length(readsPairedEnd), 3u);
    SEQAN_ASSERT_EQ(removedReads, 2u);

    removedReads += removeShortSeqs(readsPairedEnd, 6);
    SEQAN_ASSERT_EQ(length(readsPairedEnd), 2u);
    SEQAN_ASSERT_EQ(removedReads, 3u);

    removedReads += removeShortSeqs(readsPairedEnd, 24);
    SEQAN_ASSERT_EQ(length(readsPairedEnd), 1u);
    SEQAN_ASSERT_EQ(removedReads, 4u);

    removedReads += removeShortSeqs(readsPairedEnd, 50);
    SEQAN_ASSERT_EQ(length(readsPairedEnd), 0u);
    SEQAN_ASSERT_EQ(removedReads, 5u);
}


SEQAN_DEFINE_TEST(findN_test)
{
    std::vector<Read<seqan::Dna5QString>> reads(12);
    reads[0].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seq = "ATGGNGGGTACACGTGATCGTACGTAGCAGC";
    reads[2].seq = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seq = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seq = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seq = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    reads[6].seq = "ACTGTACGTGATCG.AATGCTGACTGACTGAC";
    reads[7].seq = ".ACTGTACGTGATCG.AATGCTGACTGACTGAC";
    reads[8].seq = "ACT.GTACGTGATCG.AATGCTGACTGA.CTGAC";
    reads[9].seq = "ACT.GTACGTGATCG.AATGC.TGACTGACTG.AC";
    reads[10].seq = ".......";
    reads[11].seq = "";
    unsigned allowed = 3;
    StringSet<int> exspectedNoSub;
    appendValue(exspectedNoSub, 0);
    appendValue(exspectedNoSub, 1);
    appendValue(exspectedNoSub, 2);
    appendValue(exspectedNoSub, 3);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, 1);
    appendValue(exspectedNoSub, 2);
    appendValue(exspectedNoSub, 3);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, -1);
    appendValue(exspectedNoSub, 0);
    std::vector<int> res(length(reads));
    for (unsigned i = 0; i < length(reads); ++i)
    {
        res[i] = findN(reads[i], allowed, NoSubstitute());           //Test WITHOUT substitutions
    }
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], res[i]);
    }

    auto reads2 = reads;
    reads2[1].seq[4] = 'A';
    reads2[2].seq[0] = 'A';
    reads2[2].seq[31] = 'A';
    reads2[3].seq[0] = 'A';
    reads2[3].seq[12] = 'A';
    reads2[3].seq[31] = 'A';
    reads2[4].seq[0] = 'A';
    reads2[4].seq[12] = 'A';
    reads2[4].seq[20] = 'A';
    reads2[5].seq[0] = 'A';
    reads2[5].seq[1] = 'A';
    reads2[5].seq[2] = 'A';

     for (unsigned i = 0; i < length(reads); ++i)
    {
        res[i] = findN(reads[i], allowed, 'A');      //Test WITH substitutions
    }
    for (unsigned i = 0; i < 5; ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], res[i]);
    }
}

SEQAN_DEFINE_TEST(processN_test)
{
    GeneralStats stats;
    std::vector<Read<seqan::Dna5QString>> reads(6);
    reads[0].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seq = "ATGGNGGGTACACGTGATCGTACGTAGCAGC";
    reads[2].seq = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seq = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seq = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].id = "Null";
    reads[1].id = "Eins";
    reads[2].id = "Zwei";
    reads[3].id = "Drei";
    reads[4].id = "Vier";
    reads[5].id = "Null2";
    auto reads2 = reads;

    unsigned allowed = 3;
    StringSet<String<Dna5> > exspectedNoSub;
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "ATGGNGGGTACACGTGATCGTACGTAGCAGC");
    appendValue(exspectedNoSub, "NGGGACTGTACACGTGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "NATGACTGTAGGNGGGATCGTACGTAGCAGGN");
    appendValue(exspectedNoSub, "ATGACTGTACACGTGATCGTACGTAGCAGC");

    processN(reads2, allowed, NoSubstitute(), stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedN, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedNoSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedNoSub[i], reads2[i].seq);
    }

    stats.removedN = 0;
    stats.uncalledBases = 0;
    reads2 = reads;
    Dna substitute = 'A';
    StringSet<String<Dna5> > exspectedSub = exspectedNoSub;
    exspectedSub[1][4] = substitute;
    exspectedSub[2][0] = substitute;
    exspectedSub[2][31] = substitute;
    exspectedSub[3][0] = substitute;
    exspectedSub[3][12] = substitute;
    exspectedSub[3][31] = substitute;

    processN(reads2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedN, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(exspectedSub); ++i)
    {
        SEQAN_ASSERT_EQ(exspectedSub[i], reads2[i].seq);
    }
}

SEQAN_DEFINE_TEST(processN_paired_test)
{
    GeneralStats stats;
    
    std::vector<ReadPairedEnd<seqan::Dna5QString>> reads(6);
    reads[0].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seq = "ATGGNGGGTACACGTGATCGTACGTAGCAGC";
    reads[2].seq = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seq = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seq = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].seqRev = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seqRev = "ATGNGNGGGTANCACGTGATCGTNACGTAGCANGC";
    reads[2].seqRev = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seqRev = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seqRev = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seqRev = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].id = "Null";
    reads[1].id = "Eins";
    reads[2].id = "Zwei";
    reads[3].id = "Drei";
    reads[4].id = "Vier";
    reads[5].id = "Null2";

    auto reads2 = reads;
    auto expectedReads = reads;

    expectedReads.erase(expectedReads.begin() + 1);
    expectedReads.erase(expectedReads.begin() + 3);

    unsigned allowed = 3;
    
    processN(reads2, allowed, NoSubstitute(), stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedN, 2u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(reads2[i].id, expectedReads[i].id);
        SEQAN_ASSERT_EQ(reads2[i].seq, expectedReads[i].seq);
        SEQAN_ASSERT_EQ(reads2[i].seqRev, expectedReads[i].seqRev);
    }

    stats.removedN = 0;
    stats.uncalledBases = 0;
    Dna substitute = 'A';
    reads2 = reads;
    expectedReads[1].seq[0] = substitute;
    expectedReads[1].seq[31] = substitute;
    expectedReads[2].seq[0] = substitute;
    expectedReads[2].seq[12] = substitute;
    expectedReads[2].seq[31] = substitute;

    processN(reads2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedN, 2u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(expectedReads[i].seq, reads2[i].seq);
    }
}

SEQAN_DEFINE_TEST(processN_multiplex_test)
{
    GeneralStats stats;
    std::vector<ReadMultiplex<seqan::Dna5QString>> reads(6);
    reads[0].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seq = "ATGGNGGGTACACGTGATCGTACGTAGCAGC";
    reads[2].seq = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seq = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seq = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].demultiplex = "ACTGTA";
    reads[1].demultiplex = "TGACGT";
    reads[2].demultiplex = "GTACGA";
    reads[3].demultiplex = "GTACTG";
    reads[4].demultiplex = "AAAAAA";
    reads[5].demultiplex = "GGGTAC";

    reads[0].id = "Null";
    reads[1].id = "Eins";
    reads[2].id = "Zwei";
    reads[3].id = "Drei";
    reads[4].id = "Vier";
    reads[5].id = "Null2";
    auto reads2 = reads;
    auto expectedReads = reads;
    expectedReads.erase(expectedReads.begin() + 4);

    unsigned allowed = 3;

    processN(reads2, allowed, NoSubstitute(), stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedN, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(expectedReads[i].seq, reads2[i].seq);
        SEQAN_ASSERT_EQ(expectedReads[i].demultiplex, reads2[i].demultiplex);
    }

    stats.removedN = 0;
    stats.uncalledBases = 0;
    reads2 = reads;
    Dna substitute = 'A';
    expectedReads[1].seq[4] = substitute;
    expectedReads[2].seq[0] = substitute;
    expectedReads[2].seq[31] = substitute;
    expectedReads[3].seq[0] = substitute;
    expectedReads[3].seq[12] = substitute;
    expectedReads[3].seq[31] = substitute;

    processN(reads2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedN, 1u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 6u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(expectedReads[i].seq, reads2[i].seq);
        SEQAN_ASSERT_EQ(expectedReads[i].demultiplex, reads2[i].demultiplex);
    }
}

SEQAN_DEFINE_TEST(processN_paired_multiplex_test)
{
    GeneralStats stats;
    std::vector<ReadMultiplexPairedEnd<seqan::Dna5QString>> reads(6);
    reads[0].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seq = "ATGGNGGGTACACGTGATCGTACGTAGCAGC";
    reads[2].seq = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seq = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seq = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seq = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].seqRev = "ATGACTGTACACGTGATCGTACGTAGCAGC";
    reads[1].seqRev = "ATGNGNGGGTANCACGTGATCGTNACGTAGCANGC";
    reads[2].seqRev = "NGGGACTGTACACGTGATCGTACGTAGCAGGN";
    reads[3].seqRev = "NATGACTGTAGGNGGGATCGTACGTAGCAGGN";
    reads[4].seqRev = "NGGGACTGTAGGNGGGATGGNGGGTAGCAGGN";
    reads[5].seqRev = "ATGACTGTACACGTGATCGTACGTAGCAGC";

    reads[0].demultiplex = "ACTGTA";
    reads[1].demultiplex = "TGACGT";
    reads[2].demultiplex = "GTACGA";
    reads[3].demultiplex = "GTACTG";
    reads[4].demultiplex = "AAAAAA";
    reads[5].demultiplex = "GGGTAC";

    reads[0].id = "Null";
    reads[1].id = "Eins";
    reads[2].id = "Zwei";
    reads[3].id = "Drei";
    reads[4].id = "Vier";
    reads[5].id = "Null2";
    auto reads2 = reads;
    auto expectedReads = reads;
    expectedReads.erase(expectedReads.begin() + 1);
    expectedReads.erase(expectedReads.begin() + 3);

    unsigned allowed = 3;
    
    processN(reads2, allowed, NoSubstitute(), stats);        //No Substitutions
    SEQAN_ASSERT_EQ(stats.removedN, 2u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(reads2[i].seq, expectedReads[i].seq);
        SEQAN_ASSERT_EQ(reads2[i].id, expectedReads[i].id);
        SEQAN_ASSERT_EQ(reads2[i].seqRev, expectedReads[i].seqRev);
        SEQAN_ASSERT_EQ(reads2[i].demultiplex, expectedReads[i].demultiplex);
    }

    stats.removedN = 0;
    stats.uncalledBases = 0;
    reads2 = reads;
    Dna substitute = 'A';
    expectedReads[1].seq[0] = substitute;
    expectedReads[1].seq[31] = substitute;
    expectedReads[2].seq[0] = substitute;
    expectedReads[2].seq[12] = substitute;
    expectedReads[2].seq[31] = substitute;

    processN(reads2, allowed, substitute, stats);
    SEQAN_ASSERT_EQ(stats.removedN, 2u);
    SEQAN_ASSERT_EQ(stats.uncalledBases, 10u);
    for (unsigned i = 0; i < length(expectedReads); ++i)
    {
        SEQAN_ASSERT_EQ(reads2[i].seq, expectedReads[i].seq);
        SEQAN_ASSERT_EQ(reads2[i].demultiplex, expectedReads[i].demultiplex);
    }
}

SEQAN_DEFINE_TEST(preTrim_test)
{
    GeneralStats stats;
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(7);

    reads[0].seq = "ACGTAACTGA";
    reads[1].seq = "AAAAAACTTTTT";
    reads[2].seq = "AAAAAAG";
    reads[3].seq = "TACGG";
    reads[4].seq = "TAAAAAA";
    reads[5].seq = "";
    reads[6].seq = "GATTACAGATTACA";

    reads[0].id = "loeschenTrim";
    reads[1].id = "Head/Tail";
    reads[2].id = "Head";
    reads[3].id = "loeschenNone";
    reads[4].id = "Tail";
    reads[5].id = "loeschenLeer";
    reads[6].id = "None";

    auto reads2 = reads;

    preTrim(reads2, 3, 3, 4, false, stats);
    SEQAN_ASSERT_EQ(reads2[0].seq, "TAAC");
    SEQAN_ASSERT_EQ(reads2[0].id, "loeschenTrim");
    SEQAN_ASSERT_EQ(reads2[1].seq, "AAACTT");
    SEQAN_ASSERT_EQ(reads2[1].id, "Head/Tail");
    SEQAN_ASSERT_EQ(reads2[2].seq, "TACAGATT");
    SEQAN_ASSERT_EQ(reads2[2].id, "None");
    SEQAN_ASSERT_EQ(length(reads2), 3u);

    reads2 = reads;
    preTrim(reads2, 6, 5, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[0].seq, "C");
    SEQAN_ASSERT_EQ(reads2[0].id, "Head/Tail");
    SEQAN_ASSERT_EQ(reads2[1].seq, "AGA");
    SEQAN_ASSERT_EQ(reads2[1].id, "None");
    SEQAN_ASSERT_EQ(length(reads2), 2u);

    reads2 = reads;
    preTrim(reads2, 6, 0, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[2].seq, "G");
    SEQAN_ASSERT_EQ(reads2[2].id, "Head");

    reads2 = reads;
    preTrim(reads2, 0, 6, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[3].seq, "T");
    SEQAN_ASSERT_EQ(reads2[3].id, "Tail");

    reads2 = reads;
    preTrim(reads2, 0, 0, 6, false, stats);
    SEQAN_ASSERT_EQ(reads2[3].seq, "TAAAAAA");
    SEQAN_ASSERT_EQ(reads2[3].id, "Tail");
    SEQAN_ASSERT_EQ(reads2[4].seq, "GATTACAGATTACA");
    SEQAN_ASSERT_EQ(reads2[4].id, "None");
    SEQAN_ASSERT_EQ(length(reads2), 5u);
}

SEQAN_DEFINE_TEST(preTrim_paired_test)
{
    GeneralStats stats;
    using TRead = ReadPairedEnd<seqan::Dna5QString>;
    std::vector<TRead> reads(7);
    std::vector<TRead> reads2(7);

    reads[0].seq = "ACGTAACTGA";
    reads[1].seq = "AAAAAACTTTTT";
    reads[2].seq = "AAAAAAG";
    reads[3].seq = "TACGG";
    reads[4].seq = "TAAAAAA";
    reads[5].seq = "";
    reads[6].seq = "GATTACAGATTACA";

    reads[0].id = "loeschenTrim";
    reads[1].id = "Head/Tail";
    reads[2].id = "Head";
    reads[3].id = "loeschenNone";
    reads[4].id = "Tail";
    reads[5].id = "loeschenLeer";
    reads[6].id = "None";
    for (unsigned int i = 0;i < 7;++i)
    {
        reads[i].seqRev = reads[i].seq;
        reads[i].idRev = reads[i].id;
    }
    reads[6].seqRev = "GTCA";

    reads2 = reads;
    preTrim(reads2, 3, 3, 4, false, stats);
    SEQAN_ASSERT_EQ(reads2[0].seq, "TAAC");
    SEQAN_ASSERT_EQ(reads2[0].id, "loeschenTrim");
    SEQAN_ASSERT_EQ(reads2[1].seq, "AAACTT");
    SEQAN_ASSERT_EQ(reads2[1].id, "Head/Tail");
    SEQAN_ASSERT_EQ(length(reads2), 2u);

    reads2 = reads;
    preTrim(reads2, 6, 5, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[0].seq, "C");
    SEQAN_ASSERT_EQ(reads2[0].id, "Head/Tail");
    SEQAN_ASSERT_EQ(length(reads2), 1u);

    reads2 = reads;
    preTrim(reads2, 6, 0, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[2].seq, "G");
    SEQAN_ASSERT_EQ(reads2[2].id, "Head");

    reads2 = reads;
    preTrim(reads2, 0, 6, 1, false, stats);
    SEQAN_ASSERT_EQ(reads2[3].seq, "T");
    SEQAN_ASSERT_EQ(reads2[3].id, "Tail");

    reads2 = reads;
    preTrim(reads2, 0, 0, 6, false, stats);
    SEQAN_ASSERT_EQ(reads2[3].seq, "TAAAAAA");
    SEQAN_ASSERT_EQ(reads2[3].id, "Tail");
    SEQAN_ASSERT_EQ(length(reads2), 4u);
}

SEQAN_DEFINE_TEST(trimTo_test)
{
    GeneralStats stats;
    using TRead = Read<seqan::Dna5QString>;
    std::vector<TRead> reads(4);
    
    reads[0].seq = "123456789";
    reads[1].seq = "123456";
    reads[2].seq = "1234567";
    reads[3].seq = "";

    reads[0].id = "neun";
    reads[1].id = "sechs";
    reads[2].id = "sieben";
    reads[3].id = "null";

    trimTo(reads, 7, stats);

    SEQAN_ASSERT_EQ(reads[0].seq, "1234567");
    SEQAN_ASSERT_EQ(reads[1].seq, "1234567");
    SEQAN_ASSERT_EQ(reads[0].id, "neun");
    SEQAN_ASSERT_EQ(reads[1].id, "sieben");
    SEQAN_ASSERT_EQ(length(reads), 2u);
}

SEQAN_DEFINE_TEST(trimTo_paired_test)
{
    GeneralStats stats;
    using TRead = ReadPairedEnd<seqan::Dna5QString>;
    std::vector<TRead> reads(4);

    reads[0].seq = "123456789";
    reads[1].seq = "123456";
    reads[2].seq = "1234567";
    reads[3].seq = "";

    reads[0].id = "neun";
    reads[1].id = "sechs";
    reads[2].id = "sieben";
    reads[3].id = "null";
    for (unsigned int i = 0;i < 4;++i)
    {
        reads[i].seqRev = reads[i].seq;
        reads[i].idRev = reads[i].id;
    }
    reads[0].seqRev = "123456";

    trimTo(reads, 7, stats);

    SEQAN_ASSERT_EQ(reads[0].seq, "1234567");
    SEQAN_ASSERT_EQ(reads[0].id, "sieben");
    SEQAN_ASSERT_EQ(length(reads), 1u);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(removeShortSeqs_test);
    SEQAN_CALL_TEST(findN_test);
    SEQAN_CALL_TEST(processN_test);
    SEQAN_CALL_TEST(processN_paired_test);
    SEQAN_CALL_TEST(processN_multiplex_test);
    SEQAN_CALL_TEST(processN_paired_multiplex_test);
    SEQAN_CALL_TEST(preTrim_test);
    SEQAN_CALL_TEST(preTrim_paired_test);
    SEQAN_CALL_TEST(trimTo_test);
    SEQAN_CALL_TEST(trimTo_paired_test);
}
SEQAN_END_TESTSUITE
