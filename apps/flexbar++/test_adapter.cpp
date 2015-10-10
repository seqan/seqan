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
// Author: Benjamin Strauch <b.strauch@fu-berlin.de>
// Author: Sebastian Roskosch <serosko@zedat.fu-berlin.de>
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/file.h>
#include "adapter_trimming.h"
#include "read.h"

SEQAN_DEFINE_TEST(get_overlap_test)
{
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAAAAAA";
    TSequence seq2 = "TTTTTTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Overlap of 4.
    // ---AAAAAAAAA
    // TTTTTTT-----
    seqan::insertGaps(row1, 0, 3);
    seqan::insertGaps(row2, length(row2), 5);

    SEQAN_ASSERT_EQ(getOverlap(align), 4u);

    // Overlap of 0.
    // -------AAAAAAAAA
    // TTTTTTT---------
    seqan::clearGaps(align);
    seqan::insertGaps(row1, 0, 7);
    seqan::insertGaps(row2, length(row2), 9);

    SEQAN_ASSERT_EQ(getOverlap(align), 0u);

    // Overlap of 7.
    // AAAAAAAAA
    // -TTTTTTT-
    seqan::clearGaps(align);
    seqan::insertGap(row2, 0);
    seqan::insertGap(row2, length(row2));

    SEQAN_ASSERT_EQ(getOverlap(align), 7u);
}

SEQAN_DEFINE_TEST(count_gap_test)
{
	seqan::Dna5String seq("AAGTCTATCTA");
	seqan::Gaps<seqan::Dna5String> row(seq);

	// no gaps yet.
	SEQAN_ASSERT_EQ(countTotalGaps(row), 0u);

	// insert a few gaps at the start.
	seqan::insertGaps(row, 0, 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 10u);

	// create multiple continuous gaps

	seqan::insertGaps(row, 15, 10);
	seqan::insertGaps(row, length(row), 10);
	SEQAN_ASSERT_EQ(countTotalGaps(row), 30u);
}

SEQAN_DEFINE_TEST(insert_size_test)
{
	// The insert method is specified only for actually
	// overlapping sequences so we test those.
	typedef seqan::String<seqan::Dna5> TSequence;
	typedef seqan::Align<TSequence, seqan::ArrayGaps > TAlign;
	typedef seqan::Row<TAlign>::Type TRow;

    TSequence seq1 = "AAAAA";
    TSequence seq2 = "TTT";

    TAlign align;
    resize(rows(align), 2);
    seqan::assignSource(row(align,0),seq1);
    seqan::assignSource(row(align,1),seq2);

    TRow &row1 = row(align, 0);
    TRow &row2 = row(align, 1);

    // Insert size is 7
    // AAAAA--
    // ----TTT
    seqan::insertGaps(row1, seqan::length(row1), 2);
    seqan::insertGaps(row2, 0, 4);
    SEQAN_ASSERT_EQ(getInsertSize(align), 7u);

    // Insert size is 5
    // AAAAA
    // --TTT
    seqan::clearGaps(align);
    seqan::insertGaps(row2, 0, 2);
    SEQAN_ASSERT_EQ(getInsertSize(align), 5u);

    // Insert size is 4 (seq1 goes one base into the adapter).
	// AAAA|A
	// -TTT|-
	seqan::clearGaps(align);
	seqan::insertGap(row2, seqan::length(row2));
	seqan::insertGap(row2, 0);
	SEQAN_ASSERT_EQ(getInsertSize(align), 4u);
}

// Returns (num2 - num1 + 1) in a string "xxxx_num1_num2_xxx...".
int insertSize(char* s)
{
	int insert = 0;
	strtok(s, "_");
	insert -= atoi(strtok(NULL, "_"));
	insert += atoi(strtok(NULL, "_")) + 1;
	return insert;
}

SEQAN_DEFINE_TEST(match_test)
{
    AdapterMatchSettings u(4, 0, 0.2, 0, 1);

	SEQAN_ASSERT(isMatch(4,0,u));
	SEQAN_ASSERT_NOT(isMatch(4,1,u));

    SEQAN_ASSERT(isMatch(10,2,u));
	SEQAN_ASSERT_NOT(isMatch(10,3,u));
	// 20% errors are allowed.
	SEQAN_ASSERT(isMatch(100,20,u));
	SEQAN_ASSERT_NOT(isMatch(100,21,u));

    AdapterMatchSettings u2(7, 2, 0, 0, 1);
    // We need an overlap of length 7 and no more than 2 errors.
	SEQAN_ASSERT_NOT(isMatch(5,3,u2));
	SEQAN_ASSERT_NOT(isMatch(7,3,u2));
	SEQAN_ASSERT(isMatch(7,2,u2));
}

SEQAN_DEFINE_TEST(strip_adapter_test)
{
    typedef seqan::String<seqan::Dna5Q> TAda;
    using AdapterSet = std::vector<AdapterItem>;
    AdapterSet adapterSet;
    using TRead = Read<seqan::Dna5QString>;
    TRead read;

    read.seq = "AAAAAAAAAATTTTT";
    TAda ada = TAda(     "TTTTTTTTTTT");

	int len = length(read.seq);
    // overlap=4, error_rate = 0.2, times = 1
    AdapterMatchSettings matchSettings(4, 0, 0.2, 0, 1);
    AdapterTrimmingStats stats;
 
    int removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end3, 0,0, false, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
	SEQAN_ASSERT_EQ(removed, 5);
	SEQAN_ASSERT_EQ(len - length(read.seq), 5u);

    read.seq = "AAAAAGAAAATATATTA";
    //               ||| |||||||		   
    ada = TAda(     "GAATATATATTT");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end3, 0,0, false, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 12);
    SEQAN_ASSERT_EQ(len - length(read.seq), 12u);

    read.seq = "AAAAAGAATATATATA";
    //               ||||||||||		   
    ada = TAda(     "GAATATATATTT");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end3, 0,0, false, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 11);
    SEQAN_ASSERT_EQ(len - length(read.seq), 11u);

    read.seq = "AAAAAGAATATATATA";
    //               ||||||||||		   
    ada = TAda(     "GAATATATATTT");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end5, 0,0, false, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 0);
    SEQAN_ASSERT_EQ(len - length(read.seq), 0u);

    read.seq =      "AAAAAGAATATATATA";
    //               ||||||||		   
    ada = TAda("CCCCAAAAAAGAAC");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end5, 0,0, false, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 9);
    SEQAN_ASSERT_EQ(len - length(read.seq), 9u);

    // test anchored mode for 5-end adapters
    read.seq = "AAAAAGAATATATATA";
    //          ||||||||		   
    ada = TAda("AAAAAGAAC");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end5, 0,0, true, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 9);
    SEQAN_ASSERT_EQ(len - length(read.seq), 9u);

    read.seq =  "AAAAAGAATATATATA";
    //           ||||||||		   
    ada = TAda("CAAAAAGAAC");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end5, 0,0, true, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 0);
    SEQAN_ASSERT_EQ(len - length(read.seq), 0u);

    // test anchored mode for 3-end adapters
    read.seq = "AAAAAGAATATATATA";
    //               ||||||||||		   
    ada = TAda(     "GAATATATATT");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end3, 0,0, true, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 11);
    SEQAN_ASSERT_EQ(len - length(read.seq), 11u);

    read.seq = "AAAAAGAATATATATA";
    //               ||||||||||		   
    ada = TAda(     "GAATATATATTT");
    len = length(read.seq);
    removed = stripAdapter(read.seq, stats, AdapterSet{ AdapterItem(ada, AdapterItem::end3, 0,0, true, false) }, matchSettings, StripAdapterDirection<adapterDirection::forward>());
    SEQAN_ASSERT_EQ(removed, 0);
    SEQAN_ASSERT_EQ(len - length(read.seq), 0u);
}

SEQAN_DEFINE_TEST(align_adapter_test)
{
	typedef seqan::String<seqan::Dna5Q> TSeq;
	typedef seqan::String<seqan::Dna5> TAda;

	TSeq seq = TSeq("AAAAAAAAAATTTTT");
	TAda ada = TAda("TTTTTTTTTTT");
	std::pair<unsigned, seqan::Align<TSeq> > pair;
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>());
	SEQAN_ASSERT_EQ(pair.first, 5u);

	seq = TSeq("AAAAAAAAAATATATTA");
	//                    |||||		   
	ada = TAda(       "GGTTATATATTT"); // front and back gaps are allowed
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>());
	SEQAN_ASSERT_EQ(pair.first, 2u);

	seq = TSeq("AAAAAAAAAATATATTA");
	//                || |||||||		   
	ada = TAda(     "GAATATATATTT"); // front and back gaps are allowed
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>());
	SEQAN_ASSERT_EQ(pair.first, 6u);

    unsigned int rightOverhang = 4;
    unsigned int leftOverhang = 4;
    seq = TSeq("CATCATAAAAAATATATTA");
    //          ||||||		   
    ada = TAda("CATCAT"); 
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>(), leftOverhang, rightOverhang);
    SEQAN_ASSERT_EQ(pair.first, 6u);

    // just enough overlap
    seq = TSeq(    "CATCATAAAAAATATATTA");
    //              ||||		   
    ada = TAda("GGGGCATC"); 
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>(), leftOverhang, rightOverhang);
    SEQAN_ASSERT_EQ(pair.first, 4u);

    // not enough overlap, should report score 0
    seq = TSeq(     "CATCATAAAAAATATATTA");
    //               |||		   
    ada = TAda("GGGGGCAT");
    alignPair(pair, seq, ada, seqan::AlignConfig<true, true, true, true>(), leftOverhang, rightOverhang);
    SEQAN_ASSERT_EQ(pair.first, 0u);
}

SEQAN_DEFINE_TEST(strip_pair_test)
{
	typedef seqan::String<seqan::Dna5Q> TSeq;

	TSeq seq1("AAAAAAAAAATTTTTTT");
	TSeq seq2("CCCCCCCCCCAAAAAAA");
	TSeq seq3("AAAAAAAAAAGGGGGGG");

	unsigned insert = stripPair(seq1,seq2);
	SEQAN_ASSERT_EQ(insert, 27u);

	unsigned insert2 = stripPair(seq1,seq3);
	SEQAN_ASSERT_EQ(insert2, 0u);
}

SEQAN_BEGIN_TESTSUITE(test_my_app_funcs)
{
    SEQAN_CALL_TEST(get_overlap_test);
	SEQAN_CALL_TEST(count_gap_test);
	SEQAN_CALL_TEST(insert_size_test);
	SEQAN_CALL_TEST(match_test);
	SEQAN_CALL_TEST(strip_adapter_test);
	SEQAN_CALL_TEST(align_adapter_test);
	SEQAN_CALL_TEST(strip_pair_test);
}
SEQAN_END_TESTSUITE
