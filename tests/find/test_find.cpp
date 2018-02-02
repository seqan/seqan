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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Tests for the SeqAn module find.
// ==========================================================================

#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstring>  // size_t
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
//#define SEQAN_TEST

#define SEQAN_DEBUG_MYERSBITVECTOR
//#define SEQAN_DEBUG_MYERSBITVECTOR_DUMP
//#define SEQAN_TEST_MYERS_STRICTBANDED

#include <seqan/basic.h>
#include <seqan/find.h>

#include "test_find_hamming.h"
#include "test_find_myers_banded.h"

using namespace std;
using namespace seqan;


template <typename TAlgorithmSpec>
void Test_OnlineAlg() {
		typedef typename Position<CharString>::Type TPosition;

    String<TPosition> pos;

    //____________________________________________________________________________
    // Test1 - small needle

    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    String<char> needle("ist");
    Pattern<String<char>, TAlgorithmSpec> pattern(needle);

    while (find(finder, pattern)) {
        appendValue(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(needle));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), needle);
    }

    SEQAN_ASSERT_EQ(host(pattern), needle);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), needle);
    SEQAN_ASSERT_EQ(pos[0], 5u);
    SEQAN_ASSERT_EQ(pos[1], 31u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test2 - large needle

    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstuvwxyzabcdefg";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(needle));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), needle);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 26u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test3 - different alphabet, small needle

    String<Dna> hstk = "aaaaaaacaa";
    Finder<String<Dna> > finderDna(hstk);

    String<Dna> ndl = "aa";
    setHost(pattern, ndl);

    clear(pos);
    while (find(finderDna, pattern)) {
        append(pos,position(finderDna));
        SEQAN_ASSERT_EQ(position(finderDna), beginPosition(finderDna));
        SEQAN_ASSERT_EQ(endPosition(finderDna), beginPosition(finderDna) + length(finderDna));
        SEQAN_ASSERT_EQ(length(finderDna), length(ndl));
        SEQAN_ASSERT_EQ(begin(finderDna), begin(hstk) + beginPosition(finderDna));
        SEQAN_ASSERT_EQ(end(finderDna), begin(hstk) + endPosition(finderDna));
        SEQAN_ASSERT_EQ(infix(finderDna), ndl);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 1u);
    SEQAN_ASSERT_EQ(pos[2], 2u);
    SEQAN_ASSERT_EQ(pos[3], 3u);
    SEQAN_ASSERT_EQ(pos[4], 4u);
    SEQAN_ASSERT_EQ(pos[5], 5u);
    SEQAN_ASSERT_EQ(pos[6], 8u);
    SEQAN_ASSERT_EQ(length(pos), 7u);

    //____________________________________________________________________________
    // Test3b - different alphabet, small needle, jumping finder

    goBegin(finderDna); // That's a repositioning
    clear(finderDna);       // That's why, clear state
    clear(pos);

    bool firstHit = true;
    while (find(finderDna, pattern)) {
        if (firstHit) {
            firstHit = false;
            finderDna += 2;
            clear(finderDna);  // clear the state of the finder
        } else {
            //unsigned int p = position(finderDna);
            append(pos,position(finderDna));
        }
    }

    SEQAN_ASSERT_EQ(pos[0], 2u);
    SEQAN_ASSERT_EQ(pos[1], 3u);
    SEQAN_ASSERT_EQ(pos[2], 4u);
    SEQAN_ASSERT_EQ(pos[3], 5u);
    SEQAN_ASSERT_EQ(pos[4], 8u);
    SEQAN_ASSERT_EQ(length(pos), 5u);

    //____________________________________________________________________________
    // Test4 - different alphabet, large needle
    String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    Finder<String<Dna> > finderText(text);

    String<Dna> query = "taaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    setHost(pattern, query);

    clear(pos);
    while (find(finderText, pattern)) {
        append(pos,position(finderText));
        SEQAN_ASSERT_EQ(position(finderText), beginPosition(finderText));
        SEQAN_ASSERT_EQ(endPosition(finderText), beginPosition(finderText) + length(finderText));
        SEQAN_ASSERT_EQ(length(finderText), length(query));
        SEQAN_ASSERT_EQ(begin(finderText), begin(text) + beginPosition(finderText));
        SEQAN_ASSERT_EQ(end(finderText), begin(text) + endPosition(finderText));
        SEQAN_ASSERT_EQ(infix(finderText), query);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 5u);
    SEQAN_ASSERT_EQ(pos[2], 10u);
    SEQAN_ASSERT_EQ(pos[3], 15u);
    SEQAN_ASSERT_EQ(pos[4], 20u);
    SEQAN_ASSERT_EQ(pos[5], 25u);
    SEQAN_ASSERT_EQ(length(pos), 6u);
}


template <typename TAlgorithmSpec>
void Test_OnlineAlgMulti(bool order_by_begin_position) {
	typedef typename Position<CharString>::Type TPosition;

    String<TPosition> pos;

    //____________________________________________________________________________
    // Test1 - Single keyword
    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    typedef String<String<char> > TNeedle;
    TNeedle keywords;
    appendValue(keywords, String<char>("ist"));
    Pattern<TNeedle, TAlgorithmSpec> pattern(keywords);

    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(keywords[position(pattern)]));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), keywords[position(pattern)]);
    }

    SEQAN_ASSERT_EQ(host(pattern), keywords);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<TNeedle, TAlgorithmSpec> const &>(pattern)), keywords);
    SEQAN_ASSERT_EQ(pos[0], 5u);
    SEQAN_ASSERT_EQ(pos[1], 31u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    clear(keywords);
    appendValue(keywords, String<char>("abcdefghijklmnopqrstuvwxyzabcdefg"));
    setHost(pattern, keywords);
    clear(pos);

    while (find(finder, pattern)) {
        append(pos,position(finder));
        SEQAN_ASSERT_EQ(position(finder), beginPosition(finder));
        SEQAN_ASSERT_EQ(endPosition(finder), beginPosition(finder) + length(finder));
        SEQAN_ASSERT_EQ(length(finder), length(keywords[position(pattern)]));
        SEQAN_ASSERT_EQ(begin(finder), begin(haystack) + beginPosition(finder));
        SEQAN_ASSERT_EQ(end(finder), begin(haystack) + endPosition(finder));
        SEQAN_ASSERT_EQ(infix(finder), keywords[position(pattern)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 26u);
    SEQAN_ASSERT_EQ(length(pos), 2u);


    String<Dna> hstk = "aaaaaaacaa";
    Finder<String<Dna> > finderDna(hstk);

    typedef String<String<Dna> > TDnaNeedle;
    Pattern<TDnaNeedle, TAlgorithmSpec> pattern_dna(keywords);

    TDnaNeedle dna_keywords;
    appendValue(dna_keywords, String<Dna>("aa"));
    setHost(pattern_dna, dna_keywords);

    clear(pos);
    while (find(finderDna, pattern_dna)) {
        append(pos,position(finderDna));
        SEQAN_ASSERT_EQ(position(finderDna), beginPosition(finderDna));
        SEQAN_ASSERT_EQ(endPosition(finderDna), beginPosition(finderDna) + length(finderDna));
        SEQAN_ASSERT_EQ(length(finderDna), length(dna_keywords[position(pattern_dna)]));
        SEQAN_ASSERT_EQ(begin(finderDna), begin(hstk) + beginPosition(finderDna));
        SEQAN_ASSERT_EQ(end(finderDna), begin(hstk) + endPosition(finderDna));
        SEQAN_ASSERT_EQ(infix(finderDna), dna_keywords[position(pattern_dna)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 1u);
    SEQAN_ASSERT_EQ(pos[2], 2u);
    SEQAN_ASSERT_EQ(pos[3], 3u);
    SEQAN_ASSERT_EQ(pos[4], 4u);
    SEQAN_ASSERT_EQ(pos[5], 5u);
    SEQAN_ASSERT_EQ(pos[6], 8u);
    SEQAN_ASSERT_EQ(length(pos), 7u);

    String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
    Finder<String<Dna> > finderText(text);

    clear(dna_keywords);
    appendValue(dna_keywords, String<Dna>("taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
    setHost(pattern_dna, dna_keywords);

    clear(pos);
    while (find(finderText, pattern_dna)) {
        append(pos,position(finderText));
        SEQAN_ASSERT_EQ(position(finderText), beginPosition(finderText));
        SEQAN_ASSERT_EQ(endPosition(finderText), beginPosition(finderText) + length(finderText));
        SEQAN_ASSERT_EQ(length(finderText), length(dna_keywords[position(pattern_dna)]));
        SEQAN_ASSERT_EQ(begin(finderText), begin(text) + beginPosition(finderText));
        SEQAN_ASSERT_EQ(end(finderText), begin(text) + endPosition(finderText));
        SEQAN_ASSERT_EQ(infix(finderText), dna_keywords[position(pattern_dna)]);
    }

    SEQAN_ASSERT_EQ(pos[0], 0u);
    SEQAN_ASSERT_EQ(pos[1], 5u);
    SEQAN_ASSERT_EQ(pos[2], 10u);
    SEQAN_ASSERT_EQ(pos[3], 15u);
    SEQAN_ASSERT_EQ(pos[4], 20u);
    SEQAN_ASSERT_EQ(pos[5], 25u);
    SEQAN_ASSERT_EQ(length(pos), 6u);

    //____________________________________________________________________________
    // Test2 - Multiple keywords
    String<char> hst("annual_announce_any_annually");
    Finder<String<char> > fd(hst);

    typedef String<String<char> > TN;
    TN kyw;
    appendValue(kyw, String<char>("announce"));
    appendValue(kyw, String<char>("annual"));
    appendValue(kyw, String<char>("annually"));
    Pattern<TN, TAlgorithmSpec> pt(kyw);

    String<TPosition> finderPos;
    String<TPosition> keywordIndex;
    while (find(fd, pt)) {
        append(finderPos,position(fd));
        append(keywordIndex,position(pt));
        SEQAN_ASSERT_EQ(position(fd), beginPosition(fd));
        SEQAN_ASSERT_EQ(endPosition(fd), beginPosition(fd) + length(fd));
        SEQAN_ASSERT_EQ(length(fd), length(kyw[position(pt)]));
        SEQAN_ASSERT_EQ(begin(fd), begin(hst) + beginPosition(fd));
        SEQAN_ASSERT_EQ(end(fd), begin(hst) + endPosition(fd));
        SEQAN_ASSERT_EQ(infix(fd), kyw[position(pt)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 4u);
    SEQAN_ASSERT_EQ(length(keywordIndex), 4u);
    SEQAN_ASSERT_EQ(finderPos[0], 0u);
    SEQAN_ASSERT_EQ(keywordIndex[0], 1u);
    SEQAN_ASSERT_EQ(finderPos[1], 7u);
    SEQAN_ASSERT_EQ(keywordIndex[1], 0u);
    SEQAN_ASSERT_EQ(finderPos[2], 20u);
    SEQAN_ASSERT_EQ(keywordIndex[2], 1u);
    SEQAN_ASSERT_EQ(finderPos[3], 20u);
    SEQAN_ASSERT_EQ(keywordIndex[3], 2u);

    String<Dna> hstDna("AGATACGATATATAC");
    Finder<String<Dna> > fdDna(hstDna);

    typedef String<String<Dna> > TNDna;
    TNDna kywDna;
    appendValue(kywDna, String<Dna>("ATATATA"));
    appendValue(kywDna, String<Dna>("TATAT"));
    appendValue(kywDna, String<Dna>("ACGATAT"));
    Pattern<TNDna, TAlgorithmSpec> ptDna(kywDna);

    clear(finderPos);
    clear(keywordIndex);
    while (find(fdDna, ptDna)) {
        appendValue(finderPos, position(fdDna));
        appendValue(keywordIndex, position(ptDna));
        SEQAN_ASSERT_EQ(position(fdDna), beginPosition(fdDna));
        SEQAN_ASSERT_EQ(endPosition(fdDna), beginPosition(fdDna) + length(fdDna));
        SEQAN_ASSERT_EQ(length(fdDna), length(kywDna[position(ptDna)]));
        SEQAN_ASSERT_EQ(begin(fdDna), begin(hstDna) + beginPosition(fdDna));
        SEQAN_ASSERT_EQ(end(fdDna), begin(hstDna) + endPosition(fdDna));
        SEQAN_ASSERT_EQ(infix(fdDna), kywDna[position(ptDna)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 3u);
    SEQAN_ASSERT_EQ(length(keywordIndex), 3u);
    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 4u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 2u);
        SEQAN_ASSERT_EQ(finderPos[1], 7u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0u);
        SEQAN_ASSERT_EQ(finderPos[2], 8u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1u);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 4u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 2u);
        SEQAN_ASSERT_EQ(finderPos[1], 8u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1u);
        SEQAN_ASSERT_EQ(finderPos[2], 7u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0u);
    }
    //____________________________________________________________________________
    // Test2 - Multiple keywords that do not fit into a machine word
    String<Dna> my_haystack("AGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATAC");
    Finder<String<Dna> > my_finder(my_haystack);

    typedef String<String<Dna> > TNeedle_My;
    TNeedle_My my_keywords;
    appendValue(my_keywords, String<Dna>("ATATATA"));
    appendValue(my_keywords, String<Dna>("ACCGATCCAT"));
    appendValue(my_keywords, String<Dna>("TATAT"));
    appendValue(my_keywords, String<Dna>("ACCGAT"));
    appendValue(my_keywords, String<Dna>("ACGATAT"));
    appendValue(my_keywords, String<Dna>("CCAA"));
    Pattern<TNeedle_My, TAlgorithmSpec> my_pattern(my_keywords);

    clear(finderPos);
    clear(keywordIndex);
    while (find(my_finder, my_pattern)) {
        //std::cout << position(my_finder) << "-" << position(my_pattern) << std::endl;
        append(finderPos,position(my_finder));
        append(keywordIndex,position(my_pattern));
        SEQAN_ASSERT_EQ(position(my_finder), beginPosition(my_finder));
        SEQAN_ASSERT_EQ(endPosition(my_finder), beginPosition(my_finder) + length(my_finder));
        SEQAN_ASSERT_EQ(length(my_finder), length(my_keywords[position(my_pattern)]));
        SEQAN_ASSERT_EQ(begin(my_finder), begin(my_haystack) + beginPosition(my_finder));
        SEQAN_ASSERT_EQ(end(my_finder), begin(my_haystack) + endPosition(my_finder));
        SEQAN_ASSERT_EQ(infix(my_finder), my_keywords[position(my_pattern)]);
    }

    SEQAN_ASSERT_EQ(length(finderPos), 15u);
    SEQAN_ASSERT_EQ(length(keywordIndex), 15u);
    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 4u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 4u);
        SEQAN_ASSERT_EQ(finderPos[1], 7u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0u);
        SEQAN_ASSERT_EQ(finderPos[2], 8u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 2u);
        SEQAN_ASSERT_EQ(finderPos[3], 19u);
        SEQAN_ASSERT_EQ(keywordIndex[3], 4u);
        SEQAN_ASSERT_EQ(finderPos[4], 22u);
        SEQAN_ASSERT_EQ(keywordIndex[4], 0u);
        SEQAN_ASSERT_EQ(finderPos[5], 23u);
        SEQAN_ASSERT_EQ(keywordIndex[5], 2u);
        SEQAN_ASSERT_EQ(finderPos[6], 34u);
        SEQAN_ASSERT_EQ(keywordIndex[6], 4u);
        SEQAN_ASSERT_EQ(finderPos[7], 37u);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0u);
        SEQAN_ASSERT_EQ(finderPos[8], 38u);
        SEQAN_ASSERT_EQ(keywordIndex[8], 2u);
        SEQAN_ASSERT_EQ(finderPos[9], 49u);
        SEQAN_ASSERT_EQ(keywordIndex[9], 4u);
        SEQAN_ASSERT_EQ(finderPos[10], 52u);
        SEQAN_ASSERT_EQ(keywordIndex[10], 0u);
        SEQAN_ASSERT_EQ(finderPos[11], 53u);
        SEQAN_ASSERT_EQ(keywordIndex[11], 2u);
        SEQAN_ASSERT_EQ(finderPos[12], 64u);
        SEQAN_ASSERT_EQ(keywordIndex[12], 4u);
        SEQAN_ASSERT_EQ(finderPos[13], 67u);
        SEQAN_ASSERT_EQ(keywordIndex[13], 0u);
        SEQAN_ASSERT_EQ(finderPos[14], 68u);
        SEQAN_ASSERT_EQ(keywordIndex[14], 2u);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 4u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 4u);
        SEQAN_ASSERT_EQ(finderPos[1], 8u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 2u);
        SEQAN_ASSERT_EQ(finderPos[2], 7u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0u);
        SEQAN_ASSERT_EQ(finderPos[3], 19u);
        SEQAN_ASSERT_EQ(keywordIndex[3], 4u);
        SEQAN_ASSERT_EQ(finderPos[4], 23u);
        SEQAN_ASSERT_EQ(keywordIndex[4], 2u);
        SEQAN_ASSERT_EQ(finderPos[5], 22u);
        SEQAN_ASSERT_EQ(keywordIndex[5], 0u);
        SEQAN_ASSERT_EQ(finderPos[6], 34u);
        SEQAN_ASSERT_EQ(keywordIndex[6], 4u);
        SEQAN_ASSERT_EQ(finderPos[7], 38u);
        SEQAN_ASSERT_EQ(keywordIndex[7], 2u);
        SEQAN_ASSERT_EQ(finderPos[8], 37u);
        SEQAN_ASSERT_EQ(keywordIndex[8], 0u);
        SEQAN_ASSERT_EQ(finderPos[9], 49u);
        SEQAN_ASSERT_EQ(keywordIndex[9], 4u);
        SEQAN_ASSERT_EQ(finderPos[10], 53u);
        SEQAN_ASSERT_EQ(keywordIndex[10], 2u);
        SEQAN_ASSERT_EQ(finderPos[11], 52u);
        SEQAN_ASSERT_EQ(keywordIndex[11], 0u);
        SEQAN_ASSERT_EQ(finderPos[12], 64u);
        SEQAN_ASSERT_EQ(keywordIndex[12], 4u);
        SEQAN_ASSERT_EQ(finderPos[13], 68u);
        SEQAN_ASSERT_EQ(keywordIndex[13], 2u);
        SEQAN_ASSERT_EQ(finderPos[14], 67u);
        SEQAN_ASSERT_EQ(keywordIndex[14], 0u);
    }

    //____________________________________________________________________________
    // Multiple keywords with overlapping matches
    String<Dna> my2_haystack("aaaacaaa");
    Finder<String<Dna> > my2_finder(my2_haystack);

    typedef String<String<Dna> > TNeedle_My2;
    TNeedle_My2 my2_keywords;
    appendValue(my2_keywords, String<Dna>("aa"));
    appendValue(my2_keywords, String<Dna>("aaa"));
    appendValue(my2_keywords, String<Dna>("ac"));
    appendValue(my2_keywords, String<Dna>("aac"));
    appendValue(my2_keywords, String<Dna>("gctccacctgacctagcccatggggcccaaatttccggccttaattcccattt"));
    Pattern<TNeedle_My2, TAlgorithmSpec> my2_pattern(my2_keywords);

    clear(finderPos);
    clear(keywordIndex);
    while (find(my2_finder, my2_pattern)) {
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
        SEQAN_ASSERT_EQ(position(my2_finder), beginPosition(my2_finder));
        SEQAN_ASSERT_EQ(endPosition(my2_finder), beginPosition(my2_finder) + length(my2_finder));
        SEQAN_ASSERT_EQ(length(my2_finder), length(my2_keywords[position(my2_pattern)]));
        SEQAN_ASSERT_EQ(begin(my2_finder), begin(my2_haystack) + beginPosition(my2_finder));
        SEQAN_ASSERT_EQ(end(my2_finder), begin(my2_haystack) + endPosition(my2_finder));
        SEQAN_ASSERT_EQ(infix(my2_finder), my2_keywords[position(my2_pattern)]);
    }

    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 0u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0u);
        SEQAN_ASSERT_EQ(finderPos[1], 0u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1u);
        SEQAN_ASSERT_EQ(finderPos[2], 1u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0u);
        SEQAN_ASSERT_EQ(finderPos[3], 1u);
        SEQAN_ASSERT_EQ(keywordIndex[3], 1u);
        SEQAN_ASSERT_EQ(finderPos[4], 2u);
        SEQAN_ASSERT_EQ(keywordIndex[4], 0u);
        SEQAN_ASSERT_EQ(finderPos[5], 2u);
        SEQAN_ASSERT_EQ(keywordIndex[5], 3u);
        SEQAN_ASSERT_EQ(finderPos[6], 3u);
        SEQAN_ASSERT_EQ(keywordIndex[6], 2u);
        SEQAN_ASSERT_EQ(finderPos[7], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0u);
        SEQAN_ASSERT_EQ(finderPos[8], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[8], 1u);
        SEQAN_ASSERT_EQ(finderPos[9], 6u);
        SEQAN_ASSERT_EQ(keywordIndex[9], 0u);
    } else{
        SEQAN_ASSERT_EQ(finderPos[0], 0u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0u);
        SEQAN_ASSERT_EQ(finderPos[1], 1u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0u);
        SEQAN_ASSERT_EQ(finderPos[2], 0u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1u);
        SEQAN_ASSERT_EQ(finderPos[3], 2u);
        SEQAN_ASSERT_EQ(keywordIndex[3], 0u);
        SEQAN_ASSERT_EQ(finderPos[4], 1u);
        SEQAN_ASSERT_EQ(keywordIndex[4], 1u);
        SEQAN_ASSERT_EQ(finderPos[5], 3u);
        SEQAN_ASSERT_EQ(keywordIndex[5], 2u);
        SEQAN_ASSERT_EQ(finderPos[6], 2u);
        SEQAN_ASSERT_EQ(keywordIndex[6], 3u);
        SEQAN_ASSERT_EQ(finderPos[7], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[7], 0u);
        SEQAN_ASSERT_EQ(finderPos[8], 6u);
        SEQAN_ASSERT_EQ(keywordIndex[8], 0u);
        SEQAN_ASSERT_EQ(finderPos[9], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[9], 1u);
    }

    //____________________________________________________________________________
    // Multiple duplicated keywords with overlapping matches, jumping finder
    goBegin(my2_finder); // That's a repositioning
    clear(my2_finder);      // That's why, clear state
    clear(finderPos);
    clear(keywordIndex);

    unsigned int hits = 0;
    while (find(my2_finder, my2_pattern)) {
        if (hits < 2) break;
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
        ++hits;
    }
    goBegin(my2_finder);
    my2_finder+=5;
    clear(my2_finder);
    clear(finderPos);
    clear(keywordIndex);
    while (find(my2_finder, my2_pattern)) {
        //std::cout << position(my2_finder) << ":" << position(my2_pattern) << std::endl;
        append(finderPos,position(my2_finder));
        append(keywordIndex,position(my2_pattern));
    }

    if (order_by_begin_position) {
        SEQAN_ASSERT_EQ(finderPos[0], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0u);
        SEQAN_ASSERT_EQ(finderPos[1], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 1u);
        SEQAN_ASSERT_EQ(finderPos[2], 6u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 0u);
    } else {
        SEQAN_ASSERT_EQ(finderPos[0], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[0], 0u);
        SEQAN_ASSERT_EQ(finderPos[1], 6u);
        SEQAN_ASSERT_EQ(keywordIndex[1], 0u);
        SEQAN_ASSERT_EQ(finderPos[2], 5u);
        SEQAN_ASSERT_EQ(keywordIndex[2], 1u);
    }
}


template <typename TAlgorithmSpec>
void Test_OnlineAlgWildcards() {
	typedef typename Position<CharString>::Type TPosition;

	String<TPosition> pos;

    //____________________________________________________________________________
    // Test1 - simple find wo wildcards
    String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
    Finder<String<char> > finder(haystack);

    String<char> needle("ist");
    Pattern<String<char>, TAlgorithmSpec> pattern(needle);
    clear(pos);
    while (find(finder, pattern))
        appendValue(pos, position(finder));

    SEQAN_ASSERT_EQ(pos[0], 7u);
    SEQAN_ASSERT_EQ(pos[1], 33u);
    SEQAN_ASSERT_EQ(host(pattern), needle);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), needle);

    //____________________________________________________________________________
    // Test - validation of patterns
    needle = "ist";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), true);
    SEQAN_ASSERT_EQ(valid(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)), true);

    needle = "i[a-z]s{3,4}t?a*a+c..\\a";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), true);

    needle = "i[st";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist\\";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist?*";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,4}";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,a}";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    needle = "ist{5,";
    setHost(pattern, needle);
    SEQAN_ASSERT_EQ(valid(pattern), false);

    //____________________________________________________________________________
    // Test - searching with invalid needles
    haystack = "Dies i[st ein Haystack. Ja, das i[st wirklich einer!";
    setHost(finder, haystack);
    clear(finder);

    needle = "i[st";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(length(pos), 0u);

    //____________________________________________________________________________
    // Test - handle needles with wildcards
    // to produce two \ in the pattern you need to escape both of them
    needle = "aa+c*[a-z]xx?aa\\\\";
    SEQAN_ASSERT_EQ(_lengthWithoutWildcards(needle), 9u);

    //____________________________________________________________________________
    // Test - optional characters (?)
    haystack = "abc__ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab?c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 2u);
    SEQAN_ASSERT_EQ(pos[1], 6u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test - repeatable characters (+)
    haystack = "abc__abbbbbc_ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab+c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 2u);
    SEQAN_ASSERT_EQ(pos[1], 11u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test - repeatable characters (*)
    haystack = "abc__abbbbbc_ac";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab*c";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 2u);
    SEQAN_ASSERT_EQ(pos[1], 11u);
    SEQAN_ASSERT_EQ(pos[2], 14u);

    SEQAN_ASSERT_EQ(length(pos), 3u);

    //____________________________________________________________________________
    // Test - wildcard matching
    haystack = "acccdfabdeeef";
    setHost(finder, haystack);
    clear(finder);

    needle = "ab?c*de+f";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 12u);
    SEQAN_ASSERT_EQ(length(pos), 1u);

    //____________________________________________________________________________
    // Test - wildcard matching (hard case)
    haystack = "aacccdfacccdeeef";
    setHost(finder, haystack);
    clear(finder);

    needle = "a*c*de+f";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 15u);
    SEQAN_ASSERT_EQ(length(pos), 1u);


    //____________________________________________________________________________
    // Test - character classes matching
    haystack = "annual_Annual_znnual_Znnual";
    setHost(finder, haystack);
    clear(finder);

    needle = "[a-zA]nnual";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 5u);
    SEQAN_ASSERT_EQ(pos[1], 12u);
    SEQAN_ASSERT_EQ(pos[2], 19u);
    SEQAN_ASSERT_EQ(length(pos), 3u);

    //____________________________________________________________________________
    // Test - long needles
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstuvwxyzabcdefg";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 32u);
    SEQAN_ASSERT_EQ(pos[1], 58u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test - long needles with character classes
    //              abcdefghijklmnopqrstuvwxyzabcdefghijkl
    //                                        abcdefghijklmnopqrstuvwxyzabcdefghijkl
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuzwxyzabcdefghzjklaabcdefhijkl";
    setHost(finder, haystack);
    clear(finder);

    needle = "abcdefghijklmnopqrstu[vz]wxyzabcdefgh[iz]jkl";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 37u);
    SEQAN_ASSERT_EQ(pos[1], 63u);
    SEQAN_ASSERT_EQ(length(pos), 2u);


    //____________________________________________________________________________
    // Test - long needles with repeating characters
    //              abcdefghijklmnopqrstuvwxyzabcdefghijkl
    //                                                                                                        abcdefghijklmnopqrstuvwxyzabcdefghijkl
    haystack = "abcdefghijklmnopqrstuvwxyzabcdefghiiiiijkl____aaaaabcdefghijklmnopqrstuvwxyzabcdeghijkl__aaaaabcdefghijklmnopqrstuvwxyzabcdefghjkl";
    setHost(finder, haystack);
    clear(finder);

    needle = "aa*bcdefghijklmnopqrstuvwxyzabcdef?g?hi+jkl";
    setHost(pattern, needle);

    clear(pos);
    while (find(finder, pattern))
        append(pos,position(finder));

    SEQAN_ASSERT_EQ(pos[0], 41u);
    SEQAN_ASSERT_EQ(pos[1], 86u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

    //____________________________________________________________________________
    // Test - handle .
    haystack = "annual_Annual_znnual";
    setHost(finder, haystack);
    clear(finder);

    needle = ".nnual";
    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 5u);
    SEQAN_ASSERT_EQ(pos[1], 12u);
    SEQAN_ASSERT_EQ(pos[2], 19u);
    SEQAN_ASSERT_EQ(length(pos), 3u);

    //____________________________________________________________________________
    // Test - handle backslash
    haystack = "annual_Annual_.nnual";
    setHost(finder, haystack);
    clear(finder);

    needle = "\\.nnual";
    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)){
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 19u);
    SEQAN_ASSERT_EQ(length(pos), 1u);



    //____________________________________________________________________________
    // Test - handle bounded length repeats {n,m}
    haystack = "aannual_aaaannual_annual";
    setHost(finder, haystack);
    clear(finder);

    needle = "a{2,5}n{2}ual";

    SEQAN_ASSERT_EQ(_lengthWithoutWildcards(needle), 10u);

    setHost(pattern, needle);
    clear(pos);
    while (find(finder, pattern)) {
        append(pos,position(finder));
    }
    SEQAN_ASSERT_EQ(pos[0], 6u);
    SEQAN_ASSERT_EQ(pos[1], 16u);
    SEQAN_ASSERT_EQ(length(pos), 2u);


    //____________________________________________________________________________
    // Test - handle different types of Pattern and Needle
    String<Dna> dna_haystack("AAACCTATGGGTTTAAAACCCTGAAACCCC");
    Finder<String<Dna> > dna_finder(dna_haystack);

    String<char> char_needle("a{3}c+t[ag].");
    Pattern<String<Dna>, TAlgorithmSpec> dna_pattern(char_needle);
    clear(pos);

    while (find(dna_finder, dna_pattern))
        append(pos,position(dna_finder));

    SEQAN_ASSERT_EQ(pos[0], 7u);
    SEQAN_ASSERT_EQ(pos[1], 23u);
    SEQAN_ASSERT_EQ(length(pos), 2u);

}


template <typename TPatternSpec>
void Test_Approx_EditDist() {
    //test DPSearch
    String<char> hstk("any_annealing");
    String<char> nl("annual");

    Finder<String<char> > fd(hstk);

    Pattern<String<char>, TPatternSpec> pt(nl, -2);

    SEQAN_ASSERT(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 8u);
    SEQAN_ASSERT_EQ(getScore(pt), -2);
    SEQAN_ASSERT(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 9u);
    SEQAN_ASSERT_EQ(getScore(pt), -1);
    SEQAN_ASSERT(find(fd, pt));
    SEQAN_ASSERT_EQ(position(fd), 10u);
    SEQAN_ASSERT_EQ(getScore(pt), -2);

    SEQAN_ASSERT_NOT(find(fd,pt));

    String<char> haystk("Dies ist der Haystack des Tests. Ja, das ist er wirklich!");
    String<char> ndl("des");

    Finder<String<char> > fnd(haystk);

    Pattern<String<char>, TPatternSpec> pat(ndl, -2);
    SEQAN_ASSERT_EQ(host(pat), ndl);
    SEQAN_ASSERT_EQ(host(reinterpret_cast<Pattern<String<char>, TPatternSpec> const &>(pat)), ndl);

    SEQAN_ASSERT_EQ(scoreLimit(pat), -2);
    setScoreLimit(pat, -1);
    SEQAN_ASSERT_EQ(scoreLimit(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 3u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 10u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 11u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 23u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 24u);
    SEQAN_ASSERT_EQ(getScore(pat), 0);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 25u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 28u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT(find(fnd, pat));
    SEQAN_ASSERT_EQ(position(fnd), 39u);
    SEQAN_ASSERT_EQ(getScore(pat), -1);

    SEQAN_ASSERT_NOT(find(fnd, pat));

    // Test with long needles and a Dna Alphabet
    String<Dna> long_haystk("taaaataaaatacaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat");
    String<Dna> long_ndl("taaaataaaatacaataaaataaaatataataaaataaaataaaat");

    Finder<String<Dna> > long_fnd(long_haystk);

    Pattern<String<Dna>, TPatternSpec> long_pat(long_ndl, -2);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 44u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 45u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -1);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 46u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 60u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 65u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT(find(long_fnd,long_pat));
    SEQAN_ASSERT_EQ(position(long_fnd), 70u);
    SEQAN_ASSERT_EQ(getScore(long_pat), -2);

    SEQAN_ASSERT_NOT(find(long_fnd,long_pat));

    //____________________________________________________________________________

    String<char> haystack_1 = "123XXXabaXXX45aba123";
    String<char> needle_1 = "XXXaba";
    Finder<String<char> > finder_1(haystack_1);
    Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 7u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXa");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 8u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXab");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 9u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), 0);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "Xaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), 0);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "23XXXaba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 10u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -1);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "3XXXabaX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 11u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXXabaXX");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 15u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXX45a");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(endPosition(finder_1), 17u);
    SEQAN_ASSERT_EQ(getScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "X45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XX45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT(findBegin(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(infix(finder_1), "XXX45aba");
    SEQAN_ASSERT_EQ(getBeginScore(pattern_1), -2);
    SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));

    SEQAN_ASSERT_NOT(find(finder_1, pattern_1));
}


// Test prefix search.
template <typename TPatternSpec>
void Test_Approx_Prefix_EditDist() {
    String<char> haystack_1 = "mississippi";
    String<char> needle_1 = "misssi";
    Finder<String<char> > finder_1(haystack_1);
    Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);
    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 3u);
    SEQAN_ASSERT_EQ(length(finder_1), 4u);
    SEQAN_ASSERT_EQ(beginPosition(finder_1), 0u);
    SEQAN_ASSERT_EQ(endPosition(finder_1), 4u);
    SEQAN_ASSERT_EQ(infix(finder_1), "miss");
    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 4u);
    SEQAN_ASSERT_EQ(infix(finder_1), "missi");
    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 5u);
    SEQAN_ASSERT_EQ(infix(finder_1), "missis");
    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 6u);
    SEQAN_ASSERT(find(finder_1, pattern_1));
    SEQAN_ASSERT_EQ(position(finder_1), 7u);
    SEQAN_ASSERT_NOT(find(finder_1, pattern_1));


    String<char> haystack_2 = "yyyXXaba";
    String<char> needle_2 = "yyyaba";
    Finder<String<char> > finder_2(haystack_2);
    Pattern<String<char>, TPatternSpec> pattern_2(needle_2, -2);
    SEQAN_ASSERT(find(finder_2, pattern_2));
    SEQAN_ASSERT_EQ(position(finder_2), 5u);
    SEQAN_ASSERT_EQ(infix(finder_2), "yyyXXa");
    SEQAN_ASSERT(find(finder_2, pattern_2));
    SEQAN_ASSERT_EQ(position(finder_2), 7u);
    SEQAN_ASSERT_EQ(infix(finder_2), "yyyXXaba");
    SEQAN_ASSERT_NOT(find(finder_2, pattern_2));


    String<char> haystack_3 = "testtexttext";
    String<char> needle_3 = "mismatch";
    Finder<String<char> > finder_3(haystack_3);
    Pattern<String<char>, TPatternSpec> pattern_3(needle_3, -2);
    SEQAN_ASSERT_NOT(find(finder_3, pattern_3));


    String<char> haystack_4 = "testtext";
    String<char> needle_4 = "a longer mismatch";
    Finder<String<char> > finder_4(haystack_4);
    Pattern<String<char>, TPatternSpec> pattern_4(needle_4, -2);
    SEQAN_ASSERT_NOT(find(finder_4, pattern_4));


    String<char> haystack_5 = "exactmatching";
    String<char> needle_5 = "exact";
    Finder<String<char> > finder_5(haystack_5);
    Pattern<String<char>, TPatternSpec> pattern_5(needle_5, 0);
    SEQAN_ASSERT(find(finder_5, pattern_5));
    SEQAN_ASSERT_EQ(position(finder_5), 4u);
    SEQAN_ASSERT_NOT(find(finder_5, pattern_5));


    String<char> haystack_6 = "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX";
    String<char> needle_6 =   "this is a text that is a bit longer than one machine word of 32 or 64 bits. XYX";
    Finder<String<char> > finder_6(haystack_6);
    Pattern<String<char>, TPatternSpec> pattern_6(needle_6, -2);
    SEQAN_ASSERT(find(finder_6, pattern_6));
    SEQAN_ASSERT_EQ(infix(finder_6), "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAX");
    SEQAN_ASSERT(find(finder_6, pattern_6));
    SEQAN_ASSERT_EQ(infix(finder_6), "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX");
    SEQAN_ASSERT_NOT(find(finder_6, pattern_6));

    // The very high score limit with short needle all-mismatching the
    // pattern.  This is used to mirror the test below that use a long
    // needle.
    {
        const char *kHaystackStr = "CTCT";
        const char *kNeedleStr =   "AAAA";
        const int kNeedleLen = strlen(kNeedleStr);
        const int kScoreLimit = -1000;
        // Haystack has a length of 100, needle of 80.
        String<char> haystack(kHaystackStr);
        String<char> needle(kNeedleStr);
        Finder<String<char> > finder(haystack);
        Pattern<String<char>, TPatternSpec> pattern(needle, kScoreLimit);

        while (find(finder, pattern)) {
            SEQAN_ASSERT_EQ(-kNeedleLen, getScore(pattern));
        }
    }

    // Test very high score limit with large needles and very high
    // maximum scores.
    {
        const char *kHaystackStr = "CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT";
        const char *kNeedleStr =   "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        const int kNeedleLen = strlen(kNeedleStr);
        const int kScoreLimit = -1000;
        // Haystack has a length of 100, needle of 80.
        String<char> haystack(kHaystackStr);
        String<char> needle(kNeedleStr);
        Finder<String<char> > finder(haystack);
        Pattern<String<char>, TPatternSpec> pattern(needle, kScoreLimit);

        while (find(finder, pattern)) {
            SEQAN_ASSERT_EQ(-kNeedleLen, getScore(pattern));
        }
    }
}


SEQAN_DEFINE_TEST(test_find_online_Simple) {
    Test_OnlineAlg<Simple>();
}


SEQAN_DEFINE_TEST(test_find_online_Horspool) {
    Test_OnlineAlg<Horspool>();
}


SEQAN_DEFINE_TEST(test_find_online_ShiftAnd) {
    Test_OnlineAlg<ShiftAnd>();
}


SEQAN_DEFINE_TEST(test_find_online_ShiftOr) {
    Test_OnlineAlg<ShiftOr>();
}


SEQAN_DEFINE_TEST(test_find_online_BndmAlgo) {
    Test_OnlineAlg<BndmAlgo>();
}


SEQAN_DEFINE_TEST(test_find_online_BFAM_Oracle) {
    Test_OnlineAlg<Bfam<Oracle> >();
}


SEQAN_DEFINE_TEST(test_find_online_BFAM_Trie) {
    Test_OnlineAlg<Bfam<Trie> >();
}


SEQAN_DEFINE_TEST(test_find_online_wildcards) {
    Test_OnlineAlgWildcards<WildShiftAnd>();
}


SEQAN_DEFINE_TEST(test_find_online_multi_AhoCorasick) {
    Test_OnlineAlgMulti<AhoCorasick>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultipleShiftAnd) {
    // TODO(holtgrew): Original comment: "leaks".
    // TODO(holtgrew): Fails, but was commented out in original code.
    // Test_OnlineAlgMulti<MultipleShiftAnd>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_SetHorspool) {
    Test_OnlineAlgMulti<SetHorspool>(false);
}


SEQAN_DEFINE_TEST(test_find_online_multi_WuManber) {
    Test_OnlineAlgMulti<WuManber>(true);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultiBFAM_Oracle) {
    Test_OnlineAlgMulti<MultiBfam<Oracle> >(true);
}


SEQAN_DEFINE_TEST(test_find_online_multi_MultiBFAM_Trie) {
    Test_OnlineAlgMulti<MultiBfam<Trie> >(true);
}


SEQAN_DEFINE_TEST(test_find_approx_prefix_edit_dist_dpsearch) {
    Test_Approx_Prefix_EditDist<DPSearch<Score<>, FindPrefix> >();
}


SEQAN_DEFINE_TEST(test_approx_prefix_edit_dist_myers) {
    Test_Approx_Prefix_EditDist<Myers<FindPrefix> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_dp_search_simple_score) {
    Test_Approx_EditDist<DPSearch<SimpleScore> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dp_search_simple_score_legacy_case) {
    // TODO(holtgrew): This was written out like this in Test_Approx() in original code.
    // Test DPSearch.
    Pattern<String<char>, DPSearch<SimpleScore> > pat1;

    SimpleScore sc;
    setScoreGap(sc, -10);
    setScoringScheme(pat1, sc);
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -10);

    setScoreGap(sc, -1);
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -10);
    setScoringScheme(pat1, sc);
    SEQAN_ASSERT_EQ(scoreGap(scoringScheme(pat1)), -1);
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_myers) {
    Test_Approx_EditDist<Myers<> >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_abndm_algo) {
    Test_Approx_EditDist<AbndmAlgo >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical) {
    Test_Approx_EditDist<PexNonHierarchical>();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_hierarchical) {
    Test_Approx_EditDist<PexHierarchical>();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical_aho_corasick) {
    Test_Approx_EditDist< Pex<NonHierarchical,AhoCorasick > >();
}


SEQAN_DEFINE_TEST(test_approx_edit_dist_pex_non_hierarchical_multi_bfam) {
    Test_Approx_EditDist< Pex<NonHierarchical,MultiBfam<> > >();
}


// Test that shows the problem with Myers-Ukkonen approximate search
// when large needles were used and the maximum score limit was set
// higher than needle length.
SEQAN_DEFINE_TEST(test_regression_rmbench) {
    // The data to test with:  Needle, haystack and score limit.
    const char *kCharNeedle = "CCATATGCTTGTGTCGCGGGTTTATTTGCATTCGACCCAGTTGACTCGGAAGTCGAAATGTTCCTGCCCCGTTTCTGCGTTCCGTGCAGTTGCGCGGTCTGGTTGGGCGGGTCCCCCCCTGA";
    const char *kCharHaystack = "ATTCCATATGCTTGTGTCGCGGGTTTATTTGCATTCGACCCAGTTGACTCGGAAGTCGAAATGTTCCTGCCCCGTTTCTGCGTTCCGTGCAGTTGCGCGGTCTGGTTGGGCGGGTCCCCCGCCTACGGATGCACGTTCTCCCGGGCTCGTAAATCC";
    const int kScoreLimit = -100000;

    // Build actual DNA needle and haystack, the pattern and the
    // finder.
    DnaString needle(kCharNeedle);
    DnaString haystack(kCharHaystack);
    Finder<DnaString> finder(haystack);
    Pattern<DnaString, MyersUkkonen> pattern(needle);
    setScoreLimit(pattern, kScoreLimit);

    // The following invariant should always hold: The pattern score
    // should be smaller than the needle size.
    while (find(finder, pattern)) {
        SEQAN_ASSERT_LEQ(pattern.errors, pattern.needleSize);
    }
    SEQAN_ASSERT_LEQ(pattern.errors, pattern.needleSize);
}


// Tests for the Simple Hamming Finder.
SEQAN_DEFINE_TEST(test_find_hamming_simple) {
    // Test that the interface works.
    {
        DnaString haystack("AC");
        DnaString needle("AA");
        // Define finder and pattern.
        Finder<DnaString> finder(haystack);
        Pattern<DnaString, HammingSimple> pattern(needle, -1);

        // Perform a search, run all functions defined on the pattern.
        bool res = find(finder, pattern);
        SEQAN_ASSERT(res);
        SEQAN_ASSERT_EQ(0u, position(finder));
        SEQAN_ASSERT_EQ(2u, endPosition(finder));
        SEQAN_ASSERT_EQ(-1, score(pattern));
        SEQAN_ASSERT_EQ(-1, getScore(pattern));
    }

    // Test for distance 0;
    {
        // TODO(holtgrew): Should be const, but finder does not allow this.
        // Define haystack and needle.
        DnaString haystack("ACA");
        DnaString needle("AA");
        // Define finder and pattern.
        Finder<DnaString> finder(haystack);
        Pattern<DnaString, HammingSimple> pattern(needle, 0);
        // Perform the searches;
        bool res;

        res = find(finder, pattern);
        SEQAN_ASSERT_NOT(res);
    }

    // Test for distance -1, -2, -3.  Should yield the same results,
    // as tested for below.  The numbers are interesting since they
    // are the first > 0, length of pattern, first greater than length
    // of pattern.
    for (int i = 1; i < 3; ++i) {
        // TODO(holtgrew): Should be const, but finder does not allow this.
        // Define haystack and needle.
        DnaString haystack("ACA");
        DnaString needle("AA");
        // Define finder and pattern.
        Finder<DnaString> finder(haystack);
        Pattern<DnaString, HammingSimple> pattern(needle, -i);
        // Perform the searches;
        bool res;

        res = find(finder, pattern);
        SEQAN_ASSERT(res);
        SEQAN_ASSERT_EQ(0u, position(finder));
        SEQAN_ASSERT_EQ(2u, endPosition(finder));
        SEQAN_ASSERT_EQ(-1, score(pattern));

        res = find(finder, pattern);
        SEQAN_ASSERT(res);
        SEQAN_ASSERT_EQ(1u, position(finder));
        SEQAN_ASSERT_EQ(3u, endPosition(finder));
        SEQAN_ASSERT_EQ(-1, score(pattern));
    }

    // Test setting the score limit.
    {
        // TODO(holtgrew): Should be const, but finder does not allow this.
        // Define haystack and needle.
        DnaString haystack("AAC");
        DnaString needle("AA");
        // Define finder and pattern.
        Finder<DnaString> finder(haystack);
        Pattern<DnaString, HammingSimple> pattern(needle, 0);
        // Perform the searches;
        bool res;

        res = find(finder, pattern);
        SEQAN_ASSERT(res);
        SEQAN_ASSERT_EQ(0u, position(finder));

        setScoreLimit(pattern, -1);

        res = find(finder, pattern);
        SEQAN_ASSERT(res);
        SEQAN_ASSERT_EQ(1u, position(finder));
    }
}


// Tests for a regression found in read mapper benchmark.
SEQAN_DEFINE_TEST(test_find_hamming_simple_regression_rmbench) {
    // TODO(holtgrew): Should be const, but finder does not allow this.
    // Define haystack and needle.
    DnaString haystack("CCCCCCCCCCCCCCCCCCCCA");
    DnaString needle("CCCCCCCCCCCCCCCCCCCC");
    // Define finder and pattern.
    Finder<DnaString> finder(haystack);
    Pattern<DnaString, HammingSimple> pattern(needle);
    setScoreLimit(pattern, -1000);
    // Perform the searches;
    bool res;
    res = find(finder, pattern);
    SEQAN_ASSERT(res);
    SEQAN_ASSERT_EQ(0u, position(finder));
    SEQAN_ASSERT_EQ(length(needle), endPosition(finder));
    SEQAN_ASSERT_EQ(0, getScore(pattern));

    res = find(finder, pattern);
    SEQAN_ASSERT(res);
    SEQAN_ASSERT_EQ(1u, position(finder));
    SEQAN_ASSERT_EQ(length(needle) + 1, endPosition(finder));
    SEQAN_ASSERT_EQ(-1, getScore(pattern));
}

/*
SEQAN_DEFINE_TEST(test_find_hamming_simple) {
    testFindApproximateHamming<HammingSimple, Dna>();
    }*/

SEQAN_DEFINE_TEST(test_myers_find_infix_find_begin_at_start) {
    String<char> haystack = "___AAA___AAA";
    String<char> needle = "___AAA";
    Finder<String<char> > finder(haystack);
    Pattern<String<char>, Myers<FindInfix> > pattern(needle, -2);

    // Find match: ___AAA___AAA
    //             ___A
    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(4u, endPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));  // TODO(holtgrew): getScore(pattern) is in book but should not be necessary
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));

    // Find match: ___AAA___AAA
    //             ___AA
    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(5u, endPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));  // TODO(holtgrew): getScore(pattern) is in book but should not be necessary
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));

    // Find match: ___AAA___AAA
    //             ___AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));  // TODO(holtgrew): getScore(pattern) is in book but should not be necessary
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(0u, beginPosition(finder));
}


SEQAN_DEFINE_TEST(test_myers_find_infix_find_begin_within) {
    String<char> haystack = "A___AAA___AAA";
    String<char> needle = "___AAA";
    Finder<String<char> > finder(haystack);
    Pattern<String<char>, Myers<FindInfix> > pattern(needle, -1);

    // Find match: A___AAA___AAA
    //              ___AA
    bool ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(6u, endPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));  // getScore(pattern) is in book but should not be necessary
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));

    // Find match: A___AAA___AAA
    //              ___AAA
    ret = find(finder, pattern);
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(7u, endPosition(finder));
    ret = findBegin(finder, pattern, getScore(pattern));  // getScore(pattern) is in book but should not be necessary
    SEQAN_ASSERT(ret);
    SEQAN_ASSERT_EQ(1u, beginPosition(finder));
}


template <typename TString, typename TSegmentOrString>
void test_find_on_segments_Helper(TString &haystack, const TSegmentOrString &needle) {
    Finder<TString> finder(haystack);
    Pattern<TSegmentOrString, Myers<FindInfix> > pattern(needle);

    bool didFind = false;
    while (find(finder, pattern)) {
        findBegin(finder, pattern);
        didFind = true;
    }
    SEQAN_ASSERT(didFind);

    // TODO(holtgrew): Some kind of assertion on the results.
}


// Test string search code on segments.  At the moment, this only
// tests whether the code compiles and runs through without any
// obvious errors.  No checks are done on the results.
SEQAN_DEFINE_TEST(test_find_on_segments) {
    // TODO(holtgrew): Should be const.
    CharString kHaystack = "CGATCGAT";
    CharString kNeedle = "GATC";

    test_find_on_segments_Helper<>(kHaystack, kNeedle);

    Segment<CharString, PrefixSegment> myPrefix(prefix(kNeedle, 1));
    test_find_on_segments_Helper<>(kHaystack, myPrefix);
    test_find_on_segments_Helper<>(kHaystack, prefix(kNeedle, 1));

    Segment<CharString, InfixSegment> myInfix(infix(kNeedle, 1, 2));
    test_find_on_segments_Helper<>(kHaystack, myInfix);
    test_find_on_segments_Helper<>(kHaystack, infix(kNeedle, 1, 2));

    Segment<CharString, SuffixSegment> mySuffix(suffix(kNeedle, 1));
    test_find_on_segments_Helper<>(kHaystack, mySuffix);
    test_find_on_segments_Helper<>(kHaystack, suffix(kNeedle, 1));
}


SEQAN_DEFINE_TEST(test_myers_trigger_bug) {
    DnaString haystackString = "TCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCCCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTC";
    DnaString needleString = "GAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAAGAAG";
    reverseComplement(needleString);
    typedef Segment<DnaString, InfixSegment> TSegment;
    typedef ModifiedString<TSegment, ModReverse> TSegmentRev;
//     std::cout << "haystack = " << haystackString << std::endl;
//     std::cout << "needle = " << needleString << std::endl;
//     std::cout << "haystack length = " << length(haystackString) << std::endl
//               << "needle length = " << length(needleString) << std::endl
//               << "machine word length = " << sizeof(unsigned) * 8 << std::endl;

    for (int maxDistance = 0; maxDistance < 10; ++maxDistance) {
        for (unsigned beginPosition = 0; beginPosition < length(haystackString) - length(needleString); ++beginPosition) {
//             std::cout << "max distance = " << maxDistance << ", begin position = " << beginPosition << std::endl;
            TSegmentRev haystackSegment(infix(haystackString, beginPosition, length(haystackString)));
//             std::cout << "haystack = " << haystackSegment << std::endl;
//             std::cout << "needle = " << needleString << std::endl;
//             std::cout << "-----------------------" << std::endl;

            String<size_t> positionsMyers;
            {
                TSegmentRev segment(infix(haystackString, beginPosition, length(haystackString)));
                Finder<TSegmentRev> finder(segment);
                Pattern<DnaString, MyersUkkonenGlobal > pattern(needleString, -maxDistance);
                while (find(finder, pattern)) {
                    appendValue(positionsMyers, endPosition(finder));
                }
            }

            String<size_t> positionsDpSearch;
            {
                TSegmentRev segment(infix(haystackString, beginPosition, length(haystackString)));
                Finder<TSegmentRev> finder(segment);
                Pattern<DnaString, DPSearch<EditDistanceScore, FindPrefix> > pattern(needleString, -maxDistance);
                while (find(finder, pattern)) {
                    appendValue(positionsDpSearch, endPosition(finder));
                }
            }

//             std::cout << "Myers end positions: " << std::endl;
//             for (size_t i = 0; i < length(positionsMyers); ++i)
//                 std::cout << positionsMyers[i] << " ";
//             std::cout << std::endl;
//             std::cout << "DPSearch end positions: " << std::endl;
//             for (size_t i = 0; i < length(positionsDpSearch); ++i)
//                 std::cout << positionsDpSearch[i] << " ";
//             std::cout << std::endl;

            SEQAN_ASSERT_EQ(length(positionsDpSearch), length(positionsMyers));
            for (unsigned i = 0; i < length(positionsMyers); ++i)
                SEQAN_ASSERT_EQ_MSG(positionsDpSearch[i], positionsMyers[i], "i = %u", i);
        }
    }
}


// Test myers algorithm: Palindrom vs. non-palindrom.
SEQAN_DEFINE_TEST(test_myers_find_begin) {
    {
        String<char> haystack_1 = "AABBAA";
        String<char> needle_1 = "ABBA";
        Finder<String<char> > finder_1(haystack_1);
        Pattern<String<char>, Myers<FindInfix> > pattern_1(needle_1, 0);
        SEQAN_ASSERT(find(finder_1, pattern_1));
        SEQAN_ASSERT_EQ(endPosition(finder_1), 5u);
        SEQAN_ASSERT_EQ(getScore(pattern_1), 0);
        SEQAN_ASSERT(findBegin(finder_1, pattern_1));
        SEQAN_ASSERT_EQ(infix(finder_1), "ABBA");
        SEQAN_ASSERT_EQ(getBeginScore(pattern_1), 0);
        SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));
    }
    {
        String<char> haystack_1 = "ABCD";
        String<char> needle_1 = "BCD";
        Finder<String<char> > finder_1(haystack_1);
        Pattern<String<char>, Myers<FindInfix> > pattern_1(needle_1, 0);
        SEQAN_ASSERT(find(finder_1, pattern_1));
        SEQAN_ASSERT_EQ(endPosition(finder_1), 4u);
        SEQAN_ASSERT_EQ(getScore(pattern_1), 0);
        SEQAN_ASSERT(findBegin(finder_1, pattern_1));
        SEQAN_ASSERT_EQ(infix(finder_1), "BCD");
        SEQAN_ASSERT_EQ(getBeginScore(pattern_1), 0);
        SEQAN_ASSERT_NOT(findBegin(finder_1, pattern_1));
    }
}

template <typename TPatternSpec>
void test_pattern_copycon() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2(p1);
    TPattern const p3(p2);
    TPattern const p4(p3);
}

template <typename TPatternSpec>
void test_pattern_assign() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2;
    TPattern const p3(p1);
    TPattern p4;
    p2 = p1;
    p4 = p3;
}

template <typename TPatternSpec>
void test_pattern_movecon() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2(std::move(p1));
    TPattern const p3(std::move(p2));
    TPattern const p4(std::move(p3));
}

template <typename TPatternSpec>
void test_pattern_moveassign() {
    typedef Pattern<CharString, TPatternSpec> TPattern;
    TPattern p1("Some needle");
    TPattern p2;
    TPattern const p3(p1);
    TPattern p4;
    p2 = std::move(p1);
    p4 = std::move(p3);
}

template <typename TPatternSpec>
void test_pattern_set_host()
{
    typedef Pattern<DnaString, TPatternSpec> TPattern;

    {  // Set to lvalue reference.
        TPattern p;
        DnaString ndl = "AGTGGATAGAGAT";
        setHost(p, ndl);

        SEQAN_ASSERT(&host(p) == &ndl);
    }

    {  // Set to lvalue with different alphabet causing the source to be changed.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT_EQ(&host(p), &ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT(&host(p) == &ndl);
        SEQAN_ASSERT_EQ(ndl, "ATACAAAA");
        SEQAN_ASSERT_EQ(ndl2, "AT CARLS");
    }

    {  // Set to const lvalue reference with different alphabet causing the source to be changed.
        DnaString const ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(ndl, "AGTGGATAGAGAT");
        SEQAN_ASSERT_EQ(ndl2, "AT CARLS");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set with rvalue reference.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, std::move(ndl));
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), "AGTGGATAGAGAT");
    }

    {  // Set to rvalue reference after setting holder dependent.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, ndl);
        SEQAN_ASSERT(&host(p) == &ndl);

        CharString ndl2 = "AT CARLS";
        setHost(p, std::move(ndl2));
        SEQAN_ASSERT(&host(p) == &ndl);
        SEQAN_ASSERT_EQ(ndl, "ATACAAAA");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set to rvalue reference before setting holder dependent.
        DnaString ndl = "AGTGGATAGAGAT";
        TPattern p;
        setHost(p, std::move(ndl));
        SEQAN_ASSERT(&host(p) != &ndl);
        SEQAN_ASSERT_EQ(host(p), "AGTGGATAGAGAT");

        CharString ndl2 = "AT CARLS";
        setHost(p, ndl2);
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }

    {  // Set to rvalue with different alphabet.
        TPattern p;
        setHost(p, "AT CARLS");
        SEQAN_ASSERT_EQ(host(p), "ATACAAAA");
    }
}

SEQAN_DEFINE_TEST(test_pattern_copycon) {
    // Test whether the needle is preserved in copying a pattern.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/318
    test_pattern_copycon<Simple>();
    test_pattern_copycon<Horspool>();
    test_pattern_copycon<ShiftAnd>();
    test_pattern_copycon<ShiftOr>();
    test_pattern_copycon<HammingSimple>();
    test_pattern_copycon<WildShiftAnd>();
    test_pattern_copycon<Bfam<Oracle> >();
    test_pattern_copycon<Bfam<Trie> >();
}

SEQAN_DEFINE_TEST(test_pattern_assign) {
    // Test whether the needle is preserved in assigning a pattern.
    // See http://trac.mi.fu-berlin.de/seqan/ticket/318
    test_pattern_assign<Simple>();
    test_pattern_assign<Horspool>();
    test_pattern_assign<ShiftAnd>();
    test_pattern_assign<ShiftOr>();
    test_pattern_assign<HammingSimple>();
    test_pattern_assign<WildShiftAnd>();
    test_pattern_assign<Bfam<Oracle> >();
    test_pattern_assign<Bfam<Trie> >();
}

SEQAN_DEFINE_TEST(test_pattern_movecon) {
    test_pattern_movecon<Simple>();
    test_pattern_movecon<Horspool>();
    test_pattern_movecon<ShiftAnd>();
    test_pattern_movecon<ShiftOr>();
    test_pattern_copycon<HammingSimple>();
    test_pattern_copycon<WildShiftAnd>();
    test_pattern_copycon<Bfam<Oracle> >();
    test_pattern_copycon<Bfam<Trie> >();
}

SEQAN_DEFINE_TEST(test_pattern_moveassign) {
    test_pattern_moveassign<Simple>();
    test_pattern_moveassign<Horspool>();
    test_pattern_moveassign<ShiftAnd>();
    test_pattern_moveassign<ShiftOr>();
    test_pattern_assign<HammingSimple>();
    test_pattern_assign<WildShiftAnd>();
    test_pattern_assign<Bfam<Oracle> >();
    test_pattern_assign<Bfam<Trie> >();
}

// TODO(rrahn): Should be a typed test for all pattern classes.
SEQAN_DEFINE_TEST(test_pattern_set_host) {
    test_pattern_set_host<Simple>();
    test_pattern_set_host<Horspool>();
    test_pattern_set_host<ShiftAnd>();
    test_pattern_set_host<ShiftOr>();
    test_pattern_set_host<HammingSimple>();
    test_pattern_set_host<WildShiftAnd>();
    test_pattern_set_host<Bfam<Oracle> >();
    test_pattern_set_host<Bfam<Trie> >();
    test_pattern_set_host<Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, void> >();
    test_pattern_set_host<Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, Myers<FindPrefix> > >();
    test_pattern_set_host<Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, void> >();
    test_pattern_set_host<Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, Myers<FindPrefix> > >();
    test_pattern_set_host<Myers<FindInfix, void> >();
    test_pattern_set_host<Myers<FindInfix, Myers<FindPrefix> > >();
    test_pattern_set_host<Myers<FindPrefix, void> >();
    test_pattern_set_host<Myers<FindPrefix, Myers<FindPrefix> > >();
}

SEQAN_BEGIN_TESTSUITE(test_find) {
//     SEQAN_CALL_TEST(test_myers_trigger_bug);
    SEQAN_CALL_TEST(test_myers_find_begin);
    SEQAN_CALL_TEST(test_myers_find_banded);
    SEQAN_CALL_TEST(test_myers_find_banded_csp);

    // Testing Myers<FindInfix> with findBegin().
    SEQAN_CALL_TEST(test_myers_find_infix_find_begin_at_start);
    SEQAN_CALL_TEST(test_myers_find_infix_find_begin_within);

    SEQAN_CALL_TEST(test_find_on_segments);

    SEQAN_CALL_TEST(test_find_hamming_simple);

    SEQAN_CALL_TEST(test_find_hamming_simple_regression_rmbench);

    // Testing MyersUkkonen with large needle and manual score limit.
    SEQAN_CALL_TEST(test_regression_rmbench);

    // Call all tests.
    SEQAN_CALL_TEST(test_find_online_Simple);
    SEQAN_CALL_TEST(test_find_online_Horspool);
    SEQAN_CALL_TEST(test_find_online_ShiftAnd);
    SEQAN_CALL_TEST(test_find_online_ShiftOr);
    SEQAN_CALL_TEST(test_find_online_BndmAlgo);
    SEQAN_CALL_TEST(test_find_online_BFAM_Oracle);
    SEQAN_CALL_TEST(test_find_online_BFAM_Trie);
    SEQAN_CALL_TEST(test_find_online_wildcards);
    SEQAN_CALL_TEST(test_find_online_multi_AhoCorasick);
    SEQAN_CALL_TEST(test_find_online_multi_MultipleShiftAnd);
    SEQAN_CALL_TEST(test_find_online_multi_SetHorspool);
    SEQAN_CALL_TEST(test_find_online_multi_WuManber);
    SEQAN_CALL_TEST(test_find_online_multi_MultiBFAM_Oracle);
    SEQAN_CALL_TEST(test_find_online_multi_MultiBFAM_Trie);
    SEQAN_CALL_TEST(test_find_approx_prefix_edit_dist_dpsearch);
    SEQAN_CALL_TEST(test_approx_prefix_edit_dist_myers);
    SEQAN_CALL_TEST(test_approx_edit_dist_dp_search_simple_score);
    SEQAN_CALL_TEST(test_approx_edit_dist_myers);
    // Test DP Serach.
    SEQAN_CALL_TEST(test_approx_edit_dp_search_simple_score_legacy_case);
    // Test other approximate search algorithms.
    SEQAN_CALL_TEST(test_approx_edit_dist_abndm_algo);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_hierarchical);
    // Tests with different multifinder.
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical_aho_corasick);
    SEQAN_CALL_TEST(test_approx_edit_dist_pex_non_hierarchical_multi_bfam);
    // Test for hamming distance approximate matching.
    SEQAN_CALL_TEST(test_find_hamming_simple);

    SEQAN_CALL_TEST(test_pattern_copycon);
    SEQAN_CALL_TEST(test_pattern_assign);

    SEQAN_CALL_TEST(test_pattern_movecon);
    SEQAN_CALL_TEST(test_pattern_moveassign);
    SEQAN_CALL_TEST(test_pattern_set_host);
}
SEQAN_END_TESTSUITE
