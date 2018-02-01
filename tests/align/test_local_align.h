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

#ifndef TESTS_ALIGN_TEST_LOCAL_ALIGN_H_
#define TESTS_ALIGN_TEST_LOCAL_ALIGN_H_

#include <iostream>
#include <cstdio>
#include <vector>

//#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/align.h>


using namespace std;
using namespace seqan;


SEQAN_DEFINE_TEST(testLocalAlign) {
    //align two sequences using Smith-Waterman-algorithm
    String<char> str0 = "ataagcgtctcg";
    String<char> str1 = "tcatagagttgc";

    Align< String<char>, ArrayGaps> ali;
    resize(rows(ali), 2);
    setSource(row(ali, 0), str0);
    setSource(row(ali, 1), str1);

    Score<int> score_type = Score<int>(2,-1,-2,0) ;
    LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>();

    int cutoff = 0;

    int score = _smithWaterman(ali,sw_finder,score_type,cutoff);

    SEQAN_ASSERT_EQ(score, 9);
    SEQAN_ASSERT(row(ali,0) == "ataagcgt");
    SEQAN_ASSERT(row(ali,1) == "ata-gagt");

    int i = 1;
    while (true){

        score = _smithWatermanGetNext(ali,sw_finder,score_type,cutoff);

        if(score==0){
    //        cout <<"No more alignments satisfying score > "<<cutoff<<"found.\n";
            break;
        }
        if(i == 1){
            SEQAN_ASSERT_EQ(score, 5);
            SEQAN_ASSERT(row(ali,0) == "tc-tcg");
            SEQAN_ASSERT(row(ali,1) == "tcatag");
        }
        if(i == 2){
            SEQAN_ASSERT_EQ(score, 4);
            SEQAN_ASSERT(row(ali,0) == "tc");
            SEQAN_ASSERT(row(ali,1) == "tc");
        }
        if(i == 3){
            SEQAN_ASSERT_EQ(score, 4);
            SEQAN_ASSERT(row(ali,0) == "gc");
            SEQAN_ASSERT(row(ali,1) == "gc");
        }
        if(i == 4){
            SEQAN_ASSERT_EQ(score, 4);
            SEQAN_ASSERT(row(ali,0) == "ag");
            SEQAN_ASSERT(row(ali,1) == "ag");
        }
        if(i == 5){
            SEQAN_ASSERT_EQ(score, 4);
            SEQAN_ASSERT(row(ali,0) == "taagcgtctcg");
            SEQAN_ASSERT(row(ali,1) == "tcatagagttg");
        }
        if(i == 6){
            SEQAN_ASSERT_EQ(score, 3);
            SEQAN_ASSERT(row(ali,0) == "ata");
            SEQAN_ASSERT(row(ali,1) == "aga");
        }
        if(i == 7){
            SEQAN_ASSERT_EQ(score, 3);
            SEQAN_ASSERT(row(ali,0) == "cgt");
            SEQAN_ASSERT(row(ali,1) == "cat");
        }
        if(i == 8){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "g");
            SEQAN_ASSERT(row(ali,1) == "g");
        }
        if(i == 9){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "t");
            SEQAN_ASSERT(row(ali,1) == "t");
        }
        if(i == 10){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "a");
            SEQAN_ASSERT(row(ali,1) == "a");
        }
        if(i == 11){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "c");
            SEQAN_ASSERT(row(ali,1) == "c");
        }
        if(i == 12){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "a");
            SEQAN_ASSERT(row(ali,1) == "a");
        }
        if(i == 13){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "t");
            SEQAN_ASSERT(row(ali,1) == "t");
        }
        if(i == 14){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "c");
            SEQAN_ASSERT(row(ali,1) == "c");
        }
        if(i == 15){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "g");
            SEQAN_ASSERT(row(ali,1) == "g");
        }
        if(i == 16){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "t");
            SEQAN_ASSERT(row(ali,1) == "t");
        }
        if(i == 17){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "a");
            SEQAN_ASSERT(row(ali,1) == "a");
        }
        if(i == 18){
            SEQAN_ASSERT_EQ(score, 2);
            SEQAN_ASSERT(row(ali,0) == "t");
            SEQAN_ASSERT(row(ali,1) == "t");
        }
        ++i;
    }

//test if every cell has been reduced to 0
//only makes sense if cutoff=0
    if(cutoff==0){
        int str0len = length(str0) + 1;
        int str1len = length(str1) + 1;
        bool check = true;
        for(int i = 0; i <str1len; ++i){
            for(int j=0;j<str0len;++j){
                if(getValue(sw_finder.matrix,(i*str0len)+j)!=0){
                    check = false;
                }
            }
        }

        SEQAN_ASSERT(check == true);

    }

//desweiteren nur so:
    push(sw_finder.pQ,LocalAlignmentFinder<int>::TPQEntry());
    SEQAN_ASSERT(empty(sw_finder.pQ) == false);
    clear(sw_finder.pQ);
    SEQAN_ASSERT(empty(sw_finder.pQ) == true);




}


SEQAN_DEFINE_TEST(testLocalAlign2) {
//new interface

    String<char> str0 = "ataagcgtctcg";
    String<char> str1 = "tcatagagttgc";

    Align< String<char>, ArrayGaps> ali;
    resize(rows(ali), 2);
    setSource(row(ali, 0), str0);
    setSource(row(ali, 1), str1);

    Score<int> score_type = Score<int>(2,-1,-2,0) ;
    LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>();

    int score = localAlignment(ali, sw_finder, score_type, 5);
    SEQAN_ASSERT_EQ(score, 9);
    SEQAN_ASSERT(row(ali,0) == "ataagcgt");
    SEQAN_ASSERT(row(ali,1) == "ata-gagt");

    score = localAlignment(ali, sw_finder, score_type, 5);
    SEQAN_ASSERT_EQ(score, 5);
    SEQAN_ASSERT(row(ali,0) == "tc-tcg");
    SEQAN_ASSERT(row(ali,1) == "tcatag");

    score = localAlignment(ali, sw_finder, score_type, 5, WatermanEggert());
    SEQAN_ASSERT_EQ(score, 0);
}


SEQAN_DEFINE_TEST(testBandedLocalAlign) {
    typedef String<Dna> TString;
    TString str0("ggggcttaagcttgggg");
    TString str1("aaaacttagctctaaaa");

    Score<int> score_type = Score<int>(2,-1,-2,-2);

    Align<TString> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), str0);
    assignSource(row(align, 1), str1);

    LocalAlignmentFinder<int> finder = LocalAlignmentFinder<int>();

    int score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 12);
    SEQAN_ASSERT(row(align, 0) == "cttaagct");
    SEQAN_ASSERT(row(align, 1) == "ctt-agct");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 10);
    SEQAN_ASSERT(row(align, 0) == "aagcttgg");
    SEQAN_ASSERT(row(align, 1) == "aaacttag");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 10);
    SEQAN_ASSERT(row(align, 0) == "gct-taa");
    SEQAN_ASSERT(row(align, 1) == "gctctaa");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 5);
    SEQAN_ASSERT(row(align, 0) == "aagcttg");
    SEQAN_ASSERT(row(align, 1) == "aacttag");

    score = localAlignment(align, finder, score_type, 5, -6, 6, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 0);



    Align<TString> align1;
    resize(rows(align1), 2);
    assignSource(row(align1, 0), "gcagaattaaggaggattacaagtgggaatttgaagagcttttgaaatcc");
    assignSource(row(align1, 1), "cggttgagcagaacttgggctacgagactccccccgaggaatttgaaggctttcttcaaatccaaaagca");

    Score<int> score_type1 = Score<int>(1,-9,-9,-9);

    LocalAlignmentFinder<int> finder1 = LocalAlignmentFinder<int>();
    score = localAlignment(align1, finder1, score_type1, 7, -20, 0, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 11);
    SEQAN_ASSERT(row(align1, 0) == "ggaatttgaag");
    SEQAN_ASSERT(row(align1, 1) == "ggaatttgaag");
    score = localAlignment(align1, finder1, score_type1, 7, -20, 0, BandedWatermanEggert());
    SEQAN_ASSERT_EQ(score, 0);
}


#endif
