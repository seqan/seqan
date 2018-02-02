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

#ifndef SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H
#define SEQAN_HEADER_TEST_GRAPH_TCOFFEE_H

#include <seqan/random.h>

using namespace seqan;



//////////////////////////////////////////////////////////////////////////////

void Test_Distances() {
    typedef String<AminoAcid> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;

    TString str1 = "GARFIELDTHELASTFATCAT";
    TString str2 = "GARFIELDTHEFASTCAT";
    TString str3 = "GARFIELDTHEVERYFASTCAT";
    TString str4 = "THEFATCAT";
    TStringSet strSet;
    assignValueById(strSet, str1);
    assignValueById(strSet, str2);
    assignValueById(strSet, str3);
    assignValueById(strSet, str4);
    TGraph g(strSet);

    String<double> distanceMatrix;
    getDistanceMatrix(g,distanceMatrix);
    SEQAN_ASSERT(distanceMatrix[3] ==  1.0 - 5.0 / 7.0);
    SEQAN_ASSERT(distanceMatrix[1 * length(strSet) + 3] == 1.0 - 5.0 / 7.0);
    SEQAN_ASSERT(distanceMatrix[2 * length(strSet) + 3] == 1.0 - 3.0 / 7.0);
    clear(distanceMatrix);
    String<unsigned int> pList;
    selectPairs(strSet, pList);
    Blosum62 score_type(-1,-11);
    String<Fragment<> > matches;
    String<int> scores;
    String<double> dist;
    appendSegmentMatches(strSet, pList, score_type, matches, scores, dist, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    getDistanceMatrix(g,distanceMatrix,LibraryDistance());
    SEQAN_ASSERT(getValue(distanceMatrix, 0 * length(strSet) + 1) < getValue(distanceMatrix, 2 * length(strSet) + 3));
    SEQAN_ASSERT(getValue(distanceMatrix, 1 * length(strSet) + 2) < getValue(distanceMatrix, 2 * length(strSet) + 3));
}


//////////////////////////////////////////////////////////////////////////////

void
testquickAlign__(Graph<Alignment<StringSet<String<AminoAcid>, Dependent<> >, unsigned int> >& g)
{
    Graph<Alignment<StringSet<String<AminoAcid>, Dependent<> >, void, WithoutEdgeId> > gOut(stringSet(g));
    tripletLibraryExtension(g);
    String<double> distForGuideTree;
    getDistanceMatrix(g,distForGuideTree,LibraryDistance());
    Graph<Tree<double> > guideTree;
    njTree(distForGuideTree, guideTree);
    progressiveAlignment(g, guideTree, gOut);
    //std::cout << gOut << std::endl;
    String<char> alignMat;
    convertAlignment(gOut,alignMat);
    unsigned int len = length(alignMat) / 4;
    SEQAN_ASSERT(String<char>(infix(alignMat, 0, 8)) == "GARFIELD");
    SEQAN_ASSERT(String<char>(infix(alignMat, 1*len + 0, 1*len+8)) == "GARFIELD");
    SEQAN_ASSERT(String<char>(infix(alignMat, 2*len + 0, 2*len+8)) == "GARFIELD");
    SEQAN_ASSERT(String<char>(infix(alignMat, 3*len + 0, 3*len+8)) == "--------");

    //std::cout << gOut << std::endl;
}

void Test_Libraries() {
    typedef String<AminoAcid> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, unsigned int> > TGraph;

    TString str1 = "GARFIELDTHELASTFATCAT";
    TString str2 = "GARFIELDTHEFASTCAT";
    TString str3 = "GARFIELDTHEVERYFASTCAT";
    TString str4 = "THEFATCAT";
    TStringSet strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);
    appendValue(strSet, str3);
    appendValue(strSet, str4);
    String<unsigned int> pList;
    selectPairs(strSet, pList);
    TGraph g(strSet);
    Blosum62 score_type(-1,-11);
    String<Fragment<> > matches;
    String<int> scores;
    appendSegmentMatches(strSet, pList, matches, scores, LcsLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
    clear(matches);
    clear(scores);
    appendSegmentMatches(strSet, matches, scores, KmerLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
    clear(matches);
    clear(scores);
    appendSegmentMatches(strSet, pList, score_type, matches, scores, LocalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
    clear(matches);
    clear(scores);
    Nothing noth;
    appendSegmentMatches(strSet, pList, score_type, matches, scores, noth, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
    clear(matches);
    clear(scores);
    String<double> distanceMatrix;
    appendSegmentMatches(strSet, pList, score_type, matches, scores, distanceMatrix, AlignConfig<false,false,false,false>(), GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
    Graph<Undirected<double> > distGraph;
    clear(matches);
    clear(scores);
    appendSegmentMatches(strSet, pList, score_type, matches, scores, distanceMatrix, AlignConfig<false,false,false,false>(), GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    testquickAlign__(g);
}

void Test_ExternalLibraries() {
    typedef String<char> TName;
    typedef StringSet<TName, Owner<> > TNameSet;
    typedef String<AminoAcid> TSequence;
    typedef StringSet<TSequence, Owner<> > TSequenceSet;

    TSequenceSet seqSet;
    appendValue(seqSet, "GARFIELDTHELASTFATCAT");
    appendValue(seqSet, "GARFIELDTHEFASTCAT");
    appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
    appendValue(seqSet, "THEFATCAT");
    String<unsigned int> pList;
    selectPairs(seqSet, pList);
    TNameSet nameSet;
    appendValue(nameSet, "seq1");
    appendValue(nameSet, "seq2");
    appendValue(nameSet, "seq3");
    appendValue(nameSet, "seq4");

    typedef Graph<Alignment<StringSet<TSequence, Dependent<> >, unsigned int> > TGraph;
    TGraph g(seqSet);
    Blosum62 score_type(-1,-11);
    String<Fragment<> > matches;
    String<int> scores;
    Nothing noth;
    appendSegmentMatches(stringSet(g), pList, score_type, matches, scores, noth, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    //std::cout << g << std::endl;

    // T-Coffee lib format
    std::string fileName = SEQAN_TEMP_FILENAME();

    // Writing
    std::fstream strm;
    strm.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
    write(strm, g, nameSet, TCoffeeLib());
    strm.close();
    //_debugRefinedMatches(g);
    testquickAlign__(g);

    // Reading
    clear(seqSet);
    clear(nameSet);
    std::fstream strmRead;
    strmRead.open(fileName.c_str(), std::ios_base::in | std::ios_base::binary);
    read(strmRead,seqSet,nameSet,TCoffeeLib());
    strmRead.close();
    //std::cout << g << std::endl;


    // Blast format
    fileName = SEQAN_TEMP_FILENAME();

    // Writing
    std::fstream strm2;
    strm2.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
    write(strm2, g, nameSet, BlastLib());
    strm2.close();

    // Reading
    std::fstream strmRead2;
    strmRead2.open(fileName.c_str(), std::ios_base::in | std::ios_base::binary);
    clearVertices(g);
    clear(matches);
    clear(scores);
    read(strmRead2, matches, scores, nameSet,BlastLib());
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    strmRead2.close();
    // TODO(holtgrew): Disabled output here, need to compare g before writing out and a g2 after reading in.
    // std::cout << g << std::endl;

    clear(g);
    assignStringSet(g, seqSet);
    clear(matches);
    clear(scores);
    appendSegmentMatches(stringSet(g), pList, score_type, matches, scores, noth, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    // Blast lib format

    // Writing
    fileName = SEQAN_TEMP_FILENAME();
    std::fstream strmBlast;
    strmBlast.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc);
    write(strmBlast, g, nameSet, BlastLib());
    strmBlast.close();
    //_debugRefinedMatches(g);
    testquickAlign__(g);

    // Reading
    clear(g);
    assignStringSet(g, seqSet);
    std::fstream strmBlastLib;
    strmBlastLib.open(fileName.c_str(), std::ios_base::in | std::ios_base::binary);
    clear(matches);
    clear(scores);
    read(strmBlastLib, matches, scores, nameSet,BlastLib());
    buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
    strmBlastLib.close();
    //_debugRefinedMatches(g);
    testquickAlign__(g);

}

void Test_TripletExtension() {
    typedef String<char> TName;
    typedef StringSet<TName, Owner<> > TNameSet;
    typedef String<AminoAcid> TSequence;
    typedef StringSet<TSequence, Owner<> > TSequenceSet;
    typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;

    TSequenceSet seqSet;//1234567890
    appendValue(seqSet, "GARFIELDTHE");
    appendValue(seqSet, "GARFIELDTHE");
    appendValue(seqSet, "GARFIELDTHE");
    appendValue(seqSet, "THE");
    TNameSet nameSet;
    appendValue(nameSet, "seq1");
    appendValue(nameSet, "seq2");
    appendValue(nameSet, "seq3");
    appendValue(nameSet, "seq4");

    // Triplet extension
    typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
    TGraph g(seqSet);
    addEdge(g, addVertex(g, 0, 0, 8), addVertex(g, 1, 0, 8), 40);
    addEdge(g, findVertex(g, 0, 0), addVertex(g, 2, 0, 8), 30);
    addEdge(g, addVertex(g, 1, 8, 3), addVertex(g, 0, 8, 3), 40);
    addEdge(g, findVertex(g, 1, 8), addVertex(g, 2, 8, 3), 30);
    addEdge(g, findVertex(g, 1, 8), addVertex(g, 3, 0, 3), 20);
    addEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0), 40);
    tripletLibraryExtension(g);
    SEQAN_ASSERT(cargo(findEdge(g, findVertex(g, 1, 0), findVertex(g, 2, 0))) == 30);
    SEQAN_ASSERT(cargo(findEdge(g, findVertex(g, 0, 8), findVertex(g, 3, 0))) == 20);
    SEQAN_ASSERT(cargo(findEdge(g, findVertex(g, 0, 8), findVertex(g, 2, 8))) == 30);
    SEQAN_ASSERT(cargo(findEdge(g, findVertex(g, 2, 8), findVertex(g, 3, 0))) == 60);
}

void Test_SumOfPairsScore()
{
    typedef String<char> TName;
    typedef StringSet<TName, Owner<> > TNameSet;
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence, Owner<> > TSequenceSet;
    typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;
    typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;

    TSequenceSet seqSet;
    appendValue(seqSet, "ACAAGTA");
    appendValue(seqSet, "AA");
    appendValue(seqSet, "ACCTA");
    TNameSet nameSet;
    appendValue(nameSet, "seq1");
    appendValue(nameSet, "seq2");
    appendValue(nameSet, "seq3");

    TGraph g(seqSet);
    Score<int> score_type = Score<int>(5,-4,-2,-10);
    String<double> distanceMatrix;
    String<Fragment<> > matches;
    String<int> scores;
    String<unsigned int> pList;
    selectPairs(seqSet, pList);
    appendSegmentMatches(stringSet(g), pList, score_type, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    Graph<Tree<double> > guideTree;
    njTree(distanceMatrix, guideTree);
    TGraph gOut(seqSet);
    progressiveAlignment(g, guideTree, gOut);
    SEQAN_ASSERT(sumOfPairsScore(gOut, score_type) == -8);
    SEQAN_ASSERT(sumOfPairsScoreInd(gOut, score_type) == 16);

    seqSet[1] = "AAG";
    Score<int> scType = Score<int>(5,-4,-1,-2);
    clear(distanceMatrix);
    clearVertices(g);
    clear(matches);
    clear(scores);
    appendSegmentMatches(stringSet(g), pList, scType, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    clear(guideTree);
    njTree(distanceMatrix, guideTree);
    clearVertices(gOut);
    progressiveAlignment(g, guideTree, gOut);
    SEQAN_ASSERT(sumOfPairsScore(gOut, scType) == 20);

    resize(seqSet, 2);
    seqSet[0] = "TTT";
    seqSet[1] = "AAA";
    clear(distanceMatrix);
    clear(g);
    assignStringSet(g, seqSet);
    clear(matches);
    clear(scores);
    selectPairs(seqSet, pList);
    appendSegmentMatches(stringSet(g), pList, scType, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    clear(guideTree);
    njTree(distanceMatrix, guideTree);
    clear(gOut);
    assignStringSet(gOut, seqSet);
    progressiveAlignment(g, guideTree, gOut);
    SEQAN_ASSERT(sumOfPairsScore(gOut, scType) == -8);

    seqSet[0] = "TTTAAATTT";
    seqSet[1] = "AAA";
    clear(distanceMatrix);
    clearVertices(g);
    clear(matches);
    clear(scores);
    selectPairs(seqSet, pList);
    appendSegmentMatches(stringSet(g), pList, scType, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    clear(guideTree);
    njTree(distanceMatrix, guideTree);
    clearVertices(gOut);
    progressiveAlignment(g, guideTree, gOut);
    SEQAN_ASSERT(sumOfPairsScore(gOut, scType) == 7);

    seqSet[0] = "AAAAAA";
    seqSet[1] = "TTTAAATTTAAATTT";
    clear(distanceMatrix);
    clearVertices(g);
    clear(matches);
    clear(scores);
    selectPairs(seqSet, pList);
    appendSegmentMatches(stringSet(g), pList, scType, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    clear(guideTree);
    njTree(distanceMatrix, guideTree);
    clearVertices(gOut);
    progressiveAlignment(g, guideTree, gOut);
    SEQAN_ASSERT(sumOfPairsScore(gOut, scType) == 18);
}


void Test_Progressive() {
    typedef String<char> TName;
    typedef StringSet<TName, Owner<> > TNameSet;
    typedef String<AminoAcid> TSequence;
    typedef StringSet<TSequence, Owner<> > TSequenceSet;
    typedef StringSet<TSequence, Dependent<> > TDependentSequenceSet;

    TSequenceSet seqSet;
    appendValue(seqSet, "GARFIELDTHELASTFATCAT");
    appendValue(seqSet, "GARFIELDTHEFASTCAT");
    appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
    appendValue(seqSet, "THEFATCAT");
    appendValue(seqSet, "GARFIELDTHELASTFATCAT");
    appendValue(seqSet, "GARFIELDTHEFASTCAT");
    appendValue(seqSet, "GARFIELDTHEVERYFASTCAT");
    appendValue(seqSet, "THEFATCAT");
    TNameSet nameSet;
    appendValue(nameSet, "seq1");
    appendValue(nameSet, "seq2");
    appendValue(nameSet, "seq3");
    appendValue(nameSet, "seq4");
    appendValue(nameSet, "seq5");
    appendValue(nameSet, "seq6");
    appendValue(nameSet, "seq7");
    appendValue(nameSet, "seq8");

    typedef Graph<Alignment<TDependentSequenceSet, unsigned int> > TGraph;
    TGraph g(seqSet);
    Blosum62 score_type(-1,-11);
    String<double> distanceMatrix;
    String<Fragment<> > matches;
    String<int> scores;
    String<unsigned int> pList;
    selectPairs(seqSet, pList);
    appendSegmentMatches(stringSet(g), pList, score_type, matches, scores, distanceMatrix, GlobalPairwiseLibrary() );
    buildAlignmentGraph(matches, scores, g, FractionalScore() );
    tripletLibraryExtension(g);
    Graph<Tree<double> > guideTree;
    njTree(distanceMatrix, guideTree);
    TGraph gOut(seqSet);
    progressiveAlignment(g, guideTree, gOut);
    sumOfPairsScore(gOut, score_type);
}


void Test_ReversableFragments() {
    typedef unsigned int TSize;
    typedef String<Dna> TSequence;
    TSequence seq1 = "AACGTT";
    TSequence seq2 = "AACGTTC";
    typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
    TDepSequenceSet strSet;
    appendValue(strSet, seq1);
    appendValue(strSet, seq2);
    typedef Fragment<TSize, ExactReversableFragment<> > TFragment;
    typedef String<TFragment> TFragmentString;
    TFragmentString matches;
    appendValue(matches, TFragment(0,0,1,0,2));
    appendValue(matches, TFragment(0,3,1,3,3));
    appendValue(matches, TFragment(0,1,1,2,3, true));
    typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
    TGraph g(strSet);
    matchRefinement(matches,strSet,g);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 2)) == 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 4)) != 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 2), findVertex(g, 1, 3)) != 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 3), findVertex(g, 1, 2)) != 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 3), findVertex(g, 1, 4)) == 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 0), findVertex(g, 1, 0)) != 0);
    SEQAN_ASSERT(findEdge(g, findVertex(g, 0, 1), findVertex(g, 1, 1)) != 0);
}

SEQAN_DEFINE_TEST(test_distances)
{
    Test_Distances();
}

SEQAN_DEFINE_TEST(test_libraries)
{
    Test_Libraries();
}
SEQAN_DEFINE_TEST(test_external_libraries)
{
    Test_ExternalLibraries();
}
SEQAN_DEFINE_TEST(test_triplet_extension)
{
    Test_TripletExtension();
}
SEQAN_DEFINE_TEST(test_sop)
{
    Test_SumOfPairsScore();
}
SEQAN_DEFINE_TEST(test_progressive)
{
    Test_Progressive();
}
SEQAN_DEFINE_TEST(test_reversable_fragments)
{
    Test_ReversableFragments();
}

#endif

