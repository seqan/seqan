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

#ifndef SEQAN_HEADER_TEST_GRAPH_MATCH_GRAPH_ALIGN_H
#define SEQAN_HEADER_TEST_GRAPH_MATCH_GRAPH_ALIGN_H

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>

using namespace std;
using namespace seqan;

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(RefineMatchesSelfEdges)
{
    typedef String<char> TString;
    typedef StringSet<TString > TStringSet;
    typedef StringSet<TString, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAliGraph;
    typedef Fragment<> TFragment;
    typedef String<TFragment> TFragString;


    TStringSet seq_set;
    appendValue(seq_set,String<char>("aaaaabbbbbbccccccaaaaaa"));
    appendValue(seq_set,String<char>("aaaabbbbbbaaaaaa"));


    TFragString matches;
    appendValue(matches,TFragment(0,0,1,0,4));
    appendValue(matches,TFragment(0,0,1,11,5));
    appendValue(matches,TFragment(0,17,1,10,6));
    appendValue(matches,TFragment(0,5,1,4,6));

    // without within-sequence-matches: 12 vertices, 7 edges
    TAliGraph ali_graph(seq_set);
    matchRefinement(matches,seq_set,ali_graph);

      SEQAN_ASSERT_EQ(numVertices(ali_graph), (unsigned)12);
    SEQAN_ASSERT_EQ(numEdges(ali_graph), (unsigned)7);

    // with within-sequence-match: 24 vertices, 20 edges
    appendValue(matches,TFragment(1,0,1,10,4));

    TAliGraph ali_graph2(seq_set);
    matchRefinement(matches,seq_set,ali_graph2);

    SEQAN_ASSERT_EQ(numVertices(ali_graph2), (unsigned) 24);
    SEQAN_ASSERT_EQ(numEdges(ali_graph2),  (unsigned)20);

//    std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
//    std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
//    std::cout << ali_graph << "\n\n";

//    std::cout << "\nnumEdges: "<<numEdges(ali_graph2)<<"\n";
//    std::cout << "\nnumVertices: "<<numVertices(ali_graph2)<<"\n";
//    std::cout << ali_graph2 << "\n\n";

}

//////////////////////////////////////////////////////////////////////////////////////////
// test the graph_align heuristic
SEQAN_DEFINE_TEST(RefineInexactFragment)
{
    typedef String<char> TString;
    typedef StringSet<TString > TStringSet;
    typedef StringSet<TString, Dependent<> > TDepStringSet;
    typedef Graph<Alignment<TDepStringSet> > TAliGraph;
    typedef Fragment<> TFragment;
    typedef String<TFragment> TFragString;


    TStringSet seq_set;
    appendValue(seq_set,String<char>("aaaaabbbbbbccccccaaaaaa"));
    appendValue(seq_set,String<char>("aaaabbbbbbaaaaaa"));


    TFragString matches;
    appendValue(matches,TFragment(0,0,1,0,4));
    appendValue(matches,TFragment(0,0,1,11,5));
    appendValue(matches,TFragment(0,17,1,10,6));
    appendValue(matches,TFragment(0,5,1,4,6));

    TAliGraph ali_graph(seq_set);
    matchRefinement(matches,seq_set,ali_graph,2); // min_frag_len = 2

//    std::cout << ali_graph << std::endl;

    SEQAN_ASSERT_EQ(numVertices(ali_graph), (unsigned) 7);
    SEQAN_ASSERT_EQ(numEdges(ali_graph), (unsigned) 4);

    TAliGraph ali_graph2(seq_set);
    matchRefinement(matches,seq_set,ali_graph2,5); // min_frag_len = 5

//    std::cout << ali_graph2 << std::endl;

    SEQAN_ASSERT_EQ(numVertices(ali_graph2), (unsigned)6);
    SEQAN_ASSERT_EQ(numEdges(ali_graph2),  (unsigned)3); // the first fragment (4bp long) disappears


}



//////////////////////////////////////////////////////////////////////////////////77

//// TODO(holtgrew): Is this a helper? What does it do?
//int Test_ConvertSequences(String<char> const in_path, String<char> const in_file, String<char> const path, String<char> const file_prefix) {
//    typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
//
//    // count sequences
//    unsigned seqCount = 0;
//
//    SeqFileIn file;
//    std::stringstream input;
//    input << in_path << in_file;
//    if (!open(file, input.str().c_str())) return 0;
//    String<char> id, seq;
//    while (!atEnd(file))
//    {
//        readRecord(id, seq, file);
//        std::cout << id << std::endl;
//        ++seqCount;
//    }
//    std::cout << "Number of sequences: " << seqCount << std::endl;
//
//    // import sequences
//    file.clear();
//    file.seekg(0, ios_base::beg);
//    for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i)
//    {
//        TString str;
//        //String<TraceBack, External<> > trace;
//        //open(trace, "D:\\seqan.dat");
//        std::stringstream s;
//        s << path << file_prefix << i;
//        open(str, s.str().c_str());
//        read(file, str, Fasta());
//    }
//    file.close();
//
//    return seqCount;
//}

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Is this a helper? What does it do?
template<typename TStringSet, typename TVal1, typename TVal2>
inline bool
Test_ReadSequences(String<char> const path, String<char> const file_prefix, TStringSet& str, TVal1 const start, TVal2 const nseq) {
    for(unsigned i = start; i < start + nseq; ++i) {
        std::stringstream s;
        s << path << file_prefix << i - start;
        bool f = open(str[i], s.str().c_str());
        if (!f) return false;
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(GraphMatchRefine)
{
    SEQAN_SKIP_TEST;
    // Sequences
    typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
    typedef StringSet<TString, Owner<> > TStringSet;
    typedef Id<TStringSet>::Type TId;
    typedef Size<TStringSet>::Type TSize;

    // Matches
    typedef String<Fragment<>, External<ExternalConfig<File<>, 64*1024> > > TFragmentString;
    //typedef String<Fragment<>, External<> > TFragmentString;

    // Windows
#ifdef STDLIB_VS
    String<char> in_path("Z:\\matches\\");
    String<char> out_path("Z:\\matches\\out\\");
    return;
#else
    // Linux
    String<char> in_path("/home/takifugu2/data/SeqAn/binary/");
    String<char> out_path("/home/takifugu2/data/SeqAn/binary/");
#endif


    TSize hSeq = 24;
    TSize wSeq = 24;
    TSize bSeq = 24;


    // Convert all sequences only once
    //TSize tmp = 0;
    //tmp = Test_ConvertSequences(in_path, "HUREF6CHROM.fasta",out_path,"H.chr.");
    //if (tmp != hSeq) {
    //    std::cerr << "Did not read all HUREF sequences." << std::endl;
    //    exit(1);
    //}
    //tmp = Test_ConvertSequences(in_path, "WGSACHROM.fasta",out_path,"W.chr.");
    //if (tmp != wSeq) {
    //    std::cerr << "Did not read all WGSA sequences." << std::endl;
    //    exit(1);
    //}
    //tmp = Test_ConvertSequences(in_path, "B36LCCHROM.fasta",out_path,"B.chr.");
    //if (tmp != bSeq) {
    //    std::cerr << "Did not read all B36LC sequences." << std::endl;
    //    exit(1);
    //}

    // Read all the sequences
    TStringSet str;
    resize(str, hSeq + wSeq + bSeq);
    bool f = Test_ReadSequences(out_path,"H.chr.", str, 0, hSeq);
    if (!f) {
        std::cerr << "Error importing HUREF sequences." << std::endl;
        exit(1);
    }
    f = Test_ReadSequences(out_path,"W.chr.", str, hSeq, wSeq);
    if (!f) {
        std::cerr << "Error importing WGSA sequences." << std::endl;
        exit(1);
    }
    f = Test_ReadSequences(out_path,"B.chr.", str, hSeq + wSeq, bSeq);
    if (!f) {
        std::cerr << "Error importing B36LC sequences." << std::endl;
        exit(1);
    }

    // Build a map:
    // SeqId -> Identifier
    typedef std::map<TId, String<char> > TIdToNameMap;
    TIdToNameMap idToName;
    for(TId i = 0; i < length(str); ++i) {
        TSize index = 0;
        std::stringstream s;
        if (i < 24) {
            s << "H";
            index = i;
        }
        else if (i < 48) {
            s << "W";
            index = i - hSeq;
        }
        else {
            s << "B";
            index = i - (hSeq + wSeq);
        }
        s << ":" << index;
        idToName.insert(std::make_pair(i, s.str().c_str()));
    }

    // Just to check that everything worked
    std::cout << "Number of sequences: " << length(str) << std::endl;
    for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
        std::cout << pos->second << ") ";
        for(TSize i=0; ((i<10) && (i <length(str[pos->first])));++i) {
            std::cout << str[pos->first][i];
        }
        std::cout << std::endl;
    }

    // Access the matches
    TFragmentString matches;
    std::stringstream strstream;
    strstream << out_path << "matchesTest.dat"; // 10 Matches
    //strstream << out_path << "matches1000.dat"; // 2001948 Matches
    //strstream << out_path << "matches10000.dat"; // 2111 Matches
    //strstream << out_path << "matches2000.dat"; // 653095 Matches
    //strstream << out_path << "matches500.dat"; // 3999176
    open(matches, strstream.str().c_str());


    // Convert the matches to an external string
    //for(TSize i = 1; i<4; ++i) {
    //    fstream strm;
    //    std::stringstream s;
    //    if (i==0 ) s << in_path << "TvsT.atac";
    //    else if (i==1 ) s << in_path << "BvsH.atac";
    //    else if (i==2 ) s << in_path << "BvsW.atac";
    //    else if (i==3 ) s << in_path << "WvsH.atac";
    //    strm.open(s.str().c_str(), ios_base::in);
    //    read(strm, matches, 500, AtacMatches());
    //    strm.close();
    //}

    // Print number of matches
    std::cout << "Number of matches: " << length(matches) << std::endl;
//    for(unsigned i = 0; i < length(matches); ++i)
//        printMatch(matches[i]);
    // Re7finement
    //typedef Infix<TString>::Type TInfix;
    typedef StringSet<TString, Dependent<> > TAlignmentStringSet;
    typedef Graph<Alignment<TAlignmentStringSet> > TAliGraph;
    //typedef VertexDescriptor<TAliGraph>::Type TVD;
    TAlignmentStringSet aliStr;
    for(TId i = 0; i<length(str); ++i) {
        assignValueById(aliStr, str, i);
    }
    Score<int> score_type = Score<int>(1,-1,-2,0) ;
    TAliGraph ali_graph(aliStr);
    matchRefinement(matches,str,score_type,ali_graph);//,StoreEdges());
    std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
    std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
    //std::cout << ali_graph <<"\n";

    // Print all the matches
    //typedef Iterator<TAliGraph, EdgeIterator>::Type TEdgeIterator;
    //TEdgeIterator it(ali_graph);
    //for(;!atEnd(it);goNext(it)) {
    //    TVD sV = sourceVertex(it);
    //    TVD tV = targetVertex(it);
    //    TId seqId1 = sequenceId(ali_graph,sV);
    //    TId seqId2 = sequenceId(ali_graph,tV);
    //    TSize seqBegin1 = fragmentBegin(ali_graph, sV);
    //    TSize seqBegin2 = fragmentBegin(ali_graph, tV);
    //    TSize len = fragmentLength(ali_graph, sV);
    //    TIdToNameMap::const_iterator pos1 =  idToName.find(seqId1);
    //    TIdToNameMap::const_iterator pos2 =  idToName.find(seqId2);
    //    std::cout << pos1->second << "," << seqBegin1  << "," << len << "," << pos2->second << "," << seqBegin2 << "," << len << std::endl;
    //}

    for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
        close(str[pos->first]);
    }
    close(matches);
}


//produce pairwise alignments (Align object)
template<typename TAlign, typename TSequence, typename TSeqSpec, typename TScore>
void
getAlignments(String<TAlign> & alis, StringSet<TSequence, TSeqSpec> & seq, TScore & score_type, int & numAlignments, int cutoff)
{

    unsigned int gesamt = 0;

    for(unsigned int i = 0; i < length(seq); ++i)
    {
        for(unsigned int j = i+1; j < length(seq); ++j)
        {
            TAlign ali;
            resize(rows(ali), 2);
            setSource(row(ali, 0), seq[i]);
            setSource(row(ali, 1), seq[j]);
            LocalAlignmentEnumerator<TScore, Unbanded> enumerator(score_type, cutoff);

            if (!nextLocalAlignment(ali, enumerator))
                continue;
            //cout << ali<<"\n";
            //cout <<"Seq "<<i<<" - Seq "<<j<<"\n"<<score<< ali;
            //cout << clippedBeginPosition(row(ali,0)) <<"   ";
            //cout << clippedBeginPosition(row(ali,1)) <<"\n";
            appendValue(alis,ali);
            ++gesamt;
            int k = 1;
            while(k<numAlignments)
            {
                if (!nextLocalAlignment(ali, enumerator))
                    break;
                //cout <<score<< ali;
                //cout << clippedBeginPosition(row(ali,0)) <<"   ";
                //cout << clippedBeginPosition(row(ali,1)) <<"\n";
                appendValue(alis,ali);
                ++gesamt;
                ++k;
            }
        }
    }

    numAlignments = gesamt;
    resize(alis,numAlignments);


}


SEQAN_DEFINE_TEST(RefineAlign)
{
    SEQAN_SKIP_TEST;
    typedef String<char> TString;
    typedef StringSet<TString> TStringSet;
    typedef StringSet<TString,Dependent<> > TDepStringSet;
    typedef Align<TString, ArrayGaps> TAlign;

    int numSequences = 4;

    TStringSet seq_set;


    appendValue(seq_set,String<char>("GARFIELDTHELASTFATCAT"));
    appendValue(seq_set,String<char>("GARFIELDTHEFASTCAT"));
    appendValue(seq_set,String<char>("GARFIELDTHEVERYFASTCAT"));
    appendValue(seq_set,String<char>("THEFATCAT"));



    int numAlignments = 1;
    int numSequencePairs = 0;
    int cutoff = 4;
    for(int i = 1 ; i < numSequences; ++i)
        numSequencePairs += i;
    String<TAlign> alis;
    reserve(alis,numSequencePairs*numAlignments);
    Score<int> score_type = Score<int>(1,-1,-2,-2) ;

    getAlignments(alis,seq_set,score_type,numAlignments,cutoff);

    typedef Graph<Alignment<TDepStringSet> > TAliGraph;
    TAliGraph ali_graph(seq_set);

    //std::cout <<"Number of Segments: "<<length(alis)<<"\n";
    matchRefinement(alis,seq_set,score_type,ali_graph);

    //std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
    //std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
    //std::cout << ali_graph << "\n\n";

    //fstream strmW; // Write the library
    //strmW.open(TEST_PATH "my_testlib1.lib", ios_base::out | ios_base::trunc);
    //write(strmW,ali_graph,TCoffeeLib());
    //strmW.close();

    VertexDescriptor<TAliGraph>::Type vd;

    vd = findVertex(ali_graph,0,0);
    SEQAN_ASSERT(vd == 0);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8);
    vd = findVertex(ali_graph,0,8);
    SEQAN_ASSERT(vd == 1);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,0,11);
    SEQAN_ASSERT(vd == 2);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,0,14);
    SEQAN_ASSERT(vd == 3);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);
    vd = findVertex(ali_graph,0,15);
    SEQAN_ASSERT(vd == 4);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,0,17);
    SEQAN_ASSERT(vd == 5);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);
    vd = findVertex(ali_graph,0,18);
    SEQAN_ASSERT(vd == 6);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,0,20);
    SEQAN_ASSERT(vd == 7);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);

    vd = findVertex(ali_graph,1,0);
    SEQAN_ASSERT(vd == 8);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8);
    vd = findVertex(ali_graph,1,8);
    SEQAN_ASSERT(vd == 9);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,1,11);
    SEQAN_ASSERT(vd == 10);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 11);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,1,14);
    SEQAN_ASSERT(vd == 11);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 14);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);
    vd = findVertex(ali_graph,1,15);
    SEQAN_ASSERT(vd == 12);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 15);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,1,17);
    SEQAN_ASSERT(vd == 13);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 17);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);

    vd = findVertex(ali_graph,2,0);
    SEQAN_ASSERT(vd == 14);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8);
    vd = findVertex(ali_graph,2,8);
    SEQAN_ASSERT(vd == 15);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,2,11);
    SEQAN_ASSERT(vd == 16);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 11);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 7);
    vd = findVertex(ali_graph,2,18);
    SEQAN_ASSERT(vd == 17);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 18);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);
    vd = findVertex(ali_graph,2,19);
    SEQAN_ASSERT(vd == 18);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 19);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,2,21);
    SEQAN_ASSERT(vd == 19);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 21);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);

    vd = findVertex(ali_graph,3,0);
    SEQAN_ASSERT(vd == 20);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3);
    vd = findVertex(ali_graph,3,3);
    SEQAN_ASSERT(vd == 21);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 3);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,3,5);
    SEQAN_ASSERT(vd == 22);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 5);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);
    vd = findVertex(ali_graph,3,6);
    SEQAN_ASSERT(vd == 23);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 6);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2);
    vd = findVertex(ali_graph,3,8);
    SEQAN_ASSERT(vd == 24);
    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8);
    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1);

    SEQAN_ASSERT(findEdge(ali_graph,0,14)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,0,8)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,1,15)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,1,9)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,2,10)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,3,11)!=0);
//    SEQAN_ASSERT(findEdge(ali_graph,4,12)!=0)  //doesnt exist if edges with score <= 0 are kicked out
    SEQAN_ASSERT(findEdge(ali_graph,4,21)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,5,22)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,5,13)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,6,23)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,7,24)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,8,14)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,9,20)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,9,15)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,11,22)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,12,23)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,13,24)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,17,22)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,18,23)!=0);
    SEQAN_ASSERT(findEdge(ali_graph,19,24)!=0);


    //clear(ali_graph);
    //assignStringSet(ali_graph,seq_set);

    //matchRefinement(alis,seq_set,score_type,ali_graph,StoreEdges());
    //std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
    //std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
    //std::cout << ali_graph << "\n\n";


    //fstream strmW1; // Write the library
    //strmW1.open(TEST_PATH "my_testlib2alledges.lib", ios_base::out | ios_base::trunc);
    //write(strmW1,ali_graph,TCoffeeLib());
    //strmW1.close();


}




////produce pairwise alignments (Graph<Alignment>)
//template<typename TAlign, typename TStringSet, typename TScore>
//void
//getGraphAlignments(String<TAlign> & alis, TStringSet & seq, TScore & score_type, int & numAlignments, int cutoff)
//{
//    typedef StringSet<typename Value<TStringSet>::Type, Dependent<> > TAliStringSet;
//
//    int gesamt = 0;
//
//    for(int i = 0; i < length(seq); ++i)
//    {
//        for(int j = i+1; j < length(seq); ++j)
//        {
//            TAliStringSet str;
//            assignValueById(str, seq[i],positionToId(seq, i));
//            assignValueById(str, seq[j],positionToId(seq, j));
//            TAlign ali_g(str);
//            typename Value<TScore>::Type score = localAlignment(ali_g, score_type, WatermanEggert());
//            if(score==0)
//                continue;
//             int k = 1;
//            while(k<numAlignments)
//            {
//                score = localAlignment(ali_g, score_type, WatermanEggert());
//                if(score==0) k = numAlignments;
//                else ++k;
//            }
//            appendValue(alis,ali_g);
//            cout << ali_g <<"\n";
//            ++gesamt;
//        }
//    }
//
//    numAlignments = gesamt;
//    resize(alis,numAlignments);
//
//
//}



//
//
//void
//Test_RefineAlignGraph(){
//
//    typedef String<char> TString;
//    typedef StringSet<TString> TStringSet;
//    //typedef Align<typename Reference<TStringSet>::Type, ArrayGaps> TAlign;
//    //typedef Graph<Alignment<TStringSet, unsigned int> > TAlign;
//    typedef Graph<Alignment<StringSet<TString,Dependent<> >, unsigned int> > TAlign;
//
//
//    TStringSet seq_set;
//
//
//    TString str = "GARFIELDTHELASTFATCAT";
//    //appendValue(seq_set,str);
//    assignValueById(seq_set,str);
//
//    str = "GARFIELDTHEFASTCAT";
//    //appendValue(seq_set,str);
//    assignValueById(seq_set,str);
//
//    str = "GARFIELDTHEVERYFASTCAT";
//    //appendValue(seq_set,str);
//    assignValueById(seq_set,str);
//
//    str = "THEFATCAT";
//    //appendValue(seq_set,str);
//    assignValueById(seq_set,str);
//
//    int numSequences = length(seq_set);
//
//
//    int numAlignments = 2;
//    int numSequencePairs = 0;
//    int cutoff = 3;
//    for(int i = 1 ; i < numSequences; ++i)
//        numSequencePairs += i;
//    String<TAlign> alis;
//    reserve(alis,numSequencePairs*numAlignments);
//    Score<int> score_type = Score<int>(1,-1,-2,-2) ;
//
//    getGraphAlignments(alis,seq_set,score_type,numAlignments,cutoff);
//
//    typedef Graph<Alignment<TStringSet> > TAliGraph;
//    //TAliGraph ali_graph;
//    TAliGraph ali_graph(seq_set);
//
//    matchRefinement(alis,seq_set,score_type,ali_graph);
//
//    cout << ali_graph << "\n";
//    VertexDescriptor<TAliGraph>::Type vd;
//
//    vd = findVertex(ali_graph,0,0);
//    SEQAN_ASSERT(vd == 0)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8)
//    vd = findVertex(ali_graph,0,8);
//    SEQAN_ASSERT(vd == 1)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,0,11);
//    SEQAN_ASSERT(vd == 2)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,0,14);
//    SEQAN_ASSERT(vd == 3)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//    vd = findVertex(ali_graph,0,15);
//    SEQAN_ASSERT(vd == 4)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,0,17);
//    SEQAN_ASSERT(vd == 5)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//    vd = findVertex(ali_graph,0,18);
//    SEQAN_ASSERT(vd == 6)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,0,20);
//    SEQAN_ASSERT(vd == 7)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//
//    vd = findVertex(ali_graph,1,0);
//    SEQAN_ASSERT(vd == 8)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8)
//    vd = findVertex(ali_graph,1,8);
//    SEQAN_ASSERT(vd == 9)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,1,11);
//    SEQAN_ASSERT(vd == 10)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 11)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,1,14);
//    SEQAN_ASSERT(vd == 11)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 14)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//    vd = findVertex(ali_graph,1,15);
//    SEQAN_ASSERT(vd == 12)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 15)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,1,17);
//    SEQAN_ASSERT(vd == 13)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 17)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//
//    vd = findVertex(ali_graph,2,0);
//    SEQAN_ASSERT(vd == 14)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 8)
//    vd = findVertex(ali_graph,2,8);
//    SEQAN_ASSERT(vd == 15)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,2,11);
//    SEQAN_ASSERT(vd == 16)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 11)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 7)
//    vd = findVertex(ali_graph,2,18);
//    SEQAN_ASSERT(vd == 17)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 18)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//    vd = findVertex(ali_graph,2,19);
//    SEQAN_ASSERT(vd == 18)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 19)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,2,21);
//    SEQAN_ASSERT(vd == 19)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 21)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//
//    vd = findVertex(ali_graph,3,0);
//    SEQAN_ASSERT(vd == 20)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 0)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 3)
//    vd = findVertex(ali_graph,3,3);
//    SEQAN_ASSERT(vd == 21)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 3)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,3,5);
//    SEQAN_ASSERT(vd == 22)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 5)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//    vd = findVertex(ali_graph,3,6);
//    SEQAN_ASSERT(vd == 23)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 6)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 2)
//    vd = findVertex(ali_graph,3,8);
//    SEQAN_ASSERT(vd == 24)
//    SEQAN_ASSERT(fragmentBegin(ali_graph,vd) == 8)
//    SEQAN_ASSERT(fragmentLength(ali_graph,vd) == 1)
//
//    SEQAN_ASSERT(findEdge(ali_graph,0,14)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,0,8)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,1,15)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,1,9)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,2,10)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,3,11)!=0)
////    SEQAN_ASSERT(findEdge(ali_graph,4,12)!=0)  //doesnt exist if edges with score <= 0 are kicked out
//    SEQAN_ASSERT(findEdge(ali_graph,4,21)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,5,22)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,5,13)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,6,23)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,7,24)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,8,14)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,9,20)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,9,15)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,11,22)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,12,23)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,13,24)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,17,22)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,18,23)!=0)
//    SEQAN_ASSERT(findEdge(ali_graph,19,24)!=0)
//
//
//}


// TODO(holtgrew): This seems to be an artefact from debugging. Remove?
SEQAN_DEFINE_TEST(GraphMatchRefinement_Problem)
{
    typedef String<char> TString;
    typedef StringSet<TString> TStringSet;
    //typedef StringSet<TString, Dependent<> > TAlignmentStringSet;
    //typedef Graph<Alignment<TStringSet> > TAlign;
    typedef Fragment<> TFragment;
    typedef String<TFragment> TFragString;


    TString str1 = "RKNLVQFGVGEKNGSVRWVMNALGVKDDWLLVPSHAYKFEKDYEMMEFYFNRGGTYYSISAGNVVIQSLDVGFQDVVLMKVPTIPKFRDITQHFIKKGDVPRALNRLATLVTTVNGTPMLISEGPLKMEEKATYVHKKNDGTTVDLTVDQAWRGKGEGLPGMCGGALVSSNQSIQNAILGIHVAGGNSILVAKLVTQE";
    TString str2 = "";
    TString str3 = "";
    TString str4 = "IAGGEAITTGGSRCSLGFNVVAHALTAGHCTNISAWSIGTRTGTSFNNDYGIIRHSNPAAADGRVYLYQDITTAGNAFVGQAVQRSGSTTGLRSGSVTGLNATVNYGSSGIVYGMIQTNVCAGDSGGSLFAGSTALGLTSGGSGNCRTGGTTFYQPVT";
    TString str5 = "";

    TStringSet strSet;
    assignValueById(strSet,str1,0u);
    assignValueById(strSet,str2,1u);
    assignValueById(strSet,str3,2u);
    assignValueById(strSet,str4,3u);
    assignValueById(strSet,str5,4u);
    //cout << length(strSet)<<"\n";
    //cout << idToPosition(strSet,0) <<"\n";
    //cout << idToPosition(strSet,3) <<"\n";
    //cout << strSet[0] <<"\n";
    //cout << strSet[1] <<"\n";
    //cout << strSet[2] <<"\n";
    //cout << strSet[3] <<"\n";
    //cout << positionToId(strSet,0) <<"\n";
    //cout << positionToId(strSet,1) <<"\n";


    TFragment frag(0,28,3,35,18);

//    typedef String<Fragment<>, External<ExternalConfig<File<>, 64*1024> > > TFragString;

    TFragString alis;
    //resize(alis,1);
    //alis[0] = frag;
    //appendValue(alis, (TFragment) frag);

    //TAlign outGraph(strSet);
    //matchRefinement(alis,strSet,outGraph);

}


void Test_GraphMatchRefinement()
{

//    Test_Problem();

    //test for graph_align on Align<TSource,TSpec>
//    Test_RefineAlign();

    //test for graph_align on Graph<Alignment<> >
//    Test_RefineAlignGraph();

    //test for graph_align on Fragment<>
//    Test_GraphMatchRefine();
}

}

#endif

