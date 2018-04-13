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


// External / STL
#include <sstream>

// SeqAn Includes
#include <seqan/basic.h>

#ifndef SEQAN_HEADER_TEST_GRAPH_IMPL_ALIGN_H
#define SEQAN_HEADER_TEST_GRAPH_IMPL_ALIGN_H


namespace seqan {

//////////////////////////////////////////////////////////////////////////////

// Alignment without edge weights
SEQAN_DEFINE_TEST(Test_Refinement_AlignmentGraphNoEdgeWeights)
{
    typedef String<Dna> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef    Id<TStringSet>::Type TId;
    TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

    // Difference between Owner and Dependent
    typedef StringSet<TString, Owner<> > TOwnerStringSet;
    TStringSet str;
    TOwnerStringSet ownStr;
    TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
    TId id0 = assignValueById(str, str0);
    appendValue(ownStr, str0);
    TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
    TId id1 = assignValueById(str, str1);
    appendValue(ownStr, str1);
    TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
    TId id2 = assignValueById(str, str2);
    appendValue(ownStr, str2);

    // Check that the graph makes a dependent StringSet
    TGraph alwaysDependStringSet(ownStr);
    SEQAN_ASSERT_EQ(length(value(stringSet(alwaysDependStringSet), 0)), length(str0));
    clear(alwaysDependStringSet);
    SEQAN_ASSERT_EQ(length(value(ownStr, 0)), length(str0));
    assignStringSet(alwaysDependStringSet, ownStr);
    SEQAN_ASSERT_EQ(length(value(stringSet(alwaysDependStringSet), 0)), length(str0));
    clear(ownStr);
    SEQAN_ASSERT(empty(ownStr));


    TGraph g(str);
    SEQAN_ASSERT_EQ(getStringSet(g)[0], str0);
    SEQAN_ASSERT_EQ(getStringSet(g)[1], str1);
    SEQAN_ASSERT_EQ(stringSet(g)[2], str2);
    assignStringSet(g, str);
    SEQAN_ASSERT_EQ(getStringSet(g)[0], str0);
    SEQAN_ASSERT_EQ(getStringSet(g)[1], str1);
    SEQAN_ASSERT_EQ(stringSet(g)[2], str2);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    SEQAN_ASSERT_EQ(numVertices(g), 0u);
    SEQAN_ASSERT(empty(g));

    TVertexDescriptor v0 = addVertex(g, id1, 0, 2);
    SEQAN_ASSERT_EQ(v0, 0u);
    SEQAN_ASSERT_EQ(outDegree(g, v0), 0u);
    SEQAN_ASSERT_EQ(inDegree(g, 0), 0u);
    SEQAN_ASSERT_EQ(degree(g, 0), 0u);
    SEQAN_ASSERT_EQ(numVertices(g), 1u);
    SEQAN_ASSERT(!empty(g));
    SEQAN_ASSERT_EQ(sequenceId(g, v0), id1);
    SEQAN_ASSERT_EQ(label(g, v0), "cc");
    SEQAN_ASSERT_EQ(fragmentBegin(g, v0), 0u);
    SEQAN_ASSERT_EQ(fragmentLength(g, v0), 2u);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 0), v0);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 1), v0);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 2), nilVertex);

    TVertexDescriptor v1 = addVertex(g, id2, 0, 5);
    TEdgeDescriptor e = addEdge(g,0,1);
    SEQAN_ASSERT_EQ(_getVertexString(g)[0], e);
    SEQAN_ASSERT_EQ(getIdUpperBound(_getVertexIdManager(g)), 2u);
    SEQAN_ASSERT_EQ(getIdUpperBound(_getEdgeIdManager(g)), 1u);
    SEQAN_ASSERT_EQ(v1, 1u);
    SEQAN_ASSERT_EQ(numVertices(g), 2u);
    SEQAN_ASSERT_EQ(targetVertex(g, e), 1u);
    SEQAN_ASSERT_EQ(sourceVertex(g, e), 0u);
    SEQAN_ASSERT_EQ(numEdges(g), 1u);
    SEQAN_ASSERT_EQ(outDegree(g, v0), 1u);
    SEQAN_ASSERT_EQ(inDegree(g, 1), 1u);
    SEQAN_ASSERT_EQ(inDegree(g, 0), 1u);
    SEQAN_ASSERT_EQ(degree(g, 0), 1u);
    SEQAN_ASSERT_EQ(findVertex(g, id2, 0), v1);
    SEQAN_ASSERT_EQ(findVertex(g, id2, 1), v1);
    SEQAN_ASSERT_EQ(findVertex(g, id2, 4), v1);
    SEQAN_ASSERT_EQ(findVertex(g, id2, 5), nilVertex);

    // Add more vertices and edges
    addVertex(g, id1, 10, 20);  // 2
    SEQAN_ASSERT_EQ(findVertex(g, id1, 0), v0);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 1), v0);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 2), nilVertex);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 10), 2u);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 19), 2u);
    SEQAN_ASSERT_EQ(findVertex(g, id1, 30), nilVertex);
    TVertexDescriptor v3 = addVertex(g, id2, 5, 2);  // 3
    addVertex(g, id1, 7, 3);  // 4
    addEdge(g, 3, 4);
    TEdgeDescriptor my_edge = addEdge(g, 3, 1);
    addEdge(g, 3, 0);
    SEQAN_ASSERT_EQ(v3, 3u);
    SEQAN_ASSERT_EQ(numVertices(g), 5u);
    SEQAN_ASSERT_EQ(targetVertex(g, my_edge), 3u);
    SEQAN_ASSERT_EQ(sourceVertex(g, my_edge), 1u);
    SEQAN_ASSERT_EQ(numEdges(g), 4u);
    SEQAN_ASSERT_EQ(outDegree(g, v3), 3u);
    SEQAN_ASSERT_EQ(inDegree(g, v3), 3u);
    SEQAN_ASSERT_EQ(degree(g, v3), 3u);

    // Remove edges
    removeEdge(g, 3, 1);
    removeEdge(g, 0, 1);
    SEQAN_ASSERT_EQ(numEdges(g), 2u);

    // Remove vertices
    addVertex(g, id2, 14, 4);  // 5
    addEdge(g, 5, 2);
    addEdge(g, 2, 3);
    addEdge(g, 1, 3);
    addEdge(g, 1, 4);
    SEQAN_ASSERT_EQ(outDegree(g, 3), 4u);
    SEQAN_ASSERT_EQ(outDegree(g, 4), 2u);
    removeVertex(g, v3);
    SEQAN_ASSERT_EQ(outDegree(g, 4), 1u);
    SEQAN_ASSERT_EQ(outDegree(g, 0), 0u);
    SEQAN_ASSERT_EQ(numVertices(g), 5u);
    SEQAN_ASSERT_EQ(numEdges(g), 2u);

    // Clear graph
    clearEdges(g);
    SEQAN_ASSERT_EQ(numVertices(g), 5u);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clearVertices(g);
    SEQAN_ASSERT_EQ(numVertices(g), 0u);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    addVertex(g, id1, 0, 1); addVertex(g, id1, 1, 1); addVertex(g, id1, 2, 1);
    addVertex(g, id1, 3, 1); addVertex(g, id1, 4, 1);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clear(g);
    assignStringSet(g, str);
    SEQAN_ASSERT_EQ(numVertices(g), 0u);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    addVertex(g, id1, 0, 1); addVertex(g, id1, 1, 1); addVertex(g, id1, 2, 1);
    addVertex(g, id1, 3, 1); addVertex(g, id1, 4, 1);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    addEdge(g, 4, 2);
    removeVertex(g, 3);
    SEQAN_ASSERT_EQ(numVertices(g), 4u);
    SEQAN_ASSERT_EQ(numEdges(g), 3u);
    SEQAN_ASSERT_EQ(outDegree(g, 4), 2u);
    SEQAN_ASSERT_EQ(inDegree(g, 4), 2u);

    // Transpose
    transpose(g);
    SEQAN_ASSERT_EQ(numVertices(g), 4u);
    SEQAN_ASSERT_EQ(numEdges(g), 3u);
    SEQAN_ASSERT_EQ(outDegree(g, 4), 2u);
    SEQAN_ASSERT_EQ(inDegree(g, 4), 2u);
    TGraph g_copy(g);
    SEQAN_ASSERT_EQ(numVertices(g_copy), 4u);
    SEQAN_ASSERT_EQ(numEdges(g_copy), 3u);
    SEQAN_ASSERT_EQ(outDegree(g_copy, 4), 2u);
    SEQAN_ASSERT_EQ(inDegree(g_copy, 4), 2u);
    addVertex(g_copy, id0, 0, 3);
    addEdge(g_copy, 3, 0);
    g_copy = g;
    SEQAN_ASSERT_EQ(numVertices(g_copy), 4u);
    SEQAN_ASSERT_EQ(numEdges(g_copy), 3u);
    SEQAN_ASSERT_EQ(outDegree(g_copy, 4), 2u);
    SEQAN_ASSERT_EQ(inDegree(g_copy, 4), 2u);
    //Copies the graph and transposes just the copy
    transpose(g, g_copy);  // g does not change!
    SEQAN_ASSERT_EQ(numVertices(g_copy), 4u);
    SEQAN_ASSERT_EQ(numEdges(g_copy), 3u);
    SEQAN_ASSERT_EQ(outDegree(g_copy, 4), 2u);
    SEQAN_ASSERT_EQ(inDegree(g_copy, 4), 2u);

    // Adjacency matrix
    String<unsigned int> mat;
    getAdjacencyMatrix(g, mat);
    unsigned int len = (unsigned int) std::sqrt((double) length(mat));
    SEQAN_ASSERT_EQ(getValue(mat, 0 * len + 2), 1u);
    SEQAN_ASSERT_EQ(getValue(mat, 3 * len + 2), 0u);
    SEQAN_ASSERT_EQ(getValue(mat, 0 * len + 2), getValue(mat, 2 * len + 0));
    SEQAN_ASSERT_EQ(getValue(mat, 1 * len + 4), getValue(mat, 4 * len + 1));
    SEQAN_ASSERT_EQ(getValue(mat, 2 * len + 4), getValue(mat, 4 * len + 2));

    // Vertex Adjacency vectors
    String<unsigned int> vectIn, vectOut;
    getVertexAdjacencyVector(vectIn, vectOut, g, 2u);
    SEQAN_ASSERT(length(vectIn) == 2);
    SEQAN_ASSERT(length(vectOut) == 2);
    SEQAN_ASSERT(vectIn[0] == 4);
    SEQAN_ASSERT(vectIn[1] == 0);
    SEQAN_ASSERT(vectOut[0] == 4);
    SEQAN_ASSERT(vectOut[1] == 0);
}

//////////////////////////////////////////////////////////////////////////////

// Alignments with edge weights
SEQAN_DEFINE_TEST(Test_Refinement_AlignmentGraphEdgeWeights)
{
    typedef String<Dna> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef    Id<TStringSet>::Type TId;

    typedef Graph<Alignment<TStringSet> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

    TStringSet str;
    TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
    appendValue(str, str0);
    TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
    TId id1 = assignValueById(str, str1);
    TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
    TId id2 = assignValueById(str, str2);

    TGraph g(str);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    SEQAN_ASSERT_EQ(numVertices(g), 0u);
    SEQAN_ASSERT(empty(g));

    TVertexDescriptor v0 = addVertex(g, id1, 0, 2);
    TVertexDescriptor v1 = addVertex(g, id2, 0, 5);
    SEQAN_ASSERT_EQ(v0, 0u);
    SEQAN_ASSERT_EQ(v1, 1u);
    TEdgeDescriptor e1 = addEdge(g, 0, 1, 100);
    SEQAN_ASSERT_EQ(getCargo(e1), 100u);
    SEQAN_ASSERT_EQ(numEdges(g), 1u);
    removeEdge(g, e1);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    e1 = addEdge(g, 0, 1, 1005);
    SEQAN_ASSERT_EQ(numEdges(g), 1u);
    removeOutEdges(g, 0);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
    e1 = addEdge(g, 0, 1, 1005);
    SEQAN_ASSERT_EQ(findEdge(g, 0, 1), e1);
    SEQAN_ASSERT_EQ(numEdges(g), 1u);
    removeInEdges(g, 0);
    SEQAN_ASSERT_EQ(numEdges(g), 0u);
}

//////////////////////////////////////////////////////////////////////////////

// Alignment Graph Iterators
SEQAN_DEFINE_TEST(Test_Refinement_AlignmentGraphIterators)
{
    typedef String<Dna> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef    Id<TStringSet>::Type TId;

    typedef Graph<Alignment<TStringSet, void> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TStringSet str;
    TString str0("acaagtaacataaaaaaaaaaaaaaaacccccccccttttttttaaaaa");
    appendValue(str, str0);
    TString str1("cccaaagggtttttccccccccccccttttttttttaaaaaaagggggggg");
    TId id1 = assignValueById(str, str1);
    TString str2("cacatgtaatcatgggggggggccccccttttaaaaaaaaaaatttt");
    appendValue(str, str2);

    TGraph g(str);
    TVertexDescriptor v0 = addVertex(g, id1, 0, 1);
    addVertex(g, id1, 1, 1); addVertex(g, id1, 2, 1);
    addVertex(g, id1, 3, 1); addVertex(g, id1, 4, 1);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    addEdge(g, 4, 2);
    removeVertex(g, 3);

    // Vertex Iterator
    typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    SEQAN_ASSERT(atBegin(itV));
    SEQAN_ASSERT_EQ(getValue(itV), 0u);
    SEQAN_ASSERT_EQ(value(itV), 0u);
    SEQAN_ASSERT_EQ(*itV, 0u);
    goNext(itV);
    SEQAN_ASSERT(!atBegin(itV));
    SEQAN_ASSERT_EQ(getValue(itV), 1u);
    ++itV;
    SEQAN_ASSERT_EQ(getValue(itV), 2u);
    SEQAN_ASSERT(!atEnd(itV));
    goPrevious(itV);
    SEQAN_ASSERT_EQ(*itV, 1u);
    SEQAN_ASSERT_NOT(atEnd(itV));
    itV--;
    SEQAN_ASSERT_EQ(getValue(itV), 0u);
    SEQAN_ASSERT(atBegin(itV));

    // OutEdge Iterator
    typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    TOutEdgeIterator it(g, v0);
    SEQAN_ASSERT_EQ(sourceVertex(g, getValue(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, getValue(it)), 2u);
    SEQAN_ASSERT_EQ(sourceVertex(g, value(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, *it), 2u);
    SEQAN_ASSERT_NOT(atEnd(it));
    SEQAN_ASSERT(atBegin(it));
    goNext(it);
    SEQAN_ASSERT(atEnd(it));
    SEQAN_ASSERT_NOT(atBegin(it));
    goPrevious(it);
    SEQAN_ASSERT_EQ(sourceVertex(g, getValue(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, getValue(it)), 2u);
    --it;
    it--;
    SEQAN_ASSERT(atBegin(it));
    TOutEdgeIterator it2(g, v0);
    TOutEdgeIterator it3;
    it3 = it;
    SEQAN_ASSERT(it == it2);
    SEQAN_ASSERT(it2 == it3);
    goEnd(it);
    SEQAN_ASSERT(it2 != it);
    goEnd(it2);
    SEQAN_ASSERT(it2 == it);
    goBegin(it2);
    SEQAN_ASSERT(it2 != it);
    SEQAN_ASSERT_EQ(&g, &hostGraph(it));

    // EdgeIterator
    typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEdge(g);
    SEQAN_ASSERT_EQ(sourceVertex(g, getValue(itEdge)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, getValue(itEdge)), 2u);
    SEQAN_ASSERT_EQ(sourceVertex(g, value(itEdge)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, *itEdge), 2u);
    SEQAN_ASSERT_NOT(atEnd(itEdge));
    SEQAN_ASSERT(atBegin(itEdge));
    goNext(itEdge);
    SEQAN_ASSERT_EQ(sourceVertex(g, value(itEdge)), 1u);
    SEQAN_ASSERT_EQ(targetVertex(g, *itEdge), 4u);
    ++itEdge;
    --itEdge;
    SEQAN_ASSERT_EQ(sourceVertex(g, value(itEdge)), 1u);
    SEQAN_ASSERT_EQ(targetVertex(g, *itEdge), 4u);
    goEnd(itEdge);
    SEQAN_ASSERT(atEnd(itEdge));
    SEQAN_ASSERT_NOT(atBegin(itEdge));
    goBegin(itEdge);
    SEQAN_ASSERT_NOT(atEnd(itEdge));
    SEQAN_ASSERT(atBegin(itEdge));

    // Adjacency Iterator
    typedef Iterator<TGraph, AdjacencyIterator>::Type TAdjacencyIterator;
    TAdjacencyIterator itAdj(g, 2);
    SEQAN_ASSERT_EQ(getValue(itAdj), 4u);
    SEQAN_ASSERT_EQ(&hostGraph(itAdj),  &g);
    SEQAN_ASSERT_EQ(value(itAdj),  4u);
    SEQAN_ASSERT_EQ(*itAdj,  4u);
    SEQAN_ASSERT_NOT(atEnd(itAdj));
    SEQAN_ASSERT(atBegin(itAdj));
    goNext(itAdj);
    SEQAN_ASSERT_EQ(*itAdj, 0u);
    SEQAN_ASSERT_NOT(atEnd(itAdj));
    SEQAN_ASSERT_NOT(atBegin(itAdj));
    ++itAdj;
    SEQAN_ASSERT(atEnd(itAdj));
    SEQAN_ASSERT_NOT(atBegin(itAdj));
    goPrevious(itAdj);
    --itAdj;
    SEQAN_ASSERT_EQ(*itAdj, 4u);
    goBegin(itAdj);
    SEQAN_ASSERT(atBegin(itAdj));
    goEnd(itAdj);
    SEQAN_ASSERT(atEnd(itAdj));

    // Bfs Iterator
    typedef Iterator<TGraph, BfsIterator>::Type TBfsIterator;
    TBfsIterator bfsIt(g, 2);
    SEQAN_ASSERT_NOT(atEnd(bfsIt));
    SEQAN_ASSERT(atBegin(bfsIt));
    ++bfsIt;
    SEQAN_ASSERT_EQ(getValue(bfsIt), 4u);
    SEQAN_ASSERT_EQ(&hostGraph(bfsIt), &g);
    SEQAN_ASSERT_EQ(value(bfsIt), 4u);
    SEQAN_ASSERT_EQ(*bfsIt, 4u);
    goNext(bfsIt);
    SEQAN_ASSERT_EQ(value(bfsIt), 0u);

    // Dfs Iterator
    typedef Iterator<TGraph, DfsPreorder>::Type TDfsPreorder;
    TDfsPreorder dfsIt(g, 2);
    SEQAN_ASSERT_NOT(atEnd(dfsIt));
    SEQAN_ASSERT(atBegin(dfsIt));
    SEQAN_ASSERT_EQ(*dfsIt, 2u);
    ++dfsIt;
    SEQAN_ASSERT_EQ(getValue(dfsIt), 0u);
    SEQAN_ASSERT_EQ(&hostGraph(dfsIt), &g);
    SEQAN_ASSERT_EQ(value(dfsIt), 0u);
    SEQAN_ASSERT_EQ(*dfsIt, 0u);
    goNext(dfsIt);
}

//////////////////////////////////////////////////////////////////////////////

// Output Alignments in Alignment Graph
SEQAN_DEFINE_TEST(Test_Refinement_AlignmentGraphOutput)
{
    // Alignments
    typedef String<char> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef    Id<TStringSet>::Type TId;

    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
    typedef VertexDescriptor<TAlignmentGraph>::Type TVD;
    //typedef EdgeDescriptor<TAlignmentGraph>::Type TED;

    TStringSet str;
    TString str0("Garfieldthelastfatcat");
    TId i0 = assignValueById(str, str0);
    TString str1("Garfieldthefastcat");
    TId i1 = assignValueById(str, str1);
    TString str2("Garfieldtheveryfastcat");
    TId i2 = assignValueById(str, str2);
    TString str3("thefatcat");
    TId i3 = assignValueById(str, str3);

    TAlignmentGraph g(str);
    TVD vH = addVertex(g, i1, 8, 3);
    SEQAN_ASSERT_EQ(getFirstCoveredPosition(g, i3), 9u);  // Not found, length of the sequence
    SEQAN_ASSERT_EQ(getLastCoveredPosition(g, i3), 0u);  // Not found, 0
    SEQAN_ASSERT_EQ(getFirstCoveredPosition(g, i1), 8u);
    SEQAN_ASSERT_EQ(getLastCoveredPosition(g, i1), 11u);
    TVD vT = addVertex(g, i1, 13, 1); TVD vS = addVertex(g, i3, 6, 3);
    TVD vW = addVertex(g, i2, 18, 1); TVD vA = addVertex(g, i0, 0, 8);
    SEQAN_ASSERT_EQ(getFirstCoveredPosition(g, i0), 0u);
    SEQAN_ASSERT_EQ(getLastCoveredPosition(g, i0), 8u);
    TVD vM = addVertex(g, i2, 11, 4); TVD vK = addVertex(g, i2, 0, 8);
    TVD vC = addVertex(g, i0, 11, 4); TVD vD = addVertex(g, i0, 15, 2);
    TVD vF = addVertex(g, i0, 18, 3); TVD vG = addVertex(g, i1, 0, 8);
    addEdge(g, vA, vG);
    TVD vI = addVertex(g, i1, 11, 2); TVD vQ = addVertex(g, i3, 3, 2); TVD vB = addVertex(g, i0, 8, 3);
    TVD vU = addVertex(g, i1, 14, 1); TVD vE = addVertex(g, i0, 17, 1); TVD vJ = addVertex(g, i1, 15, 3);
    TVD vL = addVertex(g, i2, 8, 3);
    addEdge(g, vH, vL);
    TVD vN = addVertex(g, i2, 15, 2); TVD vV = addVertex(g, i2, 17, 1);
    TVD vO = addVertex(g, i2, 19, 3); TVD vP = addVertex(g, i3, 0, 3); TVD vR = addVertex(g, i3, 5, 1);
    addEdge(g, vA, vK); addEdge(g, vG, vK); addEdge(g, vB, vH); addEdge(g, vB, vL);
    addEdge(g, vB, vP); addEdge(g, vH, vP); addEdge(g, vL, vP); addEdge(g, vC, vM);
    addEdge(g, vD, vI); addEdge(g, vD, vQ); addEdge(g, vD, vN); addEdge(g, vI, vQ);
    addEdge(g, vI, vN); addEdge(g, vQ, vN); addEdge(g, vT, vV); addEdge(g, vE, vU);
    addEdge(g, vE, vW); addEdge(g, vE, vR); addEdge(g, vU, vW); addEdge(g, vU, vR);
    addEdge(g, vW, vR); addEdge(g, vF, vJ); addEdge(g, vF, vO); addEdge(g, vF, vS);
    addEdge(g, vJ, vO); addEdge(g, vJ, vS); addEdge(g, vO, vS);
    SEQAN_ASSERT_EQ(getFirstCoveredPosition(g, i3), 0u);
    SEQAN_ASSERT_EQ(getLastCoveredPosition(g, i3), 9u);

    // Output of the alignment graph
    std::stringstream sstream, expected;

    // standard output
    sstream << g;
    expected << "Alignment matrix:\n"
             << "      0     .    :    .    :   \n"
             << "        Garfieldthelastfa-tcat\n"
             << "        |||||||||||    || ||||\n"
             << "        Garfieldthe----fastcat\n"
             << "        |||||||||||    |||||||\n"
             << "        Garfieldtheveryfastcat\n"
             << "                |||    || ||||\n"
             << "        --------the----fa-tcat\n\n\n";
    SEQAN_ASSERT_EQ(sstream.str(), expected.str());
    sstream.str(""); sstream.clear();
    expected.str(""); expected.clear();

    StringSet<String<char> > seqs;
    appendValue(seqs, "seq1");
    appendValue(seqs, "seq2");
    appendValue(seqs, "seq3");
    appendValue(seqs, "seq4");

    // fasta output
    write(sstream, g, seqs, FastaFormat());
    expected << ">seq1\nGarfieldthelastfa-tcat\n"
             << ">seq2\nGarfieldthe----fastcat\n"
             << ">seq3\nGarfieldtheveryfastcat\n"
             << ">seq4\n--------the----fa-tcat\n";
    SEQAN_ASSERT_EQ(sstream.str(), expected.str());
    sstream.str(""); sstream.clear();
    expected.str(""); expected.clear();

    // msf output
    write(sstream, g, seqs, MsfFormat());
    expected << "PileUp\n\n"
             << " MSF: 22 Type: P Check: 0 ..\n\n"
             << " Name: seq1 oo  Len:  22 Check: 0 Weight: 1.00\n"
             << " Name: seq2 oo  Len:  22 Check: 0 Weight: 1.00\n"
             << " Name: seq3 oo  Len:  22 Check: 0 Weight: 1.00\n"
             << " Name: seq4 oo  Len:  22 Check: 0 Weight: 1.00\n\n"
             << "//\n\n\n"
             << "seq1      Garfieldth elastfa.tc at\n"
             << "seq2      Garfieldth e....fastc at\n"
             << "seq3      Garfieldth everyfastc at\n"
             << "seq4      ........th e....fa.tc at\n\n\n";
    SEQAN_ASSERT_EQ(sstream.str(), expected.str());
    sstream.str(""); sstream.clear();
    expected.str(""); expected.clear();

    // cg viz output
    write(sstream, g, seqs, CgVizFormat());
    expected << "{DATA Data\n"
             << "\t[__GLOBAL__] tracks=4\n"
             << "\tfasta_id=\"seq1\" sequence=\"Garfieldthelastfatcat\" track=0 type=\"DNA\": 0 20\n"
             << "\tfasta_id=\"seq2\" sequence=\"Garfieldthefastcat\" track=1 type=\"DNA\": 0 17\n"
             << "\tfasta_id=\"seq3\" sequence=\"Garfieldtheveryfastcat\" track=2 type=\"DNA\": 0 21\n"
             << "\tfasta_id=\"seq4\" sequence=\"thefatcat\" track=3 type=\"DNA\": 0 8\n"
             << "}\n"
             << "{DATA 0-seqlen\n"
             << "\t[__GLOBAL__]\n"
             << "\tlength=21:\t0 20\n"
             << "}\n"
             << "{DATA 1-seqlen\n"
             << "\t[__GLOBAL__]\n"
             << "\tlength=18:\t0 17\n"
             << "}\n"
             << "{DATA 2-seqlen\n"
             << "\t[__GLOBAL__]\n"
             << "\tlength=22:\t0 21\n"
             << "}\n"
             << "{DATA 3-seqlen\n"
             << "\t[__GLOBAL__]\n"
             << "\tlength=9:\t0 8\n"
             << "}\n"
             << "{DATA 0-vs-1\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=0 target=13 edgeId=4 cargo=0 label=the labelOpp=the:\t8 8 11 11\n"
             << "\tsource=4 target=10 edgeId=0 cargo=0 label=Garfield labelOpp=Garfield:\t0 0 8 8\n"
             << "\tsource=8 target=11 edgeId=10 cargo=0 label=fa labelOpp=fa:\t15 11 17 13\n"
             << "\tsource=9 target=16 edgeId=23 cargo=0 label=cat labelOpp=cat:\t18 15 21 18\n"
             << "\tsource=14 target=15 edgeId=17 cargo=0 label=t labelOpp=t:\t14 17 15 18\n"
             << "}\n"
             << "{DATA 0-vs-2\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=3 target=15 edgeId=18 cargo=0 label=t labelOpp=t:\t18 17 19 18\n"
             << "\tsource=4 target=6 edgeId=2 cargo=0 label=Garfield labelOpp=Garfield:\t0 0 8 8\n"
             << "\tsource=5 target=7 edgeId=9 cargo=0 label=very labelOpp=last:\t11 11 15 15\n"
             << "\tsource=8 target=18 edgeId=12 cargo=0 label=fa labelOpp=fa:\t15 15 17 17\n"
             << "\tsource=9 target=20 edgeId=24 cargo=0 label=cat labelOpp=cat:\t18 19 21 22\n"
             << "\tsource=13 target=17 edgeId=5 cargo=0 label=the labelOpp=the:\t8 8 11 11\n"
             << "}\n"
             << "{DATA 0-vs-3\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=2 target=9 edgeId=25 cargo=0 label=cat labelOpp=cat:\t6 18 9 21\n"
             << "\tsource=8 target=12 edgeId=11 cargo=0 label=fa labelOpp=fa:\t15 3 17 5\n"
             << "\tsource=13 target=21 edgeId=6 cargo=0 label=the labelOpp=the:\t8 0 11 3\n"
             << "\tsource=15 target=22 edgeId=19 cargo=0 label=t labelOpp=t:\t17 5 18 6\n"
             << "}\n"
             << "{DATA 1-vs-2\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=0 target=17 edgeId=1 cargo=0 label=the labelOpp=the:\t8 8 11 11\n"
             << "\tsource=1 target=19 edgeId=16 cargo=0 label=s labelOpp=s:\t13 17 14 18\n"
             << "\tsource=3 target=14 edgeId=20 cargo=0 label=t labelOpp=t:\t18 14 19 15\n"
             << "\tsource=6 target=10 edgeId=3 cargo=0 label=Garfield labelOpp=Garfield:\t0 0 8 8\n"
             << "\tsource=11 target=18 edgeId=14 cargo=0 label=fa labelOpp=fa:\t11 15 13 17\n"
             << "\tsource=16 target=20 edgeId=26 cargo=0 label=cat labelOpp=cat:\t15 19 18 22\n"
             << "}\n"
             << "{DATA 1-vs-3\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=0 target=21 edgeId=7 cargo=0 label=the labelOpp=the:\t8 0 11 3\n"
             << "\tsource=2 target=16 edgeId=27 cargo=0 label=cat labelOpp=cat:\t6 15 9 18\n"
             << "\tsource=11 target=12 edgeId=13 cargo=0 label=fa labelOpp=fa:\t11 3 13 5\n"
             << "\tsource=14 target=22 edgeId=21 cargo=0 label=t labelOpp=t:\t14 5 15 6\n"
             << "}\n"
             << "{DATA 2-vs-3\n"
             << "\t[__GLOBAL__]\n"
             << "\tsource=2 target=20 edgeId=28 cargo=0 label=cat labelOpp=cat:\t6 19 9 22\n"
             << "\tsource=3 target=22 edgeId=22 cargo=0 label=t labelOpp=t:\t18 5 19 6\n"
             << "\tsource=12 target=18 edgeId=15 cargo=0 label=fa labelOpp=fa:\t3 15 5 17\n"
             << "\tsource=17 target=21 edgeId=8 cargo=0 label=the labelOpp=the:\t8 0 11 3\n"
             << "}\n";
    SEQAN_ASSERT_EQ(sstream.str(), expected.str());

    // GraphViz DOT output.
    sstream.str("");
    expected.str("");
    writeRecords(sstream, g, DotDrawing());
    expected << "graph G {\n"
             << "\n"
             << "/* Graph Attributes */\n"
             << "graph [rankdir = LR];\n"
             << "\n"
             << "/* Node Attributes */\n"
             << "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n"
             << "\n"
             << "/* Edge Attributes */\n"
             << "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n"
             << "\n"
             << "/* Nodes */\n"
             << "0 [label = \"[8,11)\", group = 1];\n"
             << "1 [label = \"[13,14)\", group = 1];\n"
             << "2 [label = \"[6,9)\", group = 3];\n"
             << "3 [label = \"[18,19)\", group = 2];\n"
             << "4 [label = \"[0,8)\", group = 0];\n"
             << "5 [label = \"[11,15)\", group = 2];\n"
             << "6 [label = \"[0,8)\", group = 2];\n"
             << "7 [label = \"[11,15)\", group = 0];\n"
             << "8 [label = \"[15,17)\", group = 0];\n"
             << "9 [label = \"[18,21)\", group = 0];\n"
             << "10 [label = \"[0,8)\", group = 1];\n"
             << "11 [label = \"[11,13)\", group = 1];\n"
             << "12 [label = \"[3,5)\", group = 3];\n"
             << "13 [label = \"[8,11)\", group = 0];\n"
             << "14 [label = \"[14,15)\", group = 1];\n"
             << "15 [label = \"[17,18)\", group = 0];\n"
             << "16 [label = \"[15,18)\", group = 1];\n"
             << "17 [label = \"[8,11)\", group = 2];\n"
             << "18 [label = \"[15,17)\", group = 2];\n"
             << "19 [label = \"[17,18)\", group = 2];\n"
             << "20 [label = \"[19,22)\", group = 2];\n"
             << "21 [label = \"[0,3)\", group = 3];\n"
             << "22 [label = \"[5,6)\", group = 3];\n"
             << "\n"
             << "/* Edges */\n"
             << "0 -- 21 [];\n"
             << "0 -- 13 [];\n"
             << "0 -- 17 [];\n"
             << "1 -- 19 [];\n"
             << "2 -- 20 [];\n"
             << "2 -- 16 [];\n"
             << "2 -- 9 [];\n"
             << "3 -- 22 [];\n"
             << "3 -- 14 [];\n"
             << "3 -- 15 [];\n"
             << "4 -- 6 [];\n"
             << "4 -- 10 [];\n"
             << "5 -- 7 [];\n"
             << "6 -- 10 [];\n"
             << "8 -- 18 [];\n"
             << "8 -- 12 [];\n"
             << "8 -- 11 [];\n"
             << "9 -- 20 [];\n"
             << "9 -- 16 [];\n"
             << "11 -- 18 [];\n"
             << "11 -- 12 [];\n"
             << "12 -- 18 [];\n"
             << "13 -- 21 [];\n"
             << "13 -- 17 [];\n"
             << "14 -- 22 [];\n"
             << "14 -- 15 [];\n"
             << "15 -- 22 [];\n"
             << "16 -- 20 [];\n"
             << "17 -- 21 [];\n"
             << "\n"
             << "4 -- 13 [len=3.0, arrowhead=vee];\n"
             << "13 -- 7 [len=3.0, arrowhead=vee];\n"
             << "7 -- 8 [len=3.0, arrowhead=vee];\n"
             << "8 -- 15 [len=3.0, arrowhead=vee];\n"
             << "15 -- 9 [len=3.0, arrowhead=vee];\n"
             << "10 -- 0 [len=3.0, arrowhead=vee];\n"
             << "0 -- 11 [len=3.0, arrowhead=vee];\n"
             << "11 -- 1 [len=3.0, arrowhead=vee];\n"
             << "1 -- 14 [len=3.0, arrowhead=vee];\n"
             << "14 -- 16 [len=3.0, arrowhead=vee];\n"
             << "6 -- 17 [len=3.0, arrowhead=vee];\n"
             << "17 -- 5 [len=3.0, arrowhead=vee];\n"
             << "5 -- 18 [len=3.0, arrowhead=vee];\n"
             << "18 -- 19 [len=3.0, arrowhead=vee];\n"
             << "19 -- 3 [len=3.0, arrowhead=vee];\n"
             << "3 -- 20 [len=3.0, arrowhead=vee];\n"
             << "21 -- 12 [len=3.0, arrowhead=vee];\n"
             << "12 -- 22 [len=3.0, arrowhead=vee];\n"
             << "22 -- 2 [len=3.0, arrowhead=vee];\n"
             << "\n"
             << "}\n";
    SEQAN_ASSERT_EQ(sstream.str(), expected.str());
}



//////////////////////////////////////////////////////////////////////////////

SEQAN_DEFINE_TEST(Test_Refinement_HeaviestCommonSubsequence)
{
    typedef String<AminoAcid> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, int> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TString s1 = "aaa";
    TString s2 = "aa";
    TStringSet strSet;
    assignValueById(strSet, s1);
    assignValueById(strSet, s2);
    TGraph g(strSet);
    TVertexDescriptor v0 = addVertex(g, 0, 0, 1);
    TVertexDescriptor v1 = addVertex(g, 0, 1, 1);
    TVertexDescriptor v2 = addVertex(g, 0, 2, 1);
    TVertexDescriptor v3 = addVertex(g, 1, 0, 1);
    TVertexDescriptor v4 = addVertex(g, 1, 1, 1);
    addEdge(g, 0, 3, 10); addEdge(g, 0, 4, 15);
    addEdge(g, 1, 3, 10); addEdge(g, 1, 4, 10);
    addEdge(g, 2, 3, 15); addEdge(g, 2, 4, 10);
    String<String<TVertexDescriptor> > str1;
    String<String<TVertexDescriptor> > str2;
    String<String<TVertexDescriptor> > align;
    String<TVertexDescriptor> tmp;
    clear(tmp); appendValue(tmp, v0); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v1); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v2); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v3); appendValue(str2, tmp);
    clear(tmp); appendValue(tmp, v4); appendValue(str2, tmp);
    SEQAN_ASSERT_EQ(heaviestCommonSubsequence(g, str1, str2, align), 20);
    // TODO(bkehr): Add more assertions.

    s1 = "aaaaa";
    s2 = "aaa";
    clear(strSet);
    assignValueById(strSet, s1);
    assignValueById(strSet, s2);
    assignStringSet(g, strSet);
    v0 = addVertex(g, 0, 0, 2);
    v1 = addVertex(g, 0, 2, 1);
    v2 = addVertex(g, 0, 3, 2);
    v3 = addVertex(g, 1, 0, 1);
    v4 = addVertex(g, 1, 1, 2);
    addEdge(g, 0, 4, 10); addEdge(g, 1, 3, 20);
    clear(align);
    clear(str1);
    clear(str2);
    clear(tmp); appendValue(tmp, v0); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v1); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v2); appendValue(str1, tmp);
    clear(tmp); appendValue(tmp, v3); appendValue(str2, tmp);
    clear(tmp); appendValue(tmp, v4); appendValue(str2, tmp);
    heaviestCommonSubsequence(g, str1, str2, align);
    // TODO(bkehr): Add assertions.
}




//////////////////////////////////////////////////////////////////////////////

// Graph AlignmentOutEdgeIterator
SEQAN_DEFINE_TEST(Test_Refinement_OutEdgeIteratorAlignment)
{
    typedef String<Dna> TString;
    typedef StringSet<TString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    //typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

    TString str1 = "aa";
    TString str2 = "ac";
    TStringSet strSet;
    appendValue(strSet, str1);
    appendValue(strSet, str2);
    TGraph g(strSet);
    TVertexDescriptor v0 = addVertex(g,0,0,1);
    addVertex(g,0,1,1);
    addVertex(g,1,0,1);
    addEdge(g,0,2,10);
    addVertex(g,1,1,1);

    typedef Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    TOutEdgeIterator it(g, v0);
    // Slow
    SEQAN_ASSERT_EQ(sourceVertex(g, getValue(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, getValue(it)), 2u);
    SEQAN_ASSERT_EQ(sourceVertex(g, value(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, *it), 2u);
    SEQAN_ASSERT_NOT(atEnd(it));
    SEQAN_ASSERT(atBegin(it));
    // Fast
    SEQAN_ASSERT_EQ(sourceVertex(it), 0u);
    SEQAN_ASSERT_EQ(targetVertex(it), 2u);
    SEQAN_ASSERT_NOT(atEnd(it));
    SEQAN_ASSERT(atBegin(it));
    ++it;
    SEQAN_ASSERT(atEnd(it));
    SEQAN_ASSERT_NOT(atBegin(it));
    goPrevious(it);
    SEQAN_ASSERT_EQ(sourceVertex(g, getValue(it)), 0u);
    SEQAN_ASSERT_EQ(targetVertex(g, getValue(it)), 2u);
    goNext(it);
    --it;

    TOutEdgeIterator it2(g, v0);
    TOutEdgeIterator it3;
    it3 = it;
    SEQAN_ASSERT(it == it2);
    SEQAN_ASSERT(it2 == it3);
    goEnd(it);
    SEQAN_ASSERT(it2 != it);
    goEnd(it2);
    SEQAN_ASSERT(it2 == it);
    goBegin(it2);
    SEQAN_ASSERT(it2 != it);
    SEQAN_ASSERT_EQ(&g,  &hostGraph(it));
}

}  // namespace seqan

#endif

