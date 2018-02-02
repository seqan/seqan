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

#include <seqan/basic.h>
#include <seqan/graph_types.h>

using namespace seqan;

SEQAN_DEFINE_TEST(test_graph_types_types_directed)
{
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

    typedef Graph<Directed<> > StandardGraph;
    typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;

    StandardGraph g;
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(empty(g) == true);

    // Add vertex
    TVertexDescriptor v0 = addVertex(g);
    SEQAN_ASSERT(v0 == 0);
    SEQAN_ASSERT(outDegree(g, v0) == 0);
    SEQAN_ASSERT(inDegree(g, 0) == 0);
    SEQAN_ASSERT(degree(g, 0) == 0);
    SEQAN_ASSERT(numVertices(g) == 1);
    SEQAN_ASSERT(empty(g) == false);

    // Add edge
    TEdgeDescriptor e1 = addEdge(g, v0, v0);
    SEQAN_ASSERT(findEdge(g, v0, v0) == e1);
    SEQAN_ASSERT(_getVertexString(g)[0] == e1);
    SEQAN_ASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1);
    SEQAN_ASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1);
    SEQAN_ASSERT(targetVertex(g, e1) == 0);
    SEQAN_ASSERT(sourceVertex(g, e1) == 0);  //Expensive in standard graph!
    SEQAN_ASSERT(numEdges(g) == 1);
    SEQAN_ASSERT(outDegree(g, v0) == 1);
    SEQAN_ASSERT(inDegree(g, v0) == 1);
    SEQAN_ASSERT(degree(g, v0) == 2);

    // Add further edges and vertices
    TVertexDescriptor v1 = addVertex(g);
    TEdgeDescriptor e2 = addEdge(g, 0, 1);
    SEQAN_ASSERT(v1 == 1);
    SEQAN_ASSERT(numVertices(g) == 2);
    SEQAN_ASSERT(targetVertex(g, e2) == 1);
    SEQAN_ASSERT(sourceVertex(g, e2) == 0);
    SEQAN_ASSERT(numEdges(g) == 2);
    SEQAN_ASSERT(outDegree(g, v0) == 2);
    SEQAN_ASSERT(inDegree(g, 1) == 1);
    SEQAN_ASSERT(inDegree(g, 0) == 1);
    SEQAN_ASSERT(degree(g, 0) == 3);

    // Add more vertices and edges
    addVertex(g);  //2
    TVertexDescriptor v3 = addVertex(g);  //3
    addVertex(g);  //4
    addEdge(g, 3, 4);
    TEdgeDescriptor my_edge = addEdge(g, 3, 1);
    addEdge(g, 3, 0);
    SEQAN_ASSERT(v3 == 3);
    SEQAN_ASSERT(numVertices(g) == 5);
    SEQAN_ASSERT(targetVertex(g, e2) == 1);
    SEQAN_ASSERT(sourceVertex(g, e2) == 0);
    SEQAN_ASSERT(targetVertex(g, my_edge) == 1);
    SEQAN_ASSERT(sourceVertex(g, my_edge) == 3);
    SEQAN_ASSERT(numEdges(g) == 5);
    SEQAN_ASSERT(outDegree(g, v3) == 3);

    // Output
    std::stringstream sstream;
    sstream << g;
    char const * EXPECTED =
        "Adjacency list:\n"
        "0 -> 1,0,\n"
        "1 -> \n"
        "2 -> \n"
        "3 -> 0,1,4,\n"
        "4 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 1 (Id: 1)\n"
        "Source: 0,Target: 0 (Id: 0)\n"
        "Source: 3,Target: 0 (Id: 4)\n"
        "Source: 3,Target: 1 (Id: 3)\n"
        "Source: 3,Target: 4 (Id: 2)\n";
    SEQAN_ASSERT(EXPECTED == sstream.str());

    // Remove edges
    removeEdge(g, my_edge);
    removeEdge(g, 0, 1);
    SEQAN_ASSERT(numEdges(g) == 3);

    // Remove vertices
    TEdgeDescriptor e3 = addEdge(g, 3, 3);
    addEdge(g, 1, 3);
    addEdge(g, 0, 3);
    addEdge(g, 0, 4);
    SEQAN_ASSERT(outDegree(g, 0) == 3);
    SEQAN_ASSERT(outDegree(g, 1) == 1);
    SEQAN_ASSERT(targetVertex(g, e3) == 3);
    SEQAN_ASSERT(sourceVertex(g, e3) == 3);
    removeVertex(g, v3);
    SEQAN_ASSERT(outDegree(g, 0) == 2);
    SEQAN_ASSERT(outDegree(g, 1) == 0);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 2);

    // Clear graph
    clearEdges(g);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 0);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clearVertices(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clear(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    addEdge(g, 4, 2);
    removeVertex(g, 3);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 2);
    SEQAN_ASSERT(inDegree(g, 4) == 0);

    // Transpose
    transpose(g);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 0);
    SEQAN_ASSERT(inDegree(g, 4) == 2);
    StandardGraph g_copy(g);
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    addVertex(g_copy);
    addEdge(g_copy, 3, 0);
    g_copy = g;
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    //Copies the graph and transposes just the copy
    transpose(g, g_copy);  // g does not change!
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 0);
    removeVertex(g, 0);



    // Adjacency matrix
    String<unsigned int> mat;
    getAdjacencyMatrix(g, mat);
    unsigned int len = (unsigned int) std::sqrt((double) length(mat));
    SEQAN_ASSERT(getValue(mat, 1 * len + 4) == 1);
    SEQAN_ASSERT(getValue(mat, 2 * len + 4) == 1);
    SEQAN_ASSERT(getValue(mat, 2 * len + 2) == 0);

    // Vertex Adjacency vectors
    String<unsigned int> vectIn, vectOut;
    getVertexAdjacencyVector(vectIn, vectOut, g, 1u);
    SEQAN_ASSERT(length(vectIn) == 0);
    SEQAN_ASSERT(vectOut[0] == 4);
    SEQAN_ASSERT(length(vectOut) == 1);
    getVertexAdjacencyVector(vectIn, vectOut, g, 4u);
    SEQAN_ASSERT(length(vectIn) == 2);
    SEQAN_ASSERT(vectIn[0] == 1);
    SEQAN_ASSERT(vectIn[1] == 2);
    SEQAN_ASSERT(length(vectOut) == 0);

//____________________________________________________________________________
//Graph with edge cargo and edge ids
    typedef Pair<char, int> TPair;
    typedef Directed<TPair> TEdges;
    typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
    typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

    Graph<TEdges> g2;
    SEQAN_ASSERT(numVertices(g2) == 0);
    SEQAN_ASSERT(numEdges(g2) == 0);
    TVertexDescriptor2 ver0 = addVertex(g2);
    SEQAN_ASSERT(ver0 == 0);
    SEQAN_ASSERT(numVertices(g2) == 1);
    TVertexDescriptor2 ver1 = addVertex(g2);
    SEQAN_ASSERT(ver1 == 1);
    SEQAN_ASSERT(numVertices(g2) == 2);
    TEdgeDescriptor2 ed1 = addEdge(g2, ver0, ver0, TPair('a', 3));
    TEdgeDescriptor2 ed2 = addEdge(g2, 0, 1);
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'a');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 3);
    SEQAN_ASSERT(targetVertex(g2, ed1) == v0);
    SEQAN_ASSERT(targetVertex(g2, ed1) == 0);
    SEQAN_ASSERT(sourceVertex(g2, ed1) == 0);
    SEQAN_ASSERT(targetVertex(g2, ed2) == 1);
    SEQAN_ASSERT(numEdges(g2) == 2);
    assignCargo(ed2, TPair('b', 4));
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'a');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 3);
    SEQAN_ASSERT((getCargo(ed2)).i1 == 'b');
    SEQAN_ASSERT((getCargo(ed2)).i2 == 4);
    cargo(ed1) = TPair('c', 1);
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'c');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 1);
    addVertex(g2);
    addVertex(g2);
    addVertex(g2);
    TEdgeDescriptor2 ed4 = addEdge(g2, 1, 4);
    cargo(ed4) = TPair('z', 100);
    removeVertex(g2, 2);
    Graph<TEdges> g2_copy(g2);
    SEQAN_ASSERT(numVertices(g2_copy) == 4);
    SEQAN_ASSERT(numEdges(g2_copy) == 3);
    clearEdges(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 4);
    SEQAN_ASSERT(numEdges(g2_copy) == 0);
    clearVertices(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 0);
    addVertex(g2_copy); addVertex(g2_copy);
    addEdge(g2_copy, 0, 1);
    clear(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 0);
    addVertex(g2_copy); addVertex(g2_copy);
    addEdge(g2_copy, 0, 1);
    SEQAN_ASSERT(numEdges(g2) == 3);
    SEQAN_ASSERT(outDegree(g2, 0) == 2);
    SEQAN_ASSERT(inDegree(g2, 0) == 1);
    transpose(g2, g2_copy);
    SEQAN_ASSERT(outDegree(g2_copy, 0) == 1);
    SEQAN_ASSERT(inDegree(g2_copy, 0) == 2);
    SEQAN_ASSERT(numEdges(g2_copy) == 3);
    TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 0, TPair('m', 3));
    SEQAN_ASSERT((getCargo(edgCargo)).i1 == 'm');
    SEQAN_ASSERT((getCargo(edgCargo)).i2 == 3);

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
    typedef Directed<void, WithoutEdgeId> TEdges3;
    typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

    Graph<TEdges3> g3;
    addVertex(g3); addVertex(g3); addVertex(g3);
    addVertex(g3); addVertex(g3);
    addEdge(g3, 1, 4);
    SEQAN_ASSERT(numVertices(g3) == 5);
    SEQAN_ASSERT(numEdges(g3) == 1);
    TEdgeDescriptor3 edge3 = addEdge(g3, 0, 4);
    //SEQAN_ASSERT(_getId(edge3) == 0);
    SEQAN_ASSERT(getCargo(edge3) == (void *) 0);
    addEdge(g3, 0, 2);
    addEdge(g3, 0, 0);
    removeEdge(g3, 0, 4);
    removeEdge(g3, 0, 2);
    SEQAN_ASSERT(numEdges(g3) == 2);

    // Multigraph
    StandardGraph multiG;
    addVertex(multiG); addVertex(multiG); addVertex(multiG);
    TEdgeDescriptor edgeD1 = addEdge(multiG, 1, 2);
    TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2);
    TEdgeDescriptor edgeD3 = addEdge(multiG, 1, 2);
    removeEdge(multiG, edgeD2);
    SEQAN_ASSERT(sourceVertex(multiG, edgeD1) == 1);
    SEQAN_ASSERT(sourceVertex(multiG, edgeD2) == 0); // EdgeDescriptor invalid
    SEQAN_ASSERT(sourceVertex(multiG, edgeD3) == 1);

//____________________________________________________________________________
//Graph with source id

    typedef Graph<Directed<void, WithSourceId> > GraphWithSourceId;
    typedef VertexDescriptor<GraphWithSourceId>::Type TVertexDescriptorSourceId;
    typedef EdgeDescriptor<GraphWithSourceId>::Type TEdgeDescriptorSourceId;

    GraphWithSourceId g4;
    SEQAN_ASSERT(numVertices(g4) == 0);
    SEQAN_ASSERT(numEdges(g4) == 0);
    SEQAN_ASSERT(empty(g4) == true);

    // Add vertex
    TVertexDescriptorSourceId v0_4 = addVertex(g4);
    SEQAN_ASSERT(v0_4 == 0);
    SEQAN_ASSERT(outDegree(g4, v0_4) == 0);
    SEQAN_ASSERT(inDegree(g4, 0) == 0);
    SEQAN_ASSERT(degree(g4, 0) == 0);
    SEQAN_ASSERT(numVertices(g4) == 1);
    SEQAN_ASSERT(empty(g4) == false);

    // Add edge
    TEdgeDescriptorSourceId e1_4 = addEdge(g4, v0_4, v0_4);
    SEQAN_ASSERT(findEdge(g4, v0_4, v0_4) == e1_4);
    SEQAN_ASSERT(_getVertexString(g4)[0] == e1_4);
    SEQAN_ASSERT(getIdUpperBound(_getVertexIdManager(g4)) == 1);
    SEQAN_ASSERT(getIdUpperBound(_getEdgeIdManager(g4)) == 1);
    SEQAN_ASSERT(targetVertex(g4, e1_4) == 0);
    SEQAN_ASSERT(sourceVertex(g4, e1_4) == 0);  // Not expensive in GraphWithSourceId!
    SEQAN_ASSERT(numEdges(g4) == 1);
    SEQAN_ASSERT(outDegree(g4, v0_4) == 1);
    SEQAN_ASSERT(inDegree(g4, v0_4) == 1);
    SEQAN_ASSERT(degree(g4, v0_4) == 2);

    // Add further edges and vertices
    TVertexDescriptorSourceId v1_4 = addVertex(g4);
    TEdgeDescriptorSourceId e2_4 = addEdge(g4, 0, 1);
    SEQAN_ASSERT(v1_4 == 1);
    SEQAN_ASSERT(numVertices(g4) == 2);
    SEQAN_ASSERT(targetVertex(g4, e2_4) == 1);
    SEQAN_ASSERT(sourceVertex(g4, e2_4) == 0);
    SEQAN_ASSERT(numEdges(g4) == 2);
    SEQAN_ASSERT(outDegree(g4, v0_4) == 2);
    SEQAN_ASSERT(inDegree(g4, 1) == 1);
    SEQAN_ASSERT(inDegree(g4, 0) == 1);
    SEQAN_ASSERT(degree(g4, 0) == 3);

    // Add more vertices and edges
    addVertex(g4);  //2
    TVertexDescriptorSourceId v3_4 = addVertex(g4);  //3
    addVertex(g4);  //4
    addEdge(g4, 3, 4);
    TEdgeDescriptorSourceId my_edge_4 = addEdge(g4, 3, 1);
    addEdge(g4, 3, 0);
    SEQAN_ASSERT(v3_4 == 3);
    SEQAN_ASSERT(numVertices(g4) == 5);
    SEQAN_ASSERT(targetVertex(g4, e2_4) == 1);
    SEQAN_ASSERT(sourceVertex(g4, e2_4) == 0);
    SEQAN_ASSERT(targetVertex(g4, my_edge_4) == 1);
    SEQAN_ASSERT(sourceVertex(g4, my_edge_4) == 3);
    SEQAN_ASSERT(numEdges(g4) == 5);
    SEQAN_ASSERT(outDegree(g4, v3_4) == 3);

    // Output
    std::stringstream sstream4;
    sstream4 << g4;
    char const * EXPECTED_4 =
        "Adjacency list:\n"
        "0 -> 1,0,\n"
        "1 -> \n"
        "2 -> \n"
        "3 -> 0,1,4,\n"
        "4 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 1 (Id: 1)\n"
        "Source: 0,Target: 0 (Id: 0)\n"
        "Source: 3,Target: 0 (Id: 4)\n"
        "Source: 3,Target: 1 (Id: 3)\n"
        "Source: 3,Target: 4 (Id: 2)\n";
    SEQAN_ASSERT(EXPECTED_4 == sstream4.str());

    // Remove edges
    SEQAN_ASSERT(outDegree(g4, v3_4) == 3);
    removeEdge(g4, my_edge_4);
    SEQAN_ASSERT(outDegree(g4, v3_4) == 2);
    removeEdge(g4, 0, 1);
    SEQAN_ASSERT(numEdges(g4) == 3);

    // Remove vertices
    TEdgeDescriptorSourceId e3_4 = addEdge(g4, 3, 3);
    addEdge(g4, 1, 3);
    addEdge(g4, 0, 3);
    addEdge(g4, 0, 4);
    SEQAN_ASSERT(outDegree(g4, 0) == 3);
    SEQAN_ASSERT(outDegree(g4, 1) == 1);
    SEQAN_ASSERT(targetVertex(g4, e3_4) == 3);
    SEQAN_ASSERT(sourceVertex(g4, e3_4) == 3);
    removeVertex(g4, v3_4);
    SEQAN_ASSERT(outDegree(g4, 0) == 2);
    SEQAN_ASSERT(outDegree(g4, 1) == 0);
    SEQAN_ASSERT(numVertices(g4) == 4);
    SEQAN_ASSERT(numEdges(g4) == 2);

    // Clear graph
    clearEdges(g4);
    SEQAN_ASSERT(numVertices(g4) == 4);
    SEQAN_ASSERT(numEdges(g4) == 0);
    addEdge(g4, 2, 0);
    addEdge(g4, 4, 1);
    clearVertices(g4);
    SEQAN_ASSERT(numVertices(g4) == 0);
    SEQAN_ASSERT(numEdges(g4) == 0);
    addVertex(g4); addVertex(g4); addVertex(g4);
    addVertex(g4); addVertex(g4);
    addEdge(g4, 2, 0);
    addEdge(g4, 4, 1);
    clear(g4);
    SEQAN_ASSERT(numVertices(g4) == 0);
    SEQAN_ASSERT(numEdges(g4) == 0);
    addVertex(g4); addVertex(g4); addVertex(g4);
    addVertex(g4); addVertex(g4);
    addEdge(g4, 2, 0);
    addEdge(g4, 4, 1);
    addEdge(g4, 4, 2);
    removeVertex(g4, 3);
    SEQAN_ASSERT(numVertices(g4) == 4);
    SEQAN_ASSERT(numEdges(g4) == 3);
    SEQAN_ASSERT(outDegree(g4, 4) == 2);
    SEQAN_ASSERT(inDegree(g4, 4) == 0);

    // Transpose
    transpose(g4);
    SEQAN_ASSERT(numVertices(g4) == 4);
    SEQAN_ASSERT(numEdges(g4) == 3);
    SEQAN_ASSERT(outDegree(g4, 4) == 0);
    SEQAN_ASSERT(inDegree(g4, 4) == 2);
    GraphWithSourceId g4_copy(g4);
    SEQAN_ASSERT(numVertices(g4_copy) == 4);
    SEQAN_ASSERT(numEdges(g4_copy) == 3);
    SEQAN_ASSERT(outDegree(g4_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g4_copy, 4) == 2);
    addVertex(g4_copy);
    addEdge(g4_copy, 3, 0);
    g4_copy = g4;
    SEQAN_ASSERT(numVertices(g4_copy) == 4);
    SEQAN_ASSERT(numEdges(g4_copy) == 3);
    SEQAN_ASSERT(outDegree(g4_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g4_copy, 4) == 2);
    //Copies the graph and transposes just the copy
    transpose(g4, g4_copy);  // g4 does not change!
    SEQAN_ASSERT(numVertices(g4_copy) == 4);
    SEQAN_ASSERT(numEdges(g4_copy) == 3);
    SEQAN_ASSERT(outDegree(g4_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g4_copy, 4) == 0);
    removeVertex(g4, 0);

    // Adjacency matrix
    String<unsigned int> mat4;
    getAdjacencyMatrix(g4, mat4);
    unsigned int len4 = (unsigned int) std::sqrt((double) length(mat4));
    SEQAN_ASSERT(getValue(mat4, 1 * len4 + 4) == 1);
    SEQAN_ASSERT(getValue(mat4, 2 * len4 + 4) == 1);
    SEQAN_ASSERT(getValue(mat4, 2 * len4 + 2) == 0);

    // Vertex Adjacency vectors
    String<unsigned int> vectIn4, vectOut4;
    getVertexAdjacencyVector(vectIn4, vectOut4, g4, 1u);
    SEQAN_ASSERT(length(vectIn4) == 0);
    SEQAN_ASSERT(vectOut4[0] == 4);
    SEQAN_ASSERT(length(vectOut4) == 1);
    getVertexAdjacencyVector(vectIn4, vectOut4, g4, 4u);
    SEQAN_ASSERT(length(vectIn4) == 2);
    SEQAN_ASSERT(vectIn4[0] == 1);
    SEQAN_ASSERT(vectIn4[1] == 2);
    SEQAN_ASSERT(length(vectOut4) == 0);
}

SEQAN_DEFINE_TEST(test_graph_types_types_undirected)
{
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

    typedef Graph<Undirected<void> > StandardGraph;
    typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;

    StandardGraph g;
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(empty(g) == true);

    // Add vertex
    TVertexDescriptor v0 = addVertex(g);
    SEQAN_ASSERT(v0 == 0);
    SEQAN_ASSERT(outDegree(g, v0) == 0);
    SEQAN_ASSERT(inDegree(g, 0) == 0);
    SEQAN_ASSERT(degree(g, 0) == 0);
    SEQAN_ASSERT(numVertices(g) == 1);
    SEQAN_ASSERT(empty(g) == false);

    // Add edge
    // TEdgeDescriptor e1 =addEdge(g,v0,v0);  // Self edges are not allowed in undirected graphs
    TVertexDescriptor v1 = addVertex(g);
    TEdgeDescriptor e = addEdge(g, 0, 1);
    SEQAN_ASSERT(findEdge(g, 0, 1) == e);
    SEQAN_ASSERT(_getVertexString(g)[0] == e);
    SEQAN_ASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2);
    SEQAN_ASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1);
    SEQAN_ASSERT(v1 == 1);
    SEQAN_ASSERT(numVertices(g) == 2);
    SEQAN_ASSERT(targetVertex(g, e) == 1);
    SEQAN_ASSERT(sourceVertex(g, e) == 0);
    SEQAN_ASSERT(numEdges(g) == 1);
    SEQAN_ASSERT(outDegree(g, v0) == 1);
    SEQAN_ASSERT(inDegree(g, 1) == 1);
    SEQAN_ASSERT(inDegree(g, 0) == 1);
    SEQAN_ASSERT(degree(g, 0) == 1);

    // Add more vertices and edges
    addVertex(g);  //2
    TVertexDescriptor v3 = addVertex(g);  //3
    addVertex(g);  //4
    addEdge(g, 3, 4);
    TEdgeDescriptor my_edge = addEdge(g, 3, 1);
    addEdge(g, 3, 0);
    SEQAN_ASSERT(v3 == 3);
    SEQAN_ASSERT(numVertices(g) == 5);
    SEQAN_ASSERT(targetVertex(g, my_edge) == 3);
    SEQAN_ASSERT(sourceVertex(g, my_edge) == 1);
    SEQAN_ASSERT(numEdges(g) == 4);
    SEQAN_ASSERT(outDegree(g, v3) == 3);
    SEQAN_ASSERT(inDegree(g, v3) == 3);
    SEQAN_ASSERT(degree(g, v3) == 3);

    // Output
    std::stringstream sstream;
    sstream << g;
    char const * EXPECTED =
        "Adjacency list:\n"
        "0 -> 3,1,\n"
        "1 -> 3,0,\n"
        "2 -> \n"
        "3 -> 0,1,4,\n"
        "4 -> 3,\n"
        "Edge list:\n"
        "Source: 0,Target: 3 (Id: 3)\n"
        "Source: 0,Target: 1 (Id: 0)\n"
        "Source: 1,Target: 3 (Id: 2)\n"
        "Source: 3,Target: 4 (Id: 1)\n";
    SEQAN_ASSERT(EXPECTED == sstream.str());

    // Remove edges
    removeEdge(g, my_edge);
    removeEdge(g, 0, 1);
    SEQAN_ASSERT(numEdges(g) == 2);


    // Remove vertices
    addVertex(g);  //5
    addEdge(g, 5, 2);
    addEdge(g, 2, 3);
    addEdge(g, 1, 3);
    addEdge(g, 1, 4);
    SEQAN_ASSERT(outDegree(g, 3) == 4);
    SEQAN_ASSERT(outDegree(g, 4) == 2);
    removeVertex(g, v3);
    SEQAN_ASSERT(outDegree(g, 4) == 1);
    SEQAN_ASSERT(outDegree(g, 0) == 0);
    SEQAN_ASSERT(numVertices(g) == 5);
    SEQAN_ASSERT(numEdges(g) == 2);

    // Clear graph
    clearEdges(g);
    SEQAN_ASSERT(numVertices(g) == 5);
    SEQAN_ASSERT(numEdges(g) == 0);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clearVertices(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    clear(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0);
    addEdge(g, 4, 1);
    addEdge(g, 4, 2);
    removeVertex(g, 3);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 2);
    SEQAN_ASSERT(inDegree(g, 4) == 2);

    // Transpose
    transpose(g);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 2);
    SEQAN_ASSERT(inDegree(g, 4) == 2);
    StandardGraph g_copy(g);
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    addVertex(g_copy);
    addEdge(g_copy, 3, 0);
    g_copy = g;
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    //Copies the graph and transposes just the copy
    transpose(g, g_copy);  // g does not change!
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);

    // Adjacency matrix
    String<unsigned int> mat;
    getAdjacencyMatrix(g, mat);
    unsigned int len = (unsigned int) std::sqrt((double) length(mat));
    SEQAN_ASSERT(getValue(mat, 0 * len + 2) == 1);
    SEQAN_ASSERT(getValue(mat, 3 * len + 2) == 0);
    SEQAN_ASSERT(getValue(mat, 0 * len + 2) == getValue(mat, 2 * len + 0));
    SEQAN_ASSERT(getValue(mat, 1 * len + 4) == getValue(mat, 4 * len + 1));
    SEQAN_ASSERT(getValue(mat, 2 * len + 4) == getValue(mat, 4 * len + 2));

    // Vertex Adjacency vectors
    String<unsigned int> vectIn, vectOut;
    getVertexAdjacencyVector(vectIn, vectOut, g, 2u);
    SEQAN_ASSERT(length(vectIn) == 2);
    SEQAN_ASSERT(vectIn[0] == 4);
    SEQAN_ASSERT(vectIn[1] == 0);
    SEQAN_ASSERT(vectOut[0] == 4);
    SEQAN_ASSERT(vectOut[1] == 0);
    SEQAN_ASSERT(length(vectOut) == 2);

//____________________________________________________________________________
//Graph with edge cargo and edge ids
    typedef Pair<char, int> TPair;
    typedef Undirected<TPair> TEdges;
    typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
    typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

    Graph<TEdges> g2;
    SEQAN_ASSERT(numVertices(g2) == 0);
    SEQAN_ASSERT(numEdges(g2) == 0);
    TVertexDescriptor2 ver0 = addVertex(g2);
    SEQAN_ASSERT(ver0 == 0);
    SEQAN_ASSERT(numVertices(g2) == 1);
    TVertexDescriptor2 ver1 = addVertex(g2);
    SEQAN_ASSERT(ver1 == 1);
    SEQAN_ASSERT(numVertices(g2) == 2);
    TEdgeDescriptor2 ed1 = addEdge(g2, 0, 1);
    SEQAN_ASSERT(targetVertex(g2, ed1) == 1);
    SEQAN_ASSERT(sourceVertex(g2, ed1) == 0);
    SEQAN_ASSERT(numEdges(g2) == 1);
    assignCargo(ed1, TPair('a', 3));
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'a');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 3);
    cargo(ed1) = TPair('c', 1);
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'c');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 1);
    addVertex(g2);
    addVertex(g2);
    addVertex(g2);
    TEdgeDescriptor2 ed4 = addEdge(g2, 1, 4);
    cargo(ed4) = TPair('z', 100);
    removeVertex(g2, 2);
    Graph<TEdges> g2_copy(g2);
    SEQAN_ASSERT(numVertices(g2_copy) == 4);
    SEQAN_ASSERT(numEdges(g2_copy) == 2);
    clearEdges(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 4);
    SEQAN_ASSERT(numEdges(g2_copy) == 0);
    clearVertices(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 0);
    addVertex(g2_copy); addVertex(g2_copy);
    addEdge(g2_copy, 0, 1);
    clear(g2_copy);
    SEQAN_ASSERT(numVertices(g2_copy) == 0);
    addVertex(g2_copy); addVertex(g2_copy);
    addEdge(g2_copy, 0, 1);
    transpose(g2, g2_copy);
    SEQAN_ASSERT(outDegree(g2_copy, 0) == 1);
    SEQAN_ASSERT(inDegree(g2_copy, 0) == 1);
    SEQAN_ASSERT(numEdges(g2_copy) == 2);
    TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 3, TPair('m', 3));
    SEQAN_ASSERT((getCargo(edgCargo)).i1 == 'm');
    SEQAN_ASSERT((getCargo(edgCargo)).i2 == 3);

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
    typedef Undirected<void, WithoutEdgeId> TEdges3;
    typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

    Graph<TEdges3> g3;
    addVertex(g3); addVertex(g3); addVertex(g3);
    addVertex(g3); addVertex(g3);
    addEdge(g3, 1, 4);
    SEQAN_ASSERT(numVertices(g3) == 5);
    SEQAN_ASSERT(numEdges(g3) == 1);
    TEdgeDescriptor3 edge3 = addEdge(g3, 0, 4);
    SEQAN_ASSERT(_getId(edge3) == 0);
    SEQAN_ASSERT(getCargo(edge3) == (void *) 0);
    addEdge(g3, 0, 2);
    addEdge(g3, 0, 1);
    removeEdge(g3, 0, 4);
    removeEdge(g3, 0, 2);
    SEQAN_ASSERT(numEdges(g3) == 2);
    removeInEdges(g3, 1);
    SEQAN_ASSERT(numEdges(g3) == 0);


//____________________________________________________________________________
// Undirected graph iterators
    typedef Graph<Undirected<> > TGraphIter;

    TGraphIter gIter;
    addVertex(gIter); addVertex(gIter); addVertex(gIter); addVertex(gIter);
    addVertex(gIter); addVertex(gIter); addVertex(gIter); addVertex(gIter);
    removeVertex(gIter, 0);
    removeVertex(gIter, 5);
    addEdge(gIter, 2, 7);
    addEdge(gIter, 2, 3);
    addEdge(gIter, 2, 4);
    addEdge(gIter, 4, 3);
    addEdge(gIter, 3, 6);
    addEdge(gIter, 4, 6);

    typedef Iterator<TGraphIter, OutEdgeIterator>::Type TOutEdgeIterator;
    TOutEdgeIterator itOutEdge(gIter, 3);
    // Both ways are fast for undirected graphs
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itOutEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itOutEdge)) == 6);
    SEQAN_ASSERT(sourceVertex(itOutEdge) == 3);
    SEQAN_ASSERT(targetVertex(itOutEdge) == 6);
    SEQAN_ASSERT(sourceVertex(gIter, value(itOutEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, *itOutEdge) == 6);
    SEQAN_ASSERT(atEnd(itOutEdge) == false);
    SEQAN_ASSERT(atBegin(itOutEdge) == true);
    goNext(itOutEdge);
    SEQAN_ASSERT(atEnd(itOutEdge) == false);
    SEQAN_ASSERT(atBegin(itOutEdge) == false);
    SEQAN_ASSERT(sourceVertex(itOutEdge) == 3);
    SEQAN_ASSERT(targetVertex(itOutEdge) == 4);
    ++itOutEdge;
    itOutEdge++;
    SEQAN_ASSERT(atEnd(itOutEdge) == true);
    SEQAN_ASSERT(atBegin(itOutEdge) == false);
    goPrevious(itOutEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itOutEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itOutEdge)) == 3);
    --itOutEdge;
    itOutEdge--;
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itOutEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itOutEdge)) == 6);
    itOutEdge--;
    itOutEdge--;
    SEQAN_ASSERT(atBegin(itOutEdge) == true);
    TOutEdgeIterator itEdge2(itOutEdge);
    TOutEdgeIterator itEdge3;
    itEdge3 = itOutEdge;
    SEQAN_ASSERT(itOutEdge == itEdge2);
    SEQAN_ASSERT(itEdge2 == itEdge3);
    goEnd(itOutEdge);
    SEQAN_ASSERT(itEdge2 != itOutEdge);
    goEnd(itEdge2);
    SEQAN_ASSERT(itEdge2 == itOutEdge);
    goBegin(itEdge2);
    SEQAN_ASSERT(atBegin(itEdge2) == true);
    SEQAN_ASSERT(&gIter == &hostGraph(itOutEdge));


    typedef Iterator<TGraphIter, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEdge(gIter);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 4);
    SEQAN_ASSERT(atBegin(itEdge) == true);
    SEQAN_ASSERT(atEnd(itEdge) == false);
    goNext(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 3);
    SEQAN_ASSERT(atBegin(itEdge) == false);
    SEQAN_ASSERT(atEnd(itEdge) == false);
    goNext(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 7);
    ++itEdge;
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 6);
    itEdge++;
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 4);
    goNext(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 4);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 6);
    goNext(itEdge);
    SEQAN_ASSERT(atBegin(itEdge) == false);
    SEQAN_ASSERT(atEnd(itEdge) == true);
    goPrevious(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 4);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 6);
    --itEdge;
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 4);
    itEdge--;
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 3);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 6);
    goPrevious(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 7);
    goPrevious(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 3);
    goPrevious(itEdge);
    SEQAN_ASSERT(sourceVertex(gIter, getValue(itEdge)) == 2);
    SEQAN_ASSERT(targetVertex(gIter, getValue(itEdge)) == 4);
    SEQAN_ASSERT(atBegin(itEdge) == true);
    SEQAN_ASSERT(atEnd(itEdge) == false);

    // Multigraph
    StandardGraph multiG;
    addVertex(multiG); addVertex(multiG); addVertex(multiG);
    addEdge(multiG, 1, 2);
    TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2);
    addEdge(multiG, 1, 2);
    removeEdge(multiG, edgeD2);
    SEQAN_ASSERT(numEdges(multiG) == 2);
}

SEQAN_DEFINE_TEST(test_graph_types_types_automaton)
{
//____________________________________________________________________________
// Standard automaton: No edge cargo

    typedef Graph<Automaton<Dna> > StandardAutomaton;
    typedef VertexDescriptor<StandardAutomaton>::Type TVertexDescriptor;
    typedef EdgeDescriptor<StandardAutomaton>::Type TEdgeDescriptor;

    StandardAutomaton g;
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(empty(g) == true);

    // Add vertex
    createRoot(g);
    TVertexDescriptor v0 = getRoot(g);
    SEQAN_ASSERT(v0 == 0);
    SEQAN_ASSERT(outDegree(g, v0) == 0);
    SEQAN_ASSERT(inDegree(g, 0) == 0);
    SEQAN_ASSERT(degree(g, 0) == 0);
    SEQAN_ASSERT(numVertices(g) == 1);
    SEQAN_ASSERT(empty(g) == false);

    // Add edge
    TEdgeDescriptor e1 = addEdge(g, v0, v0, 'a');
    SEQAN_ASSERT(findEdge(g, 0, 'a') == e1);
    SEQAN_ASSERT(&_getVertexString(g)[0].data_edge[0] == e1);
    SEQAN_ASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1);
    SEQAN_ASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1);
    SEQAN_ASSERT(_getId(e1) == 0);
    SEQAN_ASSERT(_getId(e1) == 0);
    SEQAN_ASSERT(targetVertex(g, e1) == 0);
    SEQAN_ASSERT(sourceVertex(g, e1) == 0);
    SEQAN_ASSERT(numEdges(g) == 1);
    SEQAN_ASSERT(outDegree(g, v0) == 1);
    SEQAN_ASSERT(inDegree(g, v0) == 1);
    SEQAN_ASSERT(degree(g, v0) == 2);

    // Add further edges and vertices
    TVertexDescriptor v1 = addVertex(g);
    TEdgeDescriptor e2 = addEdge(g, 0, 1, 'g');
    SEQAN_ASSERT(_getId(e2) == 1);
    SEQAN_ASSERT(v1 == 1);
    SEQAN_ASSERT(numVertices(g) == 2);
    SEQAN_ASSERT(targetVertex(g, e2) == 1);
    SEQAN_ASSERT(sourceVertex(g, e2) == 0);
    SEQAN_ASSERT(numEdges(g) == 2);
    SEQAN_ASSERT(outDegree(g, v0) == 2);
    SEQAN_ASSERT(inDegree(g, 1) == 1);
    SEQAN_ASSERT(inDegree(g, 0) == 1);
    SEQAN_ASSERT(degree(g, 0) == 3);

    // Add more vertices and edges
    addVertex(g);  //2
    TVertexDescriptor v3 = addVertex(g);  //3
    addVertex(g);  //4
    addEdge(g, 3, 4, 'g');
    TEdgeDescriptor my_edge = addEdge(g, 3, 1, 'c');
    SEQAN_ASSERT(_getId(my_edge) == 3);
    addEdge(g, 3, 0, 't');
    SEQAN_ASSERT(v3 == 3);
    SEQAN_ASSERT(numVertices(g) == 5);
    SEQAN_ASSERT(targetVertex(g, e2) == 1);
    SEQAN_ASSERT(sourceVertex(g, e2) == 0);
    SEQAN_ASSERT(targetVertex(g, my_edge) == 1);
    SEQAN_ASSERT(sourceVertex(g, my_edge) == 3);
    SEQAN_ASSERT(numEdges(g) == 5);
    SEQAN_ASSERT(outDegree(g, v3) == 3);

    // Output
    std::stringstream sstream;
    sstream << g;
    char const * EXPECTED =
        "Automaton - State: (Input / NextState)\n"
        "0:  (A / 0)  (C / nil)  (G / 1)  (T / nil) \n"
        "1:  (A / nil)  (C / nil)  (G / nil)  (T / nil) \n"
        "2:  (A / nil)  (C / nil)  (G / nil)  (T / nil) \n"
        "3:  (A / nil)  (C / 1)  (G / 4)  (T / 0) \n"
        "4:  (A / nil)  (C / nil)  (G / nil)  (T / nil) \n";
    SEQAN_ASSERT(EXPECTED == sstream.str());

    // Remove edges
    removeEdge(g, 3, 1, 'c');
    removeEdge(g, 0, 1, 'g');
    SEQAN_ASSERT(numEdges(g) == 3);

    // Remove vertices
    TEdgeDescriptor e3 = addEdge(g, 3, 3, 'a');
    addEdge(g, 1, 3, 'a');
    addEdge(g, 0, 3, 'c');
    addEdge(g, 0, 4, 't');
    SEQAN_ASSERT(outDegree(g, 0) == 3);
    SEQAN_ASSERT(outDegree(g, 1) == 1);
    SEQAN_ASSERT(targetVertex(g, e3) == 3);
    SEQAN_ASSERT(sourceVertex(g, e3) == 3);
    SEQAN_ASSERT(numEdges(g) == 7);
    removeVertex(g, v3);
    SEQAN_ASSERT(outDegree(g, 0) == 2);
    SEQAN_ASSERT(outDegree(g, 1) == 0);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 2);

    // Clear graph
    clearEdges(g);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 0);
    addEdge(g, 2, 0, 'a');
    addEdge(g, 4, 1, 'c');
    clearVertices(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0, 't');
    addEdge(g, 4, 1, 'g');
    clear(g);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    addVertex(g); addVertex(g); addVertex(g);
    addVertex(g); addVertex(g);
    addEdge(g, 2, 0, 'c');
    addEdge(g, 4, 1, 'g');
    addEdge(g, 4, 2, 't');
    removeVertex(g, 3);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 2);
    SEQAN_ASSERT(inDegree(g, 4) == 0);

    //Transposes the graph in-place
    transpose(g);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(outDegree(g, 4) == 0);
    SEQAN_ASSERT(inDegree(g, 4) == 2);
    StandardAutomaton g_copy(g);
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    addVertex(g_copy);
    addEdge(g_copy, 3, 0, 'a');
    g_copy = g;
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 0);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 2);
    //Copies the graph and transposes just the copy
    transpose(g, g_copy);  // g does not change!
    SEQAN_ASSERT(numVertices(g_copy) == 4);
    SEQAN_ASSERT(numEdges(g_copy) == 3);
    SEQAN_ASSERT(outDegree(g_copy, 4) == 2);
    SEQAN_ASSERT(inDegree(g_copy, 4) == 0);
    removeVertex(g, 0);

    // Adjacency matrix
    String<unsigned int> mat;
    getAdjacencyMatrix(g, mat);
    unsigned int len = (unsigned int) std::sqrt((double) length(mat));
    SEQAN_ASSERT(getValue(mat, 1 * len + 4) == 1);
    SEQAN_ASSERT(getValue(mat, 2 * len + 4) == 1);
    SEQAN_ASSERT(getValue(mat, 0 * len + 2) == 0);

    // Vertex Adjacency vectors
    String<unsigned int> vectIn, vectOut;
    getVertexAdjacencyVector(vectIn, vectOut, g, 1u);
    SEQAN_ASSERT(length(vectIn) == 0);
    SEQAN_ASSERT(length(vectOut) == 1);
    SEQAN_ASSERT(vectOut[0] == 4);
    getVertexAdjacencyVector(vectIn, vectOut, g, 4u);
    SEQAN_ASSERT(length(vectIn) == 2);
    SEQAN_ASSERT(length(vectOut) == 0);
    SEQAN_ASSERT(vectIn[0] == 1);
    SEQAN_ASSERT(vectIn[1] == 2);

    // Test iterators
    typedef Iterator<StandardAutomaton, VertexIterator>::Type TVertexIterator;
    TVertexIterator itVert(g);
    SEQAN_ASSERT(getValue(itVert) == 1);
    ++itVert;
    SEQAN_ASSERT(getValue(itVert) == 2);
    itVert++;
    SEQAN_ASSERT(getValue(itVert) == 4);
    goNext(itVert);
    SEQAN_ASSERT(atEnd(itVert) == true);

    addEdge(g, 1, 2, 'T');
    typedef Iterator<StandardAutomaton, OutEdgeIterator>::Type TOutEdgeIterator;
    TOutEdgeIterator itEdge(g, 1);
    // Slow
    SEQAN_ASSERT(sourceVertex(g, getValue(itEdge)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEdge)) == 4);
    // Fast
    SEQAN_ASSERT(sourceVertex(itEdge) == 1);
    SEQAN_ASSERT(targetVertex(itEdge) == 4);

    SEQAN_ASSERT(sourceVertex(g, value(itEdge)) == 1);
    SEQAN_ASSERT(targetVertex(g, *itEdge) == 4);
    SEQAN_ASSERT(atEnd(itEdge) == false);
    SEQAN_ASSERT(atBegin(itEdge) == true);
    goNext(itEdge);
    SEQAN_ASSERT(sourceVertex(itEdge) == 1);
    SEQAN_ASSERT(targetVertex(itEdge) == 2);
    ++itEdge;
    itEdge++;
    SEQAN_ASSERT(atEnd(itEdge) == true);
    SEQAN_ASSERT(atBegin(itEdge) == false);
    goPrevious(itEdge);
    SEQAN_ASSERT(sourceVertex(g, getValue(itEdge)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEdge)) == 2);
    --itEdge;
    itEdge++;
    SEQAN_ASSERT(sourceVertex(g, getValue(itEdge)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEdge)) == 2);
    itEdge--;
    itEdge--;
    SEQAN_ASSERT(atBegin(itEdge) == true);
    TOutEdgeIterator itEdge2(itEdge);
    TOutEdgeIterator itEdge3;
    itEdge3 = itEdge;
    SEQAN_ASSERT(itEdge == itEdge2);
    SEQAN_ASSERT(itEdge2 == itEdge3);
    goEnd(itEdge);
    SEQAN_ASSERT(itEdge2 != itEdge);
    goEnd(itEdge2);
    SEQAN_ASSERT(itEdge2 == itEdge);
    goBegin(itEdge2);
    SEQAN_ASSERT(atBegin(itEdge2) == true);
    SEQAN_ASSERT(&g == &hostGraph(itEdge));

    typedef Iterator<StandardAutomaton, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEd(g);
    SEQAN_ASSERT(sourceVertex(g, getValue(itEd)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEd)) == 4);
    SEQAN_ASSERT(sourceVertex(g, value(itEd)) == 1);
    SEQAN_ASSERT(targetVertex(g, *itEd) == 4);
    SEQAN_ASSERT(atEnd(itEd) == false);
    SEQAN_ASSERT(atBegin(itEd) == true);
    goNext(itEd);
    SEQAN_ASSERT(sourceVertex(g, getValue(itEd)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEd)) == 2);
    SEQAN_ASSERT(atEnd(itEd) == false);
    SEQAN_ASSERT(atBegin(itEd) == false);
    ++itEd;
    SEQAN_ASSERT(atEnd(itEd) == false);
    SEQAN_ASSERT(atBegin(itEd) == false);
    // Slow
    SEQAN_ASSERT(sourceVertex(g, getValue(itEd)) == 2);
    SEQAN_ASSERT(targetVertex(g, getValue(itEd)) == 4);
    // Fast
    SEQAN_ASSERT(sourceVertex(itEd) == 2);
    SEQAN_ASSERT(targetVertex(itEd) == 4);
    itEd++;
    itEd++;
    SEQAN_ASSERT(atEnd(itEd) == true);
    SEQAN_ASSERT(atBegin(itEd) == false);
    goPrevious(itEd);
    SEQAN_ASSERT(sourceVertex(g, getValue(itEd)) == 2);
    SEQAN_ASSERT(targetVertex(g, getValue(itEd)) == 4);
    --itEd;
    SEQAN_ASSERT(sourceVertex(g, getValue(itEd)) == 1);
    SEQAN_ASSERT(targetVertex(g, getValue(itEd)) == 2);
    TEdgeIterator itEd2(g);
    TEdgeIterator itEd3;
    goBegin(itEd);
    itEd3 = itEd;
    SEQAN_ASSERT(itEd == itEd2);
    SEQAN_ASSERT(itEd2 == itEd3);
    goEnd(itEd);
    SEQAN_ASSERT(itEd2 != itEd);
    goEnd(itEd2);
    SEQAN_ASSERT(itEd2 == itEd);
    goBegin(itEd2);
    SEQAN_ASSERT(itEd2 != itEd);
    SEQAN_ASSERT(&hostGraph(itEd) == &g);

    typedef Iterator<StandardAutomaton, AdjacencyIterator>::Type TAdjacencyIterator;
    TAdjacencyIterator itAd(g, 1);
    SEQAN_ASSERT(getValue(itAd) == 4);
    SEQAN_ASSERT(&hostGraph(itAd) == &g);
    SEQAN_ASSERT(value(itAd) == 4);
    SEQAN_ASSERT(*itAd == 4);
    SEQAN_ASSERT(atEnd(itAd) == false);
    SEQAN_ASSERT(atBegin(itAd) == true);
    goNext(itAd);
    SEQAN_ASSERT(getValue(itAd) == 2);
    SEQAN_ASSERT(atEnd(itAd) == false);
    SEQAN_ASSERT(atBegin(itAd) == false);
    ++itAd;
    SEQAN_ASSERT(atEnd(itAd) == true);
    SEQAN_ASSERT(atBegin(itAd) == false);
    goBegin(itAd);
    itAd++;
    itAd++;
    itAd++;
    SEQAN_ASSERT(atEnd(itAd) == true);
    SEQAN_ASSERT(atBegin(itAd) == false);
    goPrevious(itAd);
    SEQAN_ASSERT(getValue(itAd) == 2);
    --itAd;
    SEQAN_ASSERT(getValue(itAd) == 4);
    SEQAN_ASSERT(atEnd(itAd) == false);
    SEQAN_ASSERT(atBegin(itAd) == true);
    goEnd(itAd);
    itAd--;
    SEQAN_ASSERT(getValue(itAd) == 2);
    goBegin(itAd);
    TAdjacencyIterator itAd2(itAd);
    TAdjacencyIterator itAd3;
    itAd3 = itAd;
    SEQAN_ASSERT(itAd == itAd2);
    SEQAN_ASSERT(itAd2 == itAd3);
    goEnd(itAd);
    SEQAN_ASSERT(itAd2 != itAd);
    goEnd(itAd2);
    SEQAN_ASSERT(itAd2 == itAd);
    goBegin(itAd2);
    SEQAN_ASSERT(itAd2 != itAd);



//____________________________________________________________________________
// Automaton - Different alphabet
    typedef VertexDescriptor<Graph<Automaton<char> > >::Type VertexDescriptorType;
    Graph<Automaton<char> > automaton;
    VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
    addVertex(automaton); // B = 1
    addVertex(automaton); // C = 2
    addVertex(automaton); // D = 3
    addVertex(automaton); // E = 4
    addVertex(automaton); // F = 5
    addEdge(automaton, 0, 1, '2');
    addEdge(automaton, 1, 0, '1');
    addEdge(automaton, 4, 0, '6');
    addEdge(automaton, 0, 3, '7');
    addEdge(automaton, 1, 1, '3');
    addEdge(automaton, 1, 2, '4');
    addEdge(automaton, 5, 1, '8');
    addEdge(automaton, 2, 5, '5');
    addEdge(automaton, 3, 4, '2');
    addEdge(automaton, 5, 3, '7');

    VertexDescriptorType succ;
    succ = getSuccessor(automaton, rootVertex, '7');
    SEQAN_ASSERT(succ == 3);
    // Throws an error in debug mode because edge does not exist
    //succ = getSuccessor(automaton,rootVertex,'6');
    succ = getSuccessor(automaton, succ, '2');
    SEQAN_ASSERT(succ == 4);
    succ = getSuccessor(automaton, succ, '6');
    SEQAN_ASSERT(succ == 0);
    // If no map is specified it is assumed that an edge cargo exists!!!
    succ = getSuccessor(automaton, succ, '2');
    SEQAN_ASSERT(succ == 1);

    // Now using shortcuts
    SEQAN_ASSERT(canParseString(automaton, rootVertex, "7262"));
    SEQAN_ASSERT(!canParseString(automaton, rootVertex, "726C"));
    SEQAN_ASSERT(canParseString(automaton, "7262"));
    SEQAN_ASSERT(!canParseString(automaton, "726C"));
    succ = parseString(automaton, rootVertex, "7262");
    SEQAN_ASSERT(succ == 1);
    std::string str = "7262";
    succ = parseString(automaton, rootVertex, str.begin(), str.end());
    SEQAN_ASSERT(succ == 1);
    String<char> str2("7262");
    succ = parseString(automaton, rootVertex, begin(str2), end(str2));
    SEQAN_ASSERT(succ == 1);
    String<char> input("7262");
    SEQAN_ASSERT(canParseString(automaton, rootVertex, input));
    SEQAN_ASSERT(canParseString(automaton, input));
    succ = parseString(automaton, rootVertex, input);
    SEQAN_ASSERT(succ == 1);

    // Additional cargo
    typedef Graph<Automaton<Dna, short> > TGraph9;
    typedef EdgeDescriptor<TGraph9>::Type TEdgeDescriptor9;

    TGraph9 g9;
    addVertex(g9);
    addVertex(g9);
    Dna aDna('a');
    Dna gDna('g');
    TEdgeDescriptor9 edg1 = addEdge(g9, 0, 1, aDna, 12);
    TEdgeDescriptor9 edg2 = addEdge(g9, 1, 0, gDna, 21);
    TGraph9 g10;
    transpose(g9, g10);
    TEdgeDescriptor9 edg1_10 = findEdge(g10, 1, aDna);
    TEdgeDescriptor9 edg1_11 = findEdge(g10, 0, gDna);
    SEQAN_ASSERT(getCargo(edg1) == 12);
    SEQAN_ASSERT(getCargo(edg2) == 21);
    SEQAN_ASSERT(getCargo(edg1_10) == 12);
    SEQAN_ASSERT(getCargo(edg1_11) == 21);

    // Multigraph
    StandardAutomaton multiG;
    addVertex(multiG); addVertex(multiG); addVertex(multiG);
    addEdge(multiG, 1, 2, 'a');
    TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2, 'c');
    addEdge(multiG, 1, 2, 'g');
    removeEdge(multiG, edgeD2);
    SEQAN_ASSERT(numEdges(multiG) == 2);
}

SEQAN_DEFINE_TEST(test_graph_types_types_word_graph)
{
//____________________________________________________________________________
// Standard automaton: No edge cargo

    typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
    typedef VertexDescriptor<TWordGraph>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TWordGraph>::Type TEdgeDescriptor;


    TWordGraph g;
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(empty(g) == true);

    // Add vertex
    TVertexDescriptor v0 = addVertex(g);
    SEQAN_ASSERT(v0 == 0);
    SEQAN_ASSERT(outDegree(g, v0) == 0);
    SEQAN_ASSERT(inDegree(g, 0) == 0);
    SEQAN_ASSERT(degree(g, 0) == 0);
    SEQAN_ASSERT(numVertices(g) == 1);
    SEQAN_ASSERT(empty(g) == false);
    addVertex(g);
    addVertex(g);
    TVertexDescriptor v3 = addVertex(g);
    SEQAN_ASSERT(isRoot(g, 0) == true);
    SEQAN_ASSERT(getRoot(g) == 0);
    assignRoot(g, 3);
    SEQAN_ASSERT(getRoot(g) == 3);
    SEQAN_ASSERT(isRoot(g, 0) == false);
    SEQAN_ASSERT(isRoot(g, 3) == true);
    root(g) = 2;
    SEQAN_ASSERT(getRoot(g) == 2);
    SEQAN_ASSERT(isRoot(g, 3) == false);
    SEQAN_ASSERT(isRoot(g, 2) == true);

    // Add edge
    TEdgeDescriptor e1 = addEdge(g, v0, v3, "ag");
    SEQAN_ASSERT(findEdge(g, v0, 'a') == e1);
    SEQAN_ASSERT(_getId(e1) == 0);
    // First letter -> edge label, all other letters into the cargo
    SEQAN_ASSERT(getCargo(e1) == "g");
    SEQAN_ASSERT(targetVertex(g, e1) == 3);
    SEQAN_ASSERT(sourceVertex(g, e1) == 0);
    SEQAN_ASSERT(numEdges(g) == 1);
    SEQAN_ASSERT(outDegree(g, v0) == 1);
    SEQAN_ASSERT(inDegree(g, v0) == 0);
    SEQAN_ASSERT(degree(g, v0) == 1);

    // Add further edges and vertices
    addVertex(g);
    TVertexDescriptor v5 = addVertex(g);
    TEdgeDescriptor e2 = addEdge(g, 0, 5, "g");
    SEQAN_ASSERT(_getId(e2) == 1);
    SEQAN_ASSERT(v5 == 5);
    SEQAN_ASSERT(numVertices(g) == 6);
    SEQAN_ASSERT(targetVertex(g, e2) == 5);
    SEQAN_ASSERT(sourceVertex(g, e2) == 0);
    SEQAN_ASSERT(numEdges(g) == 2);
    removeEdge(g, 0, 5, String<Dna>("g"));
    SEQAN_ASSERT(numEdges(g) == 1);
    e2 = addEdge(g, 0, 5, "g");
    SEQAN_ASSERT(outDegree(g, v0) == 2);
    SEQAN_ASSERT(inDegree(g, 5) == 1);
    SEQAN_ASSERT(degree(g, 0) == 2);
    SEQAN_ASSERT(getSuccessor(g, 0, "g") == 5);
    SEQAN_ASSERT(getSuccessor(g, 0, String<Dna>("ag")) == 3);  // The whole edge label or just the first letter
    SEQAN_ASSERT(getSuccessor(g, 0, "a") == getNil<TVertexDescriptor>());
    addVertex(g);
    addVertex(g);
    addEdge(g, 3, 1, "aggg");
    addEdge(g, 3, 4, "gg");
    addEdge(g, 5, 2, "aggg");
    addEdge(g, 5, 7, "g");
    addEdge(g, 7, 6, "g");
    SEQAN_ASSERT(parseString(g, 0, "agaggg") == 1);
    SEQAN_ASSERT(parseString(g, 0, "aga") == 3);  // Does not reach 1
    SEQAN_ASSERT(parseString(g, 0, "g") == 5);
    SEQAN_ASSERT(parseString(g, 0, "ggg") == 6);
    SEQAN_ASSERT(parseString(g, 0, "gaggg") == 2);
    SEQAN_ASSERT(parseString(g, 0, "gagggg") == 2);
    assignRoot(g, 0);

    // Output
    std::stringstream sstream;
    sstream << g;
    char const * EXPECTED =
        "WordGraph - Directed:\n"
        "0->3  Label: AG\n"
        "0->5  Label: G\n"
        "3->1  Label: AGGG\n"
        "3->4  Label: GG\n"
        "5->2  Label: AGGG\n"
        "5->7  Label: G\n"
        "7->6  Label: G\n";
    SEQAN_ASSERT_EQ(EXPECTED, sstream.str());

    assignRoot(g, 2);
    TWordGraph g_tmp(g);
    SEQAN_ASSERT(numVertices(g_tmp) == 8);
    SEQAN_ASSERT(parseString(g_tmp, 0, "agaggg") == 1);
    SEQAN_ASSERT(inDegree(g_tmp, 5) == 1);
    SEQAN_ASSERT(degree(g_tmp, 0) == 2);
    SEQAN_ASSERT(isRoot(g_tmp, 2) == true);
    TWordGraph g_assign;
    g_assign = g;
    SEQAN_ASSERT(numVertices(g_assign) == 8);
    SEQAN_ASSERT(parseString(g_assign, 0, "agaggg") == 1);
    SEQAN_ASSERT(inDegree(g_assign, 5) == 1);
    SEQAN_ASSERT(degree(g_assign, 0) == 2);

    // Transpose
    transpose(g, g_tmp);
    SEQAN_ASSERT(numVertices(g_tmp) == 8);
    SEQAN_ASSERT(parseString(g_tmp, 2, "aggg") == 5);
    SEQAN_ASSERT(inDegree(g_tmp, 5) == 2);
    SEQAN_ASSERT(outDegree(g_tmp, 0) == 0);
    SEQAN_ASSERT(isRoot(g_tmp, 2) == true);
}

SEQAN_DEFINE_TEST(test_graph_types_types_word_tree)
{
//____________________________________________________________________________
// Tree without edge cargo

    typedef Graph<Tree<void> > TTree;
    typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
    typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;

    TTree g;
    SEQAN_ASSERT(empty(g) == true);
    createRoot(g);
    TVertexDescriptor rootV = getRoot(g);
    SEQAN_ASSERT(rootV == 0);
    SEQAN_ASSERT(isRoot(g, rootV) == true);
    SEQAN_ASSERT(root(g) == rootV);
    SEQAN_ASSERT(empty(g) == false);
    TVertexDescriptor childC1 = addChild(g, rootV);
    String<TVertexDescriptor> leaves;
    collectLeaves(g, rootV, leaves);
    TEdgeDescriptor childC1e = findEdge(g, rootV, childC1);
    SEQAN_ASSERT(_getVertexString(g)[0] == childC1e);
    SEQAN_ASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2);
    SEQAN_ASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 2);
    SEQAN_ASSERT(targetVertex(g, childC1e) == childC1); // Target in a tree = child
    SEQAN_ASSERT(sourceVertex(g, childC1e) == rootV);  // Source in a tree = parent
    SEQAN_ASSERT(childVertex(g, childC1e) == childC1);  // Shortcuts
    SEQAN_ASSERT(parentVertex(g, childC1e) == rootV);
    SEQAN_ASSERT(parentVertex(g, childC1) == rootV);
    SEQAN_ASSERT(empty(g) == false);
    SEQAN_ASSERT(outDegree(g, rootV) == 1);
    TVertexDescriptor childC2 = addChild(g, rootV);
    TVertexDescriptor childC3 = addChild(g, rootV);
    clear(leaves);
    collectLeaves(g, rootV, leaves);
    SEQAN_ASSERT(length(leaves) == 3);
    SEQAN_ASSERT(outDegree(g, rootV) == 3);
    SEQAN_ASSERT(childC1 == 1);
    SEQAN_ASSERT(childC2 == 2);
    SEQAN_ASSERT(childC3 == 3);
    TVertexDescriptor childC2C1 = addChild(g, childC2);
    TVertexDescriptor childC2C1C1 = addChild(g, childC2C1);
    TVertexDescriptor childC2C1C1C1 = addChild(g, childC2C1C1);
    (void)childC2C1C1C1;
    TVertexDescriptor childC2C1C1C2 = addChild(g, childC2C1C1);
    (void)childC2C1C1C2;
    TVertexDescriptor childC4 = addChild(g, rootV);
    (void)childC4;
    SEQAN_ASSERT(inDegree(g, childC2C1) == 1);
    SEQAN_ASSERT(outDegree(g, childC2C1) == 1);
    SEQAN_ASSERT(degree(g, childC2C1) == 2);
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    SEQAN_ASSERT(numTreeEdges(g) == numVertices(g) - 1);
    TEdgeDescriptor childC2C1C1e = findEdge(g, childC2C1C1, childC2C1);

    // Output
    std::stringstream sstream;
    sstream << g;
    char const * EXPECTED =
        "Adjacency list:\n"
        "0 -> 8,3,2,1,\n"
        "1 -> \n"
        "2 -> 4,\n"
        "3 -> \n"
        "4 -> 5,\n"
        "5 -> 7,6,\n"
        "6 -> \n"
        "7 -> \n"
        "8 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 8 (Id: 8)\n"
        "Source: 0,Target: 3 (Id: 3)\n"
        "Source: 0,Target: 2 (Id: 2)\n"
        "Source: 0,Target: 1 (Id: 1)\n"
        "Source: 2,Target: 4 (Id: 4)\n"
        "Source: 4,Target: 5 (Id: 5)\n"
        "Source: 5,Target: 7 (Id: 7)\n"
        "Source: 5,Target: 6 (Id: 6)\n";
    SEQAN_ASSERT(EXPECTED == sstream.str());

    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == 0);
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[4] == 2);
    SEQAN_ASSERT(g.data_parent[5] == 4);
    SEQAN_ASSERT(g.data_parent[6] == 5);
    SEQAN_ASSERT(g.data_parent[7] == 5);
    SEQAN_ASSERT(g.data_parent[8] == 0);
    _rebuildParentMap(g);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == 0);
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[4] == 2);
    SEQAN_ASSERT(g.data_parent[5] == 4);
    SEQAN_ASSERT(g.data_parent[6] == 5);
    SEQAN_ASSERT(g.data_parent[7] == 5);
    SEQAN_ASSERT(g.data_parent[8] == 0);
    SEQAN_ASSERT(childVertex(g, childC2C1C1e) == childC2C1C1);
    SEQAN_ASSERT(parentVertex(g, childC2C1C1e) == childC2C1);
    removeChild(g, rootV, childC2);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[4] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[5] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[6] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[7] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[8] == 0);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(numVertices(g) == 4);
    SEQAN_ASSERT(empty(g) == false);
    SEQAN_ASSERT(inDegree(g, rootV) == 0);
    SEQAN_ASSERT(outDegree(g, rootV) == 3);
    SEQAN_ASSERT(degree(g, rootV) == 3);
    childC2 = addChild(g, rootV);
    childC2C1 = addChild(g, childC2);
    childC2C1C1 = addChild(g, childC2C1);
    childC2C1C1C1 = addChild(g, childC2C1C1);
    childC2C1C1C2 = addChild(g, childC2C1C1);
    removeAllChildren(g, rootV);
    SEQAN_ASSERT(empty(g) == false);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(numTreeEdges(g) == 0);
    SEQAN_ASSERT(numVertices(g) == 1); // Just the root
    SEQAN_ASSERT(inDegree(g, rootV) == 0);
    SEQAN_ASSERT(outDegree(g, rootV) == 0);
    SEQAN_ASSERT(degree(g, rootV) == 0);
    addChild(g, rootV); addChild(g, rootV);
    SEQAN_ASSERT(empty(g) == false);
    SEQAN_ASSERT(numEdges(g) == 2);
    clearEdges(g);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(numVertices(g) == 3);
    SEQAN_ASSERT(empty(g) == false);
    addChild(g, rootV); addChild(g, rootV);
    clearVertices(g);
    SEQAN_ASSERT(empty(g) == true);
    SEQAN_ASSERT(numVertices(g) == 0);
    SEQAN_ASSERT(numEdges(g) == 0);
    createRoot(g);
    childC1 = addChild(g, rootV);
    SEQAN_ASSERT(empty(g) == false);
    SEQAN_ASSERT(numEdges(g) == 1);
    childC3 = addChild(g, rootV);
    childC2 = addChild(g, rootV);
    childC2C1 = addChild(g, childC2);
    childC2C1C1 = addChild(g, childC2C1);
    childC2C1C1C1 = addChild(g, childC2C1C1);
    childC2C1C1C2 = addChild(g, childC2C1C1);
    childC4 = addChild(g, rootV);
    String<unsigned int> mat;   
    // Adjacency matrix
    getAdjacencyMatrix(g, mat);
    unsigned int len = (unsigned int) std::sqrt((double) length(mat));
    SEQAN_ASSERT(getValue(mat, 0 * len + 8) == 1);
    SEQAN_ASSERT(getValue(mat, 8 * len + 0) == 0);
    SEQAN_ASSERT(getValue(mat, 3 * len + 0) == 0);
    SEQAN_ASSERT(getValue(mat, 0 * len + 3) == 1);
    SEQAN_ASSERT(getValue(mat, 0 * len + 4) == 0);
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    // Vertex Adjacency vectors
    String<unsigned int> vectIn, vectOut;
    getVertexAdjacencyVector(vectIn, vectOut, g, 5u);
    SEQAN_ASSERT(length(vectIn) == 1);
    SEQAN_ASSERT(length(vectOut) == 2);
    SEQAN_ASSERT(vectIn[0] == 4);
    SEQAN_ASSERT(vectOut[0] == 7);
    SEQAN_ASSERT(vectOut[1] == 6);
    // Transpose the graph
    transpose(g);
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    TTree g_copy(g);
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    clear(g_copy);
    g_copy = g;
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    transpose(g, g_copy);
    g = g_copy;
    _rebuildParentMap(g);
    SEQAN_ASSERT(numEdges(g) == 8);
    SEQAN_ASSERT(numVertices(g) == 9);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == 0);
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[4] == 3);
    SEQAN_ASSERT(g.data_parent[5] == 4);
    SEQAN_ASSERT(g.data_parent[6] == 5);
    SEQAN_ASSERT(g.data_parent[7] == 5);
    SEQAN_ASSERT(g.data_parent[8] == 0);
    removeOutEdges(g, childC2C1C1);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == 0);
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[4] == 3);
    SEQAN_ASSERT(g.data_parent[5] == 4);
    SEQAN_ASSERT(g.data_parent[6] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[7] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[8] == 0);
    SEQAN_ASSERT(numVertices(g) == 9);
    SEQAN_ASSERT(numEdges(g) == 6);
    removeVertex(g, childC2C1);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == 0);
    SEQAN_ASSERT(g.data_parent[2] == 0);
    SEQAN_ASSERT(g.data_parent[3] == 0);
    SEQAN_ASSERT(g.data_parent[5] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[6] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[7] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[8] == 0);
    SEQAN_ASSERT(numEdges(g) == 4);

    SEQAN_ASSERT(numVertices(g) == 8);
    removeInEdges(g, childC2);
    SEQAN_ASSERT(numEdges(g) == 3);
    SEQAN_ASSERT(numVertices(g) == 8);
    removeOutEdges(g, rootV);
    SEQAN_ASSERT(g.data_parent[0] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[1] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[2] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[3] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[5] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[6] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[7] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(g.data_parent[8] == getNil<TVertexDescriptor>());
    SEQAN_ASSERT(numVertices(g) == 8);
    SEQAN_ASSERT(numEdges(g) == 0);
    SEQAN_ASSERT(empty(g) == false);
    addVertex(g);
    TEdgeDescriptor my_edge = addEdge(g, 0, 1);
    removeEdge(g, my_edge);

//____________________________________________________________________________
// Tree with cargo

    typedef Pair<char, int> TPair;
    typedef Tree<TPair> TEdges;
    typedef Graph<TEdges> TCargoGraph;
    typedef VertexDescriptor<TCargoGraph>::Type TVertexDescriptor2;
    typedef EdgeDescriptor<TCargoGraph>::Type TEdgeDescriptor2;

    TCargoGraph g2;
    createRoot(g2);
    SEQAN_ASSERT(numVertices(g2) == 1);
    SEQAN_ASSERT(numEdges(g2) == 0);
    TVertexDescriptor2 ver1 = addChild(g2, getRoot(g2), TPair('a', 3));
    SEQAN_ASSERT(numChildren(g2, getRoot(g2)) == 1);
    SEQAN_ASSERT(ver1 == 1);
    SEQAN_ASSERT(numVertices(g2) == 2);
    TVertexDescriptor2 ver2 = addChild(g2, getRoot(g2));
    SEQAN_ASSERT(ver2 == 2);
    SEQAN_ASSERT(numVertices(g2) == 3);
    TEdgeDescriptor2 ed1 = findEdge(g2, getRoot(g2), ver1);
    TEdgeDescriptor2 ed2 = findEdge(g2, getRoot(g2), ver2);
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'a');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 3);
    SEQAN_ASSERT(targetVertex(g2, ed1) == ver1);
    SEQAN_ASSERT(sourceVertex(g2, ed1) == getRoot(g2));
    SEQAN_ASSERT(numEdges(g2) == 2);
    assignCargo(ed2, TPair('b', 4));
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'a');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 3);
    SEQAN_ASSERT((getCargo(ed2)).i1 == 'b');
    SEQAN_ASSERT((getCargo(ed2)).i2 == 4);
    cargo(ed1) = TPair('c', 1);
    SEQAN_ASSERT((getCargo(ed1)).i1 == 'c');
    SEQAN_ASSERT((getCargo(ed1)).i2 == 1);
    assignRoot(g2, 1);
    SEQAN_ASSERT(getRoot(g2) == 1);

    //// File read
    //fstream strmKnut;
    //TTree gKnut;
    //strmKnut.open(TEST_PATH "my_tree2.dot", ios_base::in);
    //String<String<char> > nodeMap;
    //String<String<char> > edgeMap;
    //read(strmKnut, gKnut, nodeMap, edgeMap, DotDrawing());
    //strmKnut.close();

    //assignRoot(gKnut, 26);
    //typedef Iterator<TTree, DfsPreorder>::Type TDfsPreorder;
    //TDfsPreorder dfsIt(gKnut, 26);
    //for(;!atEnd(dfsIt);++dfsIt) {
    //	std::cout << *dfsIt << ": ";
    //	std::cout << getProperty(nodeMap, *dfsIt) << std::endl;
    //}
    //std::cout << gKnut << std::endl;
}

SEQAN_DEFINE_TEST(test_graph_types_types_hmm)
{
    typedef double TProbability;
    typedef Dna TAlphabet;
    typedef Size<TAlphabet>::Type TSize;
    typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
    typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
    typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
    TSize alph_size = ValueSize<TAlphabet>::VALUE;

    Dna dnaA = Dna('A');
    Dna dnaC = Dna('C');
    Dna dnaG = Dna('G');
    Dna dnaT = Dna('T');

    // Create an empty HMM
    THmm hmm;
    SEQAN_ASSERT(numVertices(hmm) == 0);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    SEQAN_ASSERT(empty(hmm) == true);
    clearEdges(hmm);
    clearVertices(hmm);
    SEQAN_ASSERT(numVertices(hmm) == 0);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    SEQAN_ASSERT(empty(hmm) == true);

    // Add state1
    TVertexDescriptor state1 = addVertex(hmm);
    SEQAN_ASSERT(length(_getVertexString(hmm)) == 1);
    SEQAN_ASSERT(empty(hmm) == false);
    SEQAN_ASSERT(outDegree(hmm, state1) == 0);
    SEQAN_ASSERT(inDegree(hmm, state1) == 0);
    SEQAN_ASSERT(degree(hmm, state1) == 0);
    SEQAN_ASSERT(numVertices(hmm) == 1);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    emissionProbability(hmm, state1, dnaA) = 0.2;
    SEQAN_ASSERT(getEmissionProbability(hmm, state1, dnaA) == 0.2);
    emissionProbability(hmm, state1, dnaC) = 0.2;
    emissionProbability(hmm, state1, dnaG) = 0.3;
    emissionProbability(hmm, state1, dnaT) = 0.3;
    SEQAN_ASSERT(getEmissionProbability(hmm, state1, dnaA) == 0.2);
    SEQAN_ASSERT(getEmissionProbability(hmm, state1, dnaG) == 0.3);

    // Add state2
    String<TProbability> emis;
    resize(emis, alph_size);
    value(emis, ordValue(dnaA)) = 0.5;
    value(emis, ordValue(dnaC)) = 0.5;
    value(emis, ordValue(dnaG)) = 0.0;
    value(emis, ordValue(dnaT)) = 0.0;
    TVertexDescriptor state2 = addVertex(hmm, emis);
    SEQAN_ASSERT(numVertices(hmm) == 2);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    SEQAN_ASSERT(getEmissionProbability(hmm, state2, dnaC) == 0.5);

    // Add state3
    TVertexDescriptor state3 = addVertex(hmm, emis);
    assignEmissionProbability(hmm, state3, dnaA, 0.3);
    assignEmissionProbability(hmm, state3, dnaC, 0.3);
    assignEmissionProbability(hmm, state3, dnaG, 0.2);
    assignEmissionProbability(hmm, state3, dnaT, 0.2);
    SEQAN_ASSERT(numVertices(hmm) == 3);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    SEQAN_ASSERT(getEmissionProbability(hmm, state3, dnaC) == 0.3);

    // Add edges (transitions)
    TEdgeDescriptor e = addEdge(hmm, state1, state1, 0.95);
    SEQAN_ASSERT(numEdges(hmm) == 1);
    removeEdge(hmm, e);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    e = addEdge(hmm, state1, state1, 0.95);
    SEQAN_ASSERT(numEdges(hmm) == 1);
    removeOutEdges(hmm, state1);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    e = addEdge(hmm, state1, state1, 0.95);
    removeEdge(hmm, sourceVertex(hmm, e), targetVertex(hmm, e));
    SEQAN_ASSERT(numEdges(hmm) == 0);
    e = addEdge(hmm, state1, state1, 0.95);
    removeInEdges(hmm, state1);
    SEQAN_ASSERT(numEdges(hmm) == 0);
    e = addEdge(hmm, state1, state1, 0.95);
    SEQAN_ASSERT(outDegree(hmm, state1) == 1);
    SEQAN_ASSERT(inDegree(hmm, state1) == 1);
    SEQAN_ASSERT(degree(hmm, state1) == 2);
    e = addEdge(hmm, state1, state3);
    assignTransitionProbability(hmm, e, 0.05);
    THmm hmm_tr = hmm;
    transpose(hmm_tr);
    SEQAN_ASSERT(getTransitionProbability(hmm_tr, state3, state1) == 0.05);
    clear(hmm_tr);
    transpose(hmm, hmm_tr);
    SEQAN_ASSERT(getTransitionProbability(hmm_tr, state3, state1) == 0.05);
    e = addEdge(hmm, state3, state3);
    transitionProbability(hmm, e) = 0.4;
    SEQAN_ASSERT(getTransitionProbability(hmm, e) == 0.4);
    e = addEdge(hmm, state3, state1);
    transitionProbability(hmm, state3, state1) = 0.1;
    e = addEdge(hmm, state2, state2);
    assignTransitionProbability(hmm, state2, state2, 1.0);
    SEQAN_ASSERT(numVertices(hmm) == 3);
    SEQAN_ASSERT(numEdges(hmm) == 5);
    SEQAN_ASSERT(getTransitionProbability(hmm, state3, state1) == 0.1);

    // Add begin and end state
    TVertexDescriptor begState = addVertex(hmm);
    TVertexDescriptor eState = addVertex(hmm);
    addEdge(hmm, begState, state1, 1.0);
    addEdge(hmm, state3, eState, 0.5);
    addEdge(hmm, eState, eState, 1.0);
    beginState(hmm) = state3;
    SEQAN_ASSERT(numVertices(hmm) == 5);
    SEQAN_ASSERT(getBeginState(hmm) == state3);
    assignBeginState(hmm, begState);
    SEQAN_ASSERT(getBeginState(hmm) == begState);
    endState(hmm) = state3;
    SEQAN_ASSERT(getEndState(hmm) == state3);
    assignEndState(hmm, eState);
    SEQAN_ASSERT(getEndState(hmm) == eState);

    // Output
    std::stringstream sstream;
    sstream << hmm;
    char const * EXPECTED =
        "Alphabet:\n"
        "{A,C,G,T}\n"
        "States:\n"
        "{0,1,2,3 (Silent),4 (Silent)}\n"
        "Begin state: 3\n"
        "End state: 4\n"
        "Transition probabilities:\n"
        "0 -> 2 (0.05) ,0 (0.95) \n"
        "1 -> 1 (1) \n"
        "2 -> 4 (0.5) ,0 (0.1) ,2 (0.4) \n"
        "3 -> 0 (1) \n"
        "4 -> 4 (1) \n"
        "Emission probabilities:\n"
        "0: A (0.2) ,C (0.2) ,G (0.3) ,T (0.3) \n"
        "1: A (0.5) ,C (0.5) ,G (0) ,T (0) \n"
        "2: A (0.3) ,C (0.3) ,G (0.2) ,T (0.2) ";
    SEQAN_ASSERT_EQ(EXPECTED, sstream.str());

    // Change model
    removeVertex(hmm, state2);
    THmm hmm_copy(hmm);
    SEQAN_ASSERT(numVertices(hmm_copy) == 4);
    SEQAN_ASSERT(getBeginState(hmm_copy) == begState);
    SEQAN_ASSERT(getEndState(hmm_copy) == eState);
    clear(hmm_copy);
    SEQAN_ASSERT(numVertices(hmm_copy) == 0);
    hmm_copy = hmm;
    SEQAN_ASSERT(numVertices(hmm_copy) == 4);
    SEQAN_ASSERT(getBeginState(hmm_copy) == begState);
    SEQAN_ASSERT(getEndState(hmm_copy) == eState);
    SEQAN_ASSERT(idCount(_getEdgeIdManager(hmm_copy)) == 7);

    // Test silent states
    TVertexDescriptor testState1 = addVertex(hmm, emis, true);
    TVertexDescriptor testState2 = addVertex(hmm, emis, false);
    SEQAN_ASSERT(isSilent(hmm, testState1) != isSilent(hmm, testState2));
    assignSilentStatus(hmm, testState1, false);
    SEQAN_ASSERT(isSilent(hmm, testState1) == isSilent(hmm, testState2));
    silentStatus(hmm, testState1) = true;
    SEQAN_ASSERT(isSilent(hmm, testState1) != isSilent(hmm, testState2));
}

SEQAN_BEGIN_TESTSUITE(test_graph_types_property_map)
{
    SEQAN_CALL_TEST(test_graph_types_types_directed);
    SEQAN_CALL_TEST(test_graph_types_types_undirected);
    SEQAN_CALL_TEST(test_graph_types_types_automaton);
    SEQAN_CALL_TEST(test_graph_types_types_word_graph);
    SEQAN_CALL_TEST(test_graph_types_types_word_tree);
    SEQAN_CALL_TEST(test_graph_types_types_hmm);
}
SEQAN_END_TESTSUITE
