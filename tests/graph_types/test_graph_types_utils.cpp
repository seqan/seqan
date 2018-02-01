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

template <typename TGraph>
void Test_GraphDrawing_Tmp(char const * EXPECTED)
{
    // Create a dummy graph
    TGraph g;
    addVertex(g); addVertex(g); addEdge(g, 0, 1);
    addVertex(g); addEdge(g, 0, 2); addVertex(g);
    addVertex(g); addVertex(g); addEdge(g, 3, 4);

    // Dot Drawing
    std::stringstream sstream;
    sstream << g;
    SEQAN_ASSERT_EQ(EXPECTED, sstream.str());
}

SEQAN_DEFINE_TEST(test_graph_types_utils_graph_drawing)
{
    typedef Position<String<char> >::Type TPosition;

    char const * EXPECTED_DIRECTED =
        "Adjacency list:\n"
        "0 -> 2,1,\n"
        "1 -> \n"
        "2 -> \n"
        "3 -> 4,\n"
        "4 -> \n"
        "5 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 2 (Id: 1)\n"
        "Source: 0,Target: 1 (Id: 0)\n"
        "Source: 3,Target: 4 (Id: 2)\n";
    char const * EXPECTED_UNDIRECTED =
        "Adjacency list:\n"
        "0 -> 2,1,\n"
        "1 -> 0,\n"
        "2 -> 0,\n"
        "3 -> 4,\n"
        "4 -> 3,\n"
        "5 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 2 (Id: 1)\n"
        "Source: 0,Target: 1 (Id: 0)\n"
        "Source: 3,Target: 4 (Id: 2)\n";
    char const * EXPECTED_TREE =
        "Adjacency list:\n"
        "0 -> 2,1,\n"
        "1 -> \n"
        "2 -> \n"
        "3 -> 4,\n"
        "4 -> \n"
        "5 -> \n"
        "Edge list:\n"
        "Source: 0,Target: 2 (Id: 2)\n"
        "Source: 0,Target: 1 (Id: 1)\n"
        "Source: 3,Target: 4 (Id: 4)\n";
    Test_GraphDrawing_Tmp<Graph<Directed<> > >(EXPECTED_DIRECTED);
    Test_GraphDrawing_Tmp<Graph<Undirected<> > >(EXPECTED_UNDIRECTED);
    Test_GraphDrawing_Tmp<Graph<Tree<> > >(EXPECTED_TREE);

    // Automat
    Graph<Automaton<Dna> > automat;
    createRoot(automat); addVertex(automat); addEdge(automat, 0, 1, 'a');
    addVertex(automat); addEdge(automat, 0, 2, 'g');
    // Dot Drawing

    {
        std::stringstream sstream;
        writeRecords(sstream, automat, DotDrawing());
        char const * EXPECTED =
            "digraph G {\n"
            "\n"
            "/* Graph Attributes */\n"
            "graph [rankdir = LR];\n"
            "\n"
            "/* Node Attributes */\n"
            "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n"
            "\n"
            "/* Edge Attributes */\n"
            "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n"
            "\n"
            "/* Nodes */\n"
            "0 [label = \"0\", shape = doublecircle];\n"
            "1 [label = \"1\"];\n"
            "2 [label = \"2\"];\n"
            "\n"
            "/* Edges */\n"
            "0 -> 1 [label = \"A\"];\n"
            "0 -> 2 [label = \"G\"];\n"
            "\n"
            "}\n";
        SEQAN_ASSERT_EQ(EXPECTED, sstream.str());
    }

    // Trie
    Graph<Automaton<char> > trie;
    String<String<TPosition> > pos;
    String<String<char> > keywords;
    appendValue(keywords, String<char>("announce"));
    appendValue(keywords, String<char>("annual"));
    appendValue(keywords, String<char>("annually"));
    createTrie(trie, pos, keywords);
    // Dot Drawing
    String<String<char> > nodeMap;
    _createTrieNodeAttributes(trie, pos, nodeMap);
    String<String<char> > edgeMap;
    _createEdgeAttributes(trie, edgeMap);
    {
        std::stringstream sstream;
        char const * EXPECTED =
            "digraph G {\n"
            "\n"
            "/* Graph Attributes */\n"
            "graph [rankdir = LR];\n"
            "\n"
            "/* Node Attributes */\n"
            "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n"
            "\n"
            "/* Edge Attributes */\n"
            "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n"
            "\n"
            "/* Nodes */\n"
            "0 [label = \"0\", shape = doublecircle];\n"
            "1 [label = \"1\"];\n"
            "2 [label = \"2\"];\n"
            "3 [label = \"3\"];\n"
            "4 [label = \"4\"];\n"
            "5 [label = \"5\"];\n"
            "6 [label = \"6\"];\n"
            "7 [label = \"7\"];\n"
            "8 [shape = box, label = \"8 {0}\"];\n"
            "9 [label = \"9\"];\n"
            "10 [label = \"10\"];\n"
            "11 [shape = box, label = \"11 {1}\"];\n"
            "12 [label = \"12\"];\n"
            "13 [shape = box, label = \"13 {2}\"];\n"
            "\n"
            "/* Edges */\n"
            "0 -> 1 [label = \"a\"];\n"
            "1 -> 2 [label = \"n\"];\n"
            "2 -> 3 [label = \"n\"];\n"
            "3 -> 4 [label = \"o\"];\n"
            "3 -> 9 [label = \"u\"];\n"
            "4 -> 5 [label = \"u\"];\n"
            "5 -> 6 [label = \"n\"];\n"
            "6 -> 7 [label = \"c\"];\n"
            "7 -> 8 [label = \"e\"];\n"
            "9 -> 10 [label = \"a\"];\n"
            "10 -> 11 [label = \"l\"];\n"
            "11 -> 12 [label = \"l\"];\n"
            "12 -> 13 [label = \"y\"];\n"
            "\n"
            "}\n";
        writeRecords(sstream, trie, nodeMap, edgeMap, DotDrawing());
        SEQAN_ASSERT_EQ(EXPECTED, sstream.str());
    }

    // WordGraph
    typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
    TWordGraph wordGr;
    addVertex(wordGr); addVertex(wordGr); addVertex(wordGr); addVertex(wordGr);
    assignRoot(wordGr, 3); root(wordGr) = 2; addEdge(wordGr, 0, 3, "ag");
    addVertex(wordGr); addVertex(wordGr); addEdge(wordGr, 0, 5, "g");
    addVertex(wordGr); addVertex(wordGr); addEdge(wordGr, 3, 1, "aggg");
    addEdge(wordGr, 3, 4, "gg"); addEdge(wordGr, 5, 2, "aggg"); addEdge(wordGr, 5, 7, "g");
    addEdge(wordGr, 7, 6, "g"); assignRoot(wordGr, 0);
    {
        std::stringstream sstream;
        char const * EXPECTED =
            "digraph G {\n"
            "\n"
            "/* Graph Attributes */\n"
            "graph [rankdir = LR];\n"
            "\n"
            "/* Node Attributes */\n"
            "node [shape = rectangle, fillcolor = white, style = filled, fontname = \"Times-Italic\"];\n"
            "\n"
            "/* Edge Attributes */\n"
            "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n"
            "\n"
            "/* Nodes */\n"
            "0 [label = \"0\", shape = doublecircle];\n"
            "1 [label = \"1\"];\n"
            "2 [label = \"2\"];\n"
            "3 [label = \"3\"];\n"
            "4 [label = \"4\"];\n"
            "5 [label = \"5\"];\n"
            "6 [label = \"6\"];\n"
            "7 [label = \"7\"];\n"
            "\n"
            "/* Edges */\n"
            "0 -> 3 [label = \"AG\"];\n"
            "0 -> 5 [label = \"G\"];\n"
            "3 -> 1 [label = \"AGGG\"];\n"
            "3 -> 4 [label = \"GG\"];\n"
            "5 -> 2 [label = \"AGGG\"];\n"
            "5 -> 7 [label = \"G\"];\n"
            "7 -> 6 [label = \"G\"];\n"
            "\n"
            "}\n";
        writeRecords(sstream, wordGr, DotDrawing());
        SEQAN_ASSERT_EQ(EXPECTED, sstream.str());
    }
}

SEQAN_BEGIN_TESTSUITE(test_graph_types_utils)
{
    SEQAN_CALL_TEST(test_graph_types_utils_graph_drawing);
}
SEQAN_END_TESTSUITE
