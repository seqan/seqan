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

#include <iostream>
#include <sstream>

#include <seqan/index.h>
#include <seqan/stream.h>

SEQAN_DEFINE_TEST(test_index_drawing_esa_dot)
{
    seqan::CharString myString = "banana";
    seqan::Index<seqan::CharString> stree(myString);
    
    std::stringstream sstream;
    writeRecords(sstream, stree, seqan::DotDrawing());
    
    std::stringstream expected;
    expected << "digraph G {\n"
             << "\n"
             << "/* Graph Attributes */\n"
             << "graph [rankdir = LR];\n"
             << "\n"
             << "/* Node Attributes */\n"
             << "node [shape = ellipse, fillcolor = lightgrey, style = filled, fontname = \"Times-Italic\"];\n"
             << "\n"
             << "/* Edge Attributes */\n"
             << "edge [fontname = \"Times-Italic\", arrowsize = 0.75, fontsize = 16];\n"
             << "\n"
             << "/* Edges */\n"
             << "\"[0:6)\" [style = dashed];\n"
             << "\"[0:3)\";\n"
             << "\"[0:6)\" -> \"[0:3)\" [label = \"a\"];\n"
             << "\"[1:3)\";\n"
             << "\"[0:3)\" -> \"[1:3)\" [label = \"na\"];\n"
             << "\"[2:3)\";\n"
             << "\"[1:3)\" -> \"[2:3)\" [label = \"na\"];\n"
             << "\"[3:4)\";\n"
             << "\"[0:6)\" -> \"[3:4)\" [label = \"banana\"];\n"
             << "\"[4:6)\";\n"
             << "\"[0:6)\" -> \"[4:6)\" [label = \"na\"];\n"
             << "\"[5:6)\";\n"
             << "\"[4:6)\" -> \"[5:6)\" [label = \"na\"];\n"
             << "\n"
             << "}\n";
    SEQAN_ASSERT_EQ(expected.str(), sstream.str());
}

SEQAN_BEGIN_TESTSUITE(test_index_drawing)
{
	SEQAN_CALL_TEST(test_index_drawing_esa_dot);
}
SEQAN_END_TESTSUITE
