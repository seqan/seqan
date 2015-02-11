// ===========================================================================
//                 SGIP - Solution of Graph Isomorphism Problem
// ===========================================================================
// Copyright (C) 2012 by Jialu Hu
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your options) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ===========================================================================
// Author: Jialu Hu <Jialu.Hu@fu-berlin.de>
// ===========================================================================

#ifndef APPS_SGIP_SGIP_H_
#define APPS_SGIP_SGIP_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/graph_types.h>

// ==========================================================================
// Forwards
// ==========================================================================

size_t const MAX_ELEMENT = 1000;
struct SivaLab_;
typedef seqan::Tag<SivaLab_> SivaLab;

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function _readWord()
// --------------------------------------------------------------------------

// Read one word from SivaLab format file.
inline unsigned short _readWord(std::ifstream & in)
{
    unsigned char c1, c2;
    c1 = static_cast<char>(in.get());
    c2 = static_cast<char>(in.get());
    return c1 | (c2 << 8);
}

// --------------------------------------------------------------------------
// Function _createGraph()
// --------------------------------------------------------------------------

// Create graph form various data format, e.g. SivaLab.
template <typename TFilename, typename TSpec, typename TTag>
inline bool _createGraph(seqan::Graph<TSpec> &, seqan::Tag<TTag> const &, TFilename &);

template <typename TSpec>
inline bool _createGraph(seqan::Graph<TSpec> & graph, SivaLab const, char const * filename)
{
    using namespace seqan;

    typedef Graph<TSpec> TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef std::vector<TVertexDescriptor> TString;
    std::ifstream is;
    TVertexDescriptor i, j;
    TString edges;

    is.open(filename, std::ios_base::in | std::ios_base::binary);
    if (!is.is_open())
    {
        std::cerr << "error in open input file" << std::endl;
        return false;
    }
    unsigned short edgeNumber;
    unsigned short nodeNumber;
    nodeNumber = _readWord(is);
    edges.reserve(MAX_ELEMENT);

    for (i = 0; i < nodeNumber; i++)
    {
        edgeNumber = _readWord(is);
        for (j = 0; j < edgeNumber; j++)
        {
            TVertexDescriptor target = _readWord(is);
            edges.push_back(i);
            edges.push_back(target);
        }
    }
    size_t len = length(edges) / 2;
    addEdges(graph, edges, len);
    is.close();
    return true;
}

#endif  // #ifndef APPS_SGIP_SGIP_H_
