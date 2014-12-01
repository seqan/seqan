// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#include "contig_graph.h"

#include <ostream>
#include <string>
#include <sstream>

#include "asm/overlapper.h"

namespace assembler {

// ----------------------------------------------------------------------------
// Class Contig
// ----------------------------------------------------------------------------

std::string posReadStr(Contig::PositionedRead const & posRead)
{
    std::stringstream ss;
    ss << "PositionedRead(readID=" << posRead.readID << ", beginPos=" << posRead.beginPos
       << ", length=" << posRead.length << ")";
    return ss.str();
}

void Contig::print(std::ostream & out) const
{
    out << ",-- UNITIG\n";
    out << "| id: " << id << "\n";
    out << "+-- APPROXIMATE LAYOUT\n";
    for (auto const & posRead : posReads)
        out << "| " << posReadStr(posRead) << "\n";
    out << "`--\n";
}

void Contig::merge(Contig const & other, int overlapLength)
{
    for (auto posRead : other.posReads)
    {
        posRead.beginPos -= overlapLength;
        posReads.push_back(posRead);
    }

    // Resort reads to yield new approximate layout.
    std::sort(posReads.begin(), posReads.end());
}

// ----------------------------------------------------------------------------
// Class ContigGraph
// ----------------------------------------------------------------------------

void ContigGraph::print(std::ostream & out) const
{
    out << "ContigGraph (#contigs=" << node.size() << ")\n"
        << "NODES\n";
    for (auto u : node)
        contig[u].print(out);
    out << "EDGES\n";
    for (lemon::SmartDigraph::ArcIt arc(graph); arc != lemon::INVALID; ++arc)
    {
        out << "\t" << contig[graph.source(arc)].id << " -> " << contig[graph.target(arc)].id << "\n";
    }
}

}  // namespace assembler
