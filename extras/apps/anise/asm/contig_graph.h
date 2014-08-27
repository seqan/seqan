// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_TEST_CONTIG_GRAPH_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_TEST_CONTIG_GRAPH_H_

// TODO(holtgrew): Keep contig length in contig.

#include <iosfwd>
#include <vector>
#include <map>
#include <utility>

#include <lemon/smart_graph.h>

#include <seqan/basic.h>

#include "asm/overlapper.h"

namespace assembler {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Contig
// ----------------------------------------------------------------------------

struct Contig
{
    // Approximately positioned read on contig.
    struct PositionedRead
    {
        unsigned readID { 0 };
        int beginPos { 0 };
        int length { 0 };

        PositionedRead() = default;

        PositionedRead(unsigned readID, int beginPos, int length) :
                readID(readID), beginPos(beginPos), length(length)
        {}

        bool operator<(PositionedRead const & other) const
        {
            return makeTuple() < other.makeTuple();
        }

    private:

        std::tuple<int, int, unsigned> makeTuple() const
        {
            return std::make_tuple(beginPos, length, readID);
        }
    };

    // ContigID.
    unsigned id;
    // Positioned reads, sorted by begin position, ties broken by smaller length and read ID.
    std::vector<PositionedRead> posReads;

    explicit Contig(unsigned id = seqan::MaxValue<unsigned>::VALUE) : id(id)
    {}

    void print(std::ostream & out) const;

    // Merge other contig (right of current), overlapLength is length of predicted overlap.
    void merge(Contig const & other, int overlapLength);
};

// ----------------------------------------------------------------------------
// Class ContigGraph
// ----------------------------------------------------------------------------

struct ContigGraph
{
    // The digraph structure of the contig graph.
    lemon::SmartDigraph graph;
    // Label each node with a Contig.
    lemon::SmartDigraph::NodeMap<Contig> contig;
    // Label each edge with predicted overlap length.
    lemon::SmartDigraph::ArcMap<int> overlapLength;
    // The node for each contig.
    std::vector<lemon::SmartDigraph::Node> node;
    // Mapping from read ID to containing contig.
    std::map<unsigned, lemon::SmartDigraph::Node> readToContig;
    // Set of contained read ids.  // TODO(holtgrew): Redundant with readToContig, remove.
    std::set<unsigned> containedReadIDs;

    ContigGraph() : contig(graph), overlapLength(graph)
    {}

    // Adds contig and returns contig ID.
    unsigned addContig(Contig && utg)
    {
        unsigned contigID = node.size();
        utg.id = contigID;
        auto u = graph.addNode();
        node.push_back(u);
        SEQAN_ASSERT_LT(contigID, node.size());

        for (auto posRead : utg.posReads)
            readToContig[posRead.readID] = u;
        contig[u] = std::move(utg);

        return contigID;
    }

    void print(std::ostream & out) const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace assembler

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_TEST_CONTIG_GRAPH_H_
