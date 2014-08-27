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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_CONTIG_GRAPH_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_CONTIG_GRAPH_H_

#include <iosfwd>
#include <map>
#include <vector>

#include <seqan/sequence.h>

#include <lemon/smart_graph.h>

#include "asm/frag_store.h"

// ============================================================================
// External Forwards
// ============================================================================

namespace assembler { struct ContigGraph; }

namespace rep_solv {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class LibraryInfo
// ----------------------------------------------------------------------------

// Stores information about a library.

struct LibraryInfo
{
    // Mean and standard deviation of template size.
    double mean { 0 };
    double sd { 0 };
    // Allow (delta * sd) for deviation.
    int k { 3 };

    LibraryInfo() = default;
    LibraryInfo(double mean, double sd, int k = 3) : mean(mean), sd(sd), k(k) {}
};

// ----------------------------------------------------------------------------
// Class MateInfo
// ----------------------------------------------------------------------------

// TODO(holtgrew): Move to module asm/assembler, as well as MateInfos.

// Information relative for one mate-pair based link.

struct MateInfo
{
    static const unsigned INVALID = (unsigned)-1;

    // ID of left/right read.
    unsigned readL { INVALID };
    unsigned readR { INVALID };
    // Left/right contigID.
    unsigned leftID { INVALID };
    unsigned rightID { INVALID };
    // Start position in left and end position in right.
    int beginPos { 0 };
    int endPos { 0 };
    // Length of left/right alignment.
    int leftLen { 0 };
    int rightLen { 0 };
    // ID of the library.
    unsigned libraryID { INVALID };
    // Number of copies for this link.
    int numCopies { 1 };

    MateInfo() = default;
    MateInfo(unsigned readL, unsigned readR, unsigned leftID, unsigned rightID, int beginPos, int endPos,
             int leftLen, int rightLen, unsigned libraryID, int numCopies = 1) :
            readL(readL), readR(readR), leftID(leftID), rightID(rightID), beginPos(beginPos), endPos(endPos),
            leftLen(leftLen), rightLen(rightLen), libraryID(libraryID), numCopies(numCopies)
    {}

    // Return weight of link, is 1/numCopies.
    double weight() const { return (1.0 / numCopies); }
};

bool lt(MateInfo const & lhs, MateInfo const & rhs);

std::ostream & operator<<(std::ostream & out, MateInfo const & info);

// ----------------------------------------------------------------------------
// Class MateInfos
// ----------------------------------------------------------------------------

class MateInfos
{
public:
    void insert(MateInfo const & info);
    void print(std::ostream & out);

    std::map<unsigned, std::vector<unsigned>> readToRecords;
    std::vector<MateInfo> records;

    std::vector<LibraryInfo> libraries;
};

// ----------------------------------------------------------------------------
// Class ContigLabel
// ----------------------------------------------------------------------------

// Label for ContigGraph nodes.

struct ContigLabel
{
    static const unsigned INVALID = (unsigned)-1;

    // The contig ID.
    unsigned id { INVALID };
    // The length of the contig.
    unsigned length { 0 };
    // The number of concordant mates on the contig.
    unsigned mateCount { 0 };

    ContigLabel() = default;
    ContigLabel(unsigned id, unsigned length, unsigned mateCount) :
            id(id), length(length), mateCount(mateCount)
    {}
};

// ----------------------------------------------------------------------------
// Class EdgeLabel
// ----------------------------------------------------------------------------

// Label for ContigGraph edges.

struct EdgeLabel
{
    struct FuzzyDistance
    {
        double mean { 0 };
        double sd { 0 };

        FuzzyDistance() = default;
        FuzzyDistance(double mean, double sd) : mean(mean), sd(sd) {}
    };

    static const unsigned INVALID = (unsigned)-1; // 4294967295
    static const unsigned SOURCE = (unsigned)-2;  // 4294967294
    static const unsigned TARGET = (unsigned)-3;  // 4294967293

    // Contig ID of left and right contig for this edge.
    unsigned leftID { INVALID };
    unsigned rightID { INVALID };
    // Number of links, duplicate links count as 1 together (i.e. 1/k for k copies).
    double uniqueLinks { 0.0 };
    // Count of links, duplicate links count as 1 each.
    double duplicateLinks { 0.0 };

    // Distance from overlap and unique/duplicate links.
    FuzzyDistance overlapDistance;
    FuzzyDistance uniqueLinkDistance;
    FuzzyDistance dupLinkDistance;

    EdgeLabel() = default;
    EdgeLabel(unsigned leftID, unsigned rightID, double uniqueLinks = 0.0, double duplicateLinks = 0.0) :
            leftID(leftID), rightID(rightID), uniqueLinks(uniqueLinks), duplicateLinks(duplicateLinks)
    {}
};

// ----------------------------------------------------------------------------
// Class ContigGraph
// ----------------------------------------------------------------------------

// Contig graph used for repeat resolution.

class ContigGraph
{
public:
    // Artificial s and t vertices, no entry in nodes.
    lemon::SmartGraph::Node s, t;

    // The underlying graph structure.
    lemon::SmartGraph graph;
    // Label for the nodes.
    lemon::SmartGraph::NodeMap<ContigLabel> contig;
    // Label for the edges.
    lemon::SmartGraph::EdgeMap<EdgeLabel> link;
    // Mapping from contigID to node in graph.
    std::vector<lemon::SmartGraph::Node> node;

    // Mapping from read to all contigs it is on.
    std::map<unsigned, std::vector<unsigned>> readToContigs;

    ContigGraph() : contig(graph), link(graph) { init(); }
    ContigGraph(ContigGraph const & other);

    // Print graph to out.
    void print(std::ostream & out) const;

    // Add edge to the graph and assign the given label.
    void addEdge(lemon::SmartGraph::Node u, lemon::SmartGraph::Node v, EdgeLabel const & label)
    { link[graph.addEdge(u, v)] = label; }

private:

    void init()
    {
        s = graph.addNode();
        t = graph.addNode();
        contig[s].id = EdgeLabel::SOURCE;
        contig[t].id = EdgeLabel::TARGET;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function buildContigGraph()
// ----------------------------------------------------------------------------

// Build repeat resolution ContigGraph from assembly ContigGraph.  All added edges will be overlap edges.

void buildContigGraph(ContigGraph & out,
                      assembler::ContigGraph const & in,
                      TFragmentStore const & fragStore);

// ----------------------------------------------------------------------------
// Function addPairedEndLinks()
// ----------------------------------------------------------------------------

// Adds paired-end links to the ContigGraph.
//
// If contig sequences are given then we remove containment links and links not verified by overlapper.
//
// TODO(holtgrew): This would also be a good place to remove links with conflicts?

void addPairedEndLinks(ContigGraph & cg,
                       MateInfos const & mateInfos,
                       seqan::StringSet<seqan::Dna5String> const * contigSeqs = nullptr);

// ----------------------------------------------------------------------------
// Function removeContigs()
// ----------------------------------------------------------------------------

// Remove the contigs that are flagged so in doRemove.
//
// We need to write the graph to a copy because lemon::SmartGraph cannot have artifacts removed.
//
// Note that this can yield unmapped reads (marked by (out.readToContigs.count(readID) == 0u)).

void removeContigs(ContigGraph & out,
                   ContigGraph const & in,
                   lemon::SmartGraph::NodeMap<bool> const & doRemove);

// ----------------------------------------------------------------------------
// Function stPathExists()
// ----------------------------------------------------------------------------

// Returns true if an s-t path exists.

bool stPathExists(ContigGraph const & in);

}  // namespace rep_solv

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SOLV_CONTIG_GRAPH_H_
