// ==========================================================================
//                                 BASIL
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

#include "cluster_matching.h"

#include <lemon/smart_graph.h>
#include <lemon/matching.h>

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Function leftOf()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Move into the SeqAn library?

// Returns true if lhs is truly left of rhs.

bool leftOf(seqan::GenomicRegion const & lhs,
            seqan::GenomicRegion const & rhs)
{
    if (lhs.rID < rhs.rID)
        return true;
    if (lhs.rID > rhs.rID)
        return false;
    if (lhs.endPos <= rhs.beginPos)
        return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function doOverlap()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Move into the SeqAn library?

// Returns true if the two regions overlap.

bool doOverlap(seqan::GenomicRegion const & lhs,
               seqan::GenomicRegion const & rhs)
{
    if (lhs.rID != rhs.rID)
        return false;
    return (lhs.beginPos < rhs.endPos && rhs.beginPos < lhs.endPos);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ClusterMatchingImpl
// ----------------------------------------------------------------------------

// Implementation of the cluster matching.

class ClusterMatchingImpl
{
public:
    ClusterMatchingImpl()
    {}

    // Perform the clustering.
    void run(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > & result,
             std::vector<OeaClusterRecord> const & oeaClusters,
             std::vector<ClippingClusterRecord> const & clippingClusters);

    bool _doOverlap(OeaClusterRecord const & lhs, OeaClusterRecord const & rhs)
    {
        return (rhs.region.beginPos < lhs.region.endPos &&
                lhs.region.beginPos < rhs.region.endPos);
    }

    bool _isInCluster(OeaClusterRecord const & oea, ClippingClusterRecord const & clipping)
    {
        return (oea.region.rID == clipping.region.rID &&
                oea.region.beginPos <= clipping.region.beginPos &&
                oea.region.endPos >= clipping.region.endPos);
    }
};

void ClusterMatchingImpl::run(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > & result,
                              std::vector<OeaClusterRecord> const & oeaClusters,
                              std::vector<ClippingClusterRecord> const & clippingClusters)
{
    typedef lemon::SmartGraph TGraph;
    typedef TGraph::Node      TNode;
    typedef TGraph::Edge      TEdge;

    typedef lemon::SparseMap<TEdge, double> TEdgeWeights;
    typedef lemon::SparseMap<TNode, unsigned> TNodeIdMap;

    // Offsets of first OEA/clipping vertex.
    TNodeIdMap nodeIds;
    std::vector<TGraph::Node> oeaNodes;
    oeaNodes.reserve(oeaClusters.size());
    std::vector<TGraph::Node> clippingNodes;
    clippingNodes.reserve(clippingClusters.size());
    std::vector<TEdge> edges;

    TEdgeWeights edgeWeights;

    // -----------------------------------------------------------------------
    // Build Bipartite Graph
    // -----------------------------------------------------------------------

    // TODO(holtgrew): LEMON 1.2.3 does not support bipartite graphs and algorithms on them yet.
    TGraph g;
    g.reserveNode(oeaClusters.size() + clippingClusters.size());
    // Add vertices.
    for (std::vector<OeaClusterRecord>::const_iterator it = oeaClusters.begin(); it != oeaClusters.end(); ++it)
    {
        oeaNodes.push_back(g.addNode());
        nodeIds[oeaNodes.back()] = oeaNodes.size() - 1;
    }
    for (std::vector<ClippingClusterRecord>::const_iterator it = clippingClusters.begin(); it != clippingClusters.end(); ++it)
    {
        clippingNodes.push_back(g.addNode());
        nodeIds[clippingNodes.back()] = oeaNodes.size() + clippingNodes.size() - 1;
    }

    // Add edges in a sweep-line fashion.  We compute overlapping intervals the OEA clusters and then all clipping
    // clusters that fall into the interval.
    //
    // We need two ranges on both the OEA and the clipping records.
    typedef std::vector<OeaClusterRecord>::const_iterator TOeaIter;
    TOeaIter itOeaBegin = oeaClusters.begin();
    TOeaIter itOeaEnd = oeaClusters.begin();
    typedef std::vector<ClippingClusterRecord>::const_iterator TClippingIter;
    TClippingIter itClippingBegin = clippingClusters.begin();
    TClippingIter itClippingEnd = clippingClusters.begin();
    // Perform sweep-line algorithm.
    // TODO(holtgrew): Handling reference switch.
    while (itOeaBegin != oeaClusters.end() && itClippingBegin != clippingClusters.end())
    {
        // Go over intervals in OEA records.
        for (TOeaIter itOea = itOeaEnd++; itOeaEnd != oeaClusters.end(); ++itOeaEnd, ++itOea)
            if (!_doOverlap(*itOea, *itOeaEnd))
                break;
        // std::cerr << "OEA interval: (size = " << (itOeaEnd - itOeaBegin) << ") " << itOeaBegin->region.rID << "\t"
        //           << (itOeaBegin->region.beginPos + 1) << "\t" << (itOeaEnd - 1)->region.endPos << "\n";
        // for (TOeaIter it = itOeaBegin; it != itOeaEnd; ++it)
        //     std::cerr << "OEA\t" << *it << "\n";

        // Find first clipping record that is not before the OEA interval beginning.
        for (; itClippingBegin != clippingClusters.end() && leftOf(itClippingBegin->region, itOeaBegin->region); ++itClippingBegin)
            continue;  // skipping
        itClippingEnd = itClippingBegin;
        // Find last clipping record that is not behind the OEA interval beginning.
        for (; itClippingEnd != clippingClusters.end() && doOverlap(itClippingEnd->region, (itOeaEnd - 1)->region); ++itClippingEnd)
            continue;  // skipping
        // for (TClippingIter it = itClippingBegin; it != itClippingEnd; ++it)
        //     std::cerr << "CLIPPING\t" << *it << "\n";

        // Process the currently active window.
        for (TOeaIter itOea = itOeaBegin; itOea != itOeaEnd; ++itOea)
        {
            unsigned idOea = itOea - oeaClusters.begin();
            double weightOea = sqrt((1.0 + itOea->leftWeight) * (1.0 + itOea->rightWeight));
            for (TClippingIter itClipping = itClippingBegin; itClipping != itClippingEnd; ++itClipping)
            {
                if (!_isInCluster(*itOea, *itClipping))
                    continue;  // Clipping cluster not in OEA cluster.

                unsigned idClipping = itClipping - clippingClusters.begin();
                double weightClipping = sqrt((1.0 + itClipping->leftWeight) * (1.0 + itClipping->rightWeight));
                edges.push_back(g.addEdge(oeaNodes[idOea], clippingNodes[idClipping]));
                edgeWeights[edges.back()] = weightOea + weightClipping;
            }
        }

        // This overlapping region is done.
        itOeaBegin = itOeaEnd;
        itClippingBegin = itClippingEnd;
    }

    // -----------------------------------------------------------------------
    // Compute Matching
    // -----------------------------------------------------------------------

    lemon::MaxWeightedMatching<TGraph, TEdgeWeights> mwMatching(g, edgeWeights);
    mwMatching.run();

    // -----------------------------------------------------------------------
    // Write Out Result
    // -----------------------------------------------------------------------

    seqan::String<bool> oeaMatches;
    resize(oeaMatches, oeaClusters.size(), false);

    // TODO(holtgrew): We could also write out the candidates not supported by both OEA and clipping signals.
    for (unsigned i = 0; i < oeaClusters.size(); ++i)
    {
        if (mwMatching.mate(oeaNodes[i]) == lemon::INVALID)
            continue;
        unsigned idOea = nodeIds[oeaNodes[i]];
        unsigned idClipping = nodeIds[mwMatching.mate(oeaNodes[i])] - oeaNodes.size();
        result.push_back(std::make_pair(oeaClusters[idOea], clippingClusters[idClipping]));
        oeaMatches[idOea] = true;
    }

    ClippingClusterRecord invalidCCR;
    for (unsigned i = 0; i < length(oeaMatches); ++i)
    {
        if (oeaMatches[i])
            continue;  // Skip, is matched.
        result.push_back(std::make_pair(oeaClusters[i], invalidCCR));
    }

    // Sort resulting clusters.
    std::sort(result.begin(), result.end());
}

// ----------------------------------------------------------------------------
// Class ClusterMatching
// ----------------------------------------------------------------------------

void ClusterMatching::run(std::vector<std::pair<OeaClusterRecord, ClippingClusterRecord> > & matchedClusters,
                          std::vector<OeaClusterRecord> const & oeaClusters,
                          std::vector<ClippingClusterRecord> const & clippingClusters)
{
    impl->run(matchedClusters, oeaClusters, clippingClusters);
}

ClusterMatching::ClusterMatching() : impl(new ClusterMatchingImpl)
{}

ClusterMatching::~ClusterMatching()
{}
