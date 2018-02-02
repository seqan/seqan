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
// Author: Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// Implementation for maximum weighted matching for general graphs.
// The implemented algorithm is the greedy algorithm with a look-ahead.
// That means the weight is maximized among the BLOCKSIZE heaviest edges.
// The performance ratio of this algorithm is 1/2.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_

#include <utility>
#include <algorithm>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function maximumWeightedMatchingGreedy()
// ----------------------------------------------------------------------------

typedef std::vector<std::pair<std::size_t, std::size_t> > TConflictVect;

template <long unsigned BLOCKSIZE, typename TCargo>
inline TCargo _evaluateConflicts(uint32_t & isUsed, std::array<TCargo, BLOCKSIZE> const & weights,
                                 TConflictVect const & conflicts)
{
    // first variant (eliminate second element of pair, which is smaller)
    std::size_t eliminate = conflicts.front().second;
    TConflictVect remainingConflicts;
    std::copy_if(conflicts.begin(), conflicts.end(), std::back_inserter(remainingConflicts),
                 [&eliminate] (std::pair<std::size_t, std::size_t> const & conflict)
    {
        return conflict.first != eliminate && conflict.second != eliminate;
    });
    TCargo excludedWeight1 = weights[eliminate];
    uint32_t isUsed1 = isUsed;
    if (!remainingConflicts.empty())
        excludedWeight1 += _evaluateConflicts(isUsed1, weights, remainingConflicts);

    // second variant (eliminate first element of pair, which is larger)
    eliminate = conflicts.front().first;
    TCargo excludedWeight2 = weights[eliminate];

    // trim traversion if weight2 is too high
    if (excludedWeight1 <= excludedWeight2)
    {
        isUsed = isUsed1 & ~(1 << conflicts.front().second);  // delete bit for 2nd element
        return excludedWeight1;
    }

    bool noDependency = conflicts.size() - remainingConflicts.size() == 1u;
    remainingConflicts.clear();
    std::copy_if(conflicts.begin(), conflicts.end(), std::back_inserter(remainingConflicts),
                 [&eliminate] (std::pair<std::size_t, std::size_t> const & conflict)
    {
        return conflict.first != eliminate && conflict.second != eliminate;
    });

    // trim traversion if the removal of the first conflict does not influence any other conflict
    if (noDependency && conflicts.size() - remainingConflicts.size() == 1u)
    {
        isUsed = isUsed1 & ~(1 << conflicts.front().second);  // delete bit for 2nd element
        return excludedWeight1;
    }

    uint32_t isUsed2 = isUsed;
    if (!remainingConflicts.empty())
        excludedWeight2 += _evaluateConflicts(isUsed2, weights, remainingConflicts);

    // evaluate variants
    if (excludedWeight1 < excludedWeight2)
    {  // use first variant
        isUsed = isUsed1 & ~(1 << conflicts.front().second);  // take isUsed1 and delete bit (eliminate 2nd of pair)
        return excludedWeight1;
    }
    else
    {  // use second variant
        isUsed = isUsed2 & ~(1 << conflicts.front().first);  // take isUsed2 and delete bit (eliminate 1st of pair)
        return excludedWeight2;
    }
}

// Compute greedy MWM (performance ratio 1/2)
// look into BLOCKSIZE edges at once and maximize their weight
template <long unsigned BLOCKSIZE = 1, typename TCargo>
TCargo maximumWeightedMatchingGreedy(Graph<Undirected<TCargo> > const & graph)
{
    typedef Graph<Undirected<TCargo> > TUGraph;
    typedef typename EdgeDescriptor<TUGraph>::Type TEdgeDescr;
    typedef typename Iterator<TUGraph, EdgeIterator>::Type TEdgeIter;
    typedef typename Iterator<TUGraph, AdjacencyIterator>::Type TAdjacIterator;
    typedef typename VertexDescriptor<TUGraph>::Type TVertexDescr;

    // set up edge vector and bit vector for conflicting edges
    std::vector<TEdgeIter> edges;
    std::vector<bool> conflictFree;
    reserve(edges, numEdges(graph));
    resize(conflictFree, numEdges(graph), true);

    for (TEdgeIter edgeIt(graph); !atEnd(edgeIt); goNext(edgeIt))
        edges.push_back(edgeIt);

    // sort edges with respect to their weight, start with the highest
    std::sort(edges.begin(), edges.end(), [] (auto a, auto b) { return getCargo(*a) >= getCargo(*b); });

    TCargo maxWeight{};

    if (BLOCKSIZE == 1)
    {
        for (std::size_t idx = 0u; idx < length(edges); ++idx)
        {
            auto const & edge = *edges[idx];
            if (!conflictFree[edge->data_id])  // skip edge if conflict with a previous edge
                continue;

            maxWeight += getCargo(edge);  // edge is contained in the matching

            // mark all adjacent edges
            TVertexDescr const & src = getSource(edge);
            for (TAdjacIterator ai(graph, src); !atEnd(ai); goNext(ai))
            {
                TEdgeDescr rmEdge = findEdge(graph, src, *ai);
                conflictFree[rmEdge->data_id] = false;
            }

            TVertexDescr const & trg = getTarget(edge);
            for (TAdjacIterator ai(graph, trg); !atEnd(ai); goNext(ai))
            {
                TEdgeDescr rmEdge = findEdge(graph, trg, *ai);
                conflictFree[rmEdge->data_id] = false;
            }

            SEQAN_ASSERT(!conflictFree[edge->data_id]);
        }
    }
    else
    {
        static_assert(BLOCKSIZE <= 32u, "BLOCKSIZE is only supported for values lower or equal 32.");
        uint32_t isUsed;
        std::array<TCargo, BLOCKSIZE> weights;
        std::vector<std::size_t> selection;
        selection.reserve(BLOCKSIZE);
        std::size_t idx = 0u;

        while (idx < length(edges))
        {
            for (selection.clear(); selection.size() < BLOCKSIZE && idx < length(edges); ++idx)
            {
                if (conflictFree[(*edges[idx])->data_id])
                {
                    weights[selection.size()] = getCargo(*edges[idx]);
                    selection.push_back(idx);
                }
            }

            // find conflicts
            isUsed = 0xffffffff;
            TConflictVect conflicts;
            for (unsigned long i = 0u; i < selection.size(); ++i)
            {
                TVertexDescr const & src = getSource(*edges[selection[i]]);
                TVertexDescr const & trg = getTarget(*edges[selection[i]]);
                for (unsigned long j = i + 1u; j < selection.size(); ++j)
                {
                    if (src == getSource(*edges[selection[j]]) || trg == getSource(*edges[selection[j]]) ||
                        src == getTarget(*edges[selection[j]]) || trg == getTarget(*edges[selection[j]]))
                    {
                        conflicts.push_back(std::make_pair(i, j));
                    }
                }
            }

            if (!conflicts.empty())
                _evaluateConflicts<BLOCKSIZE, TCargo>(isUsed, weights, conflicts);

            for (std::size_t i = 0u; i < selection.size(); ++i)
            {
                if (isUsed & (1 << i))  // i-th selection is in the MWM
                {
                    maxWeight += getCargo(*edges[selection[i]]);

                    // mark all adjacent edges
                    TVertexDescr const &src = getSource(*edges[selection[i]]);
                    for (TAdjacIterator ai(graph, src); !atEnd(ai); goNext(ai))
                        conflictFree[findEdge(graph, src, *ai)->data_id] = false;

                    TVertexDescr const &trg = getTarget(*edges[selection[i]]);
                    for (TAdjacIterator ai(graph, trg); !atEnd(ai); goNext(ai))
                        conflictFree[findEdge(graph, trg, *ai)->data_id] = false;
                }
            }
        }
    }
    return maxWeight;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_MAXIMUM_WEIGHTED_MATCHING_H_
