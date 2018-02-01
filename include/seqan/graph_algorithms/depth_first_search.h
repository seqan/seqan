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
// Author: Tobias Raussch <rausch@embl.de>
// Author: Ryan Wick <rrwick@gmail.com>
// ==========================================================================
// Implementation of Depth-First-Search algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_
#define INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_

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
// Function depthFirstSearch()
// ----------------------------------------------------------------------------
/*!
 * @fn depthFirstSearch
 * @headerfile <seqan/graph_algorithms.h>
 * @brief Implements a depth-first search on a graph.
 *
 * @signature void depthFirstSearch(predecessor, discovery, finish, g);
 *
 * @param[out] predecessor A property map.Predecessor subgraph produced by the depth-first search.
 * @param[out] discovery   A property map.The discovery time of a vertex v.
 * @param[out] finish      A property map.The time when v's adjacency list has been fully explored.
 * @param[in]  g           A graph. Types: Undirected Graph, Directed Graph
 *
 * In contrast to a breadth-first search the depth-first search is repeated from multiple sources if the graph is not
 * connected.  Hence, depth-first search produces a depth-first forest.  To ensure each vertex ends up in exactly one
 * tree we need not just a distance but a discovery and finishing time.
 *
 * @section Example
 *
 * @include demos/dox/graph_algorithms/depth_first_search.cpp
 *
 * @include demos/dox/graph_algorithms/depth_first_search.cpp.stdout
 *
 * @see breadthFirstSearch
 */
template <typename TSpec, typename TPredecessorMap, typename TDiscoveryTimeMap, typename TFinishingTimeMap>
void depthFirstSearch(TPredecessorMap & predecessor,
                      TDiscoveryTimeMap & disc,
                      TFinishingTimeMap & finish,
                      Graph<TSpec> const & g)
{
    typedef Graph<TSpec> TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPredecessorMap>::Type TPredVal;
    typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator;

    enum class _DfsTask : uint8_t
    {
        EXPLORE,
        FINISH
    };

    // Initialization - set each vertex as unvisited and with no predecessor.
    resizeVertexMap(predecessor, g);
    resizeVertexMap(disc, g);
    resizeVertexMap(finish, g);
    TPredVal nil = getNil<TVertexDescriptor>();
    String<bool> tokenMap;
    resizeVertexMap(tokenMap, g);
    TVertexIterator it(g);
    for(;!atEnd(it);goNext(it))
    {
        assignProperty(tokenMap, getValue(it), false);
        assignProperty(predecessor, getValue(it), nil);
    }
    TSize time = 0;

    // We do the DFS non-recursively using a stack. The stack holds two possible tasks:
    //   EXPLORE, which means we must follow that vertex's edges
    //   FINISH, which means the vertex just needs a finish time
    std::vector<std::pair<TVertexDescriptor, _DfsTask> > vStack;

    // The graph may not be connected, so start at every vertex.
    goBegin(it);
    for(; !atEnd(it); goNext(it))
    {
        // If the vertex has already been visited, skip it.
        TVertexDescriptor v = getValue(it);
        if (getProperty(tokenMap, v))
        {
            continue;
        }

        vStack.clear();
        vStack.emplace_back(v, _DfsTask::EXPLORE);
        while (!vStack.empty())
        {
            // Get the vertex and task from the top of the stack.
            auto stackItem = vStack.back();
            vStack.pop_back();
            TVertexDescriptor v = stackItem.first;

            if (stackItem.second == _DfsTask::FINISH)
            {
                assignProperty(finish, v, ++time);
            }
            else if (!getProperty(tokenMap, v))       // If the task is EXPLORE and the vertex is not visited...
            {
                assignProperty(tokenMap, v, true);    // label as visited
                assignProperty(disc, v, ++time);      // set discovery time
                vStack.emplace_back(v, _DfsTask::FINISH); // add a task to the stack so the vertex will get a finish time

                // Add EXPLORE tasks to the stack for each unvisited adjacent vertex. They are added in reverse order
                // so the first adjacent vertex will be the first to come off the stack. This is to mimic the behaviour
                // of the original recursive implementation of this function.
                TAdjacencyIterator itad(g, v);
                goEnd(itad);
                while (!atBegin(itad))
                {
                    goPrevious(itad);
                    TVertexDescriptor nextV = getValue(itad);
                    if (!getProperty(tokenMap, nextV))
                    {
                        vStack.emplace_back(nextV, _DfsTask::EXPLORE);
                        assignProperty(predecessor, nextV, v);
                    }
                }
            }
        }
    }
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_ALGORITHMS_DEPTH_FIRST_SEARCH_H_
