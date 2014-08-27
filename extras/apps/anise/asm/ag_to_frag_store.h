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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_AG_TO_FRAG_STORE_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_AG_TO_FRAG_STORE_H_

#include <numeric>

#define DEBUG_INCONSISTENT_LEN 0

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ---------------------------------------------------------------------------
// Class AlignmentGraphToABruijnGraphConverter_
// ---------------------------------------------------------------------------

// Helper class for the conversion.  Each substep corresponds to a method which groups the related functions more
// closely and makes the code easier to read.

// TODO(holtgrew): Fix AlignmentGraph const-correctness.

template <typename TCargo, typename TSpec, typename TStringSet, typename TCargo2, typename TSpec2>
class AlignmentGraphToABruijnGraphConverter_
{
public:
    typedef Graph<Directed<TCargo, TSpec> > TABruijnGraph;
    typedef Graph<Alignment<TStringSet, TCargo2, TSpec2> > TAlignmentGraph;

    typedef typename VertexDescriptor<TAlignmentGraph>::Type TVertexDescriptor;
    TVertexDescriptor NIL;

    AlignmentGraphToABruijnGraphConverter_() : NIL(getNil<TVertexDescriptor>())
    {}

    void run(TABruijnGraph & out,
             StringSet<String<unsigned> > & vertexMap,
             String<unsigned> & vertexLengths,
             TAlignmentGraph const & in)
    {
        // Compute connected components of AG in.
        unsigned numComponents = 0;     // number of components (including empty ones)
        String<unsigned> componentMap;  // vertices of i => component
        // StringSet<String<unsigned> > vertexMap;  // vertex in A-Brujin graph to vertex in AG
        _computeComponents(numComponents, componentMap, vertexMap, in);

        // Build the vertex lengths map.
        _buildVertexLengthsMap(vertexLengths, componentMap, in);

        // From this, we can build a set of directed edges (weighted by the multiplicity).
        std::map<std::pair<unsigned, unsigned>, unsigned> edgeMult;
        _collectEdges(edgeMult, componentMap, in);

        // Build the final A-Brujin graph.
        _buildABruijnGraph(out, numComponents, edgeMult);
    }

    void _buildVertexLengthsMap(String<unsigned> & vertexLengths,
                                String<unsigned> const & componentMap,
                                TAlignmentGraph const & in)
    {
        typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIter;
        resize(vertexLengths, length(componentMap), 0);
        for (TVertexIter it(const_cast<TAlignmentGraph &>(in)); !atEnd(it); goNext(it))
        {
            vertexLengths[componentMap[*it]] = fragmentLength(const_cast<TAlignmentGraph &>(in), *it);
            // std::cerr << "vertexLengths[" << componentMap[*it] << "] == " << vertexLengths[componentMap[*it]] << "\n";
        }
    }

    // Compute connected components and a mapping from the vertices of in to the componet.
    void _computeComponents(unsigned & numComponents,
                            String<unsigned> & componentMap,
                            StringSet<String<unsigned> > & vertexMap,
                            TAlignmentGraph const & in)
    {
        // TODO(holtgrew): Repair const-ness of connectedComponents().
        numComponents = connectedComponents(const_cast<TAlignmentGraph &>(in), componentMap);

        // Build vertexMap.
        resize(vertexMap, numComponents);
        for (unsigned u = 0; u < length(componentMap); ++u)
        {
            // std::cerr << "vertexMap[" << componentMap[u] << "] = " << u << " (" << sequenceId(in, u) << ")\n";
            appendValue(vertexMap[componentMap[u]], u);
        }
    }

    // Compute set of directed edges from this.
    void _collectEdges(std::map<std::pair<unsigned, unsigned>, unsigned> & edgeMult,
                       String<unsigned> const & componentMap,
                       TAlignmentGraph const & ag)
    {
        TStringSet /*const*/ & ss = stringSet(const_cast<TAlignmentGraph &>(ag));
        // We enumerate the sequence edges (u, v) for each sequence in the AG.
        for (unsigned idx = 0; idx < length(ss); ++idx)
        {
            unsigned id = positionToId(ss, idx);
            TVertexDescriptor u = findVertex(const_cast<TAlignmentGraph &>(ag), id, 0);
            while (u != NIL)
            {
                unsigned pos = endPosition(label(const_cast<TAlignmentGraph &>(ag), u));
                TVertexDescriptor v = findVertex(const_cast<TAlignmentGraph &>(ag), id, pos);
                if (v != NIL)
                    edgeMult[std::make_pair(componentMap[u], componentMap[v])] += 1;
                u = v;
            }
        }
        // for (std::set<std::pair<unsigned, unsigned> >::const_iterator it = edges.begin(); it != edges.end(); ++it)
        //     std::cerr << "added edge " << it->first << " -> " << it->second << "\n";
    }

    // Build the final A-Bruijn graph.
    void _buildABruijnGraph(TABruijnGraph & result,
                            unsigned numVertices,
                            std::map<std::pair<unsigned, unsigned>, unsigned> const & edgeMult)
    {
        // Clear everything from graph and insert the vertices.
        clear(result);
        for (unsigned i = 0; i < numVertices; ++i)
            addVertex(result);

        typedef std::map<std::pair<unsigned, unsigned>, unsigned> TMap;
        for (typename TMap::const_iterator it = edgeMult.begin(); it != edgeMult.end(); ++it)
            addEdge(result, it->first.first, it->first.second, it->second);
    }
};

// ---------------------------------------------------------------------------
// Class SmallBulgeRemover_
// ---------------------------------------------------------------------------

// Edge weights are multiplicities, vertexLengths are vertex weights/number of characters.

template <typename TCargo, typename TSpec>
class SmallBulgeRemover_
{
public:
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    // For each edge, we store (vertex_id, vertex length).  Adjacency edge if vertex_id == -1.
    typedef Graph<Directed<std::pair<int, unsigned> > > TEnredoGraph;
    typedef typename EdgeDescriptor<TEnredoGraph>::Type TEnredoEdge;

    unsigned girth;

    explicit
    SmallBulgeRemover_(unsigned girth) : girth(girth)
    {}

    void run(TGraph & outGraph,
             StringSet<String<unsigned> > & outVertexMap,
             TGraph const & inGraph,
             StringSet<String<unsigned> > const & inVertexMap,
             String<unsigned> const & vertexLengths)
    {
        // { std::ofstream af("aba.dot"); write(af, inGraph, DotDrawing()); }

        // Build Enredo graph from inGraph for cycle enumeration.
        TEnredoGraph eg;
        _buildEnredoGraph(eg, inGraph, vertexLengths);

        // { std::ofstream ef("eg.dot"); write(ef, eg, DotDrawing()); }

        // Enumerate all small undirected cycles in Enredo graph.
        StringSet<String<unsigned> > smallCycles;  // vertices in eg
        StringSet<String<TEnredoEdge> > smallCycleEdges;
        _enumerateSmallCycles(smallCycles, smallCycleEdges, eg);
        SEQAN_ASSERT_EQ(length(smallCycles), length(smallCycleEdges));
        // for (unsigned i = 0; i < length(smallCycles); ++i)
        // {
        //     std::cerr << "PATH\t" << i << "\t::";
        //     for (unsigned j = 0; j < length(smallCycles[i]); ++j)
        //         std::cerr << "\t" << smallCycles[i][j] / 2;
        //     std::cerr << "\n";
        // }

        // Contract all cycles to the rightmost vertex (last one in smallCycles[i]).
        TGraph tmpGraph;
        StringSet<String<unsigned> > tmpVertexMap;
        _contractCycles(tmpGraph, tmpVertexMap, inGraph, inVertexMap, eg, smallCycles);

        SEQAN_ASSERT_EQ(lengthSum(tmpVertexMap), lengthSum(inVertexMap));

        // outVertexMap = tmpVertexMap;
        // outGraph = tmpGraph;
        // return;

        clear(outGraph);
        for (unsigned i = 0; i < numVertices(tmpGraph); ++i)
            addVertex(outGraph);
        clear(outVertexMap);
        resize(outVertexMap, numVertices(outGraph));

        // Compute SCCs (directed cycles) and contract them, too.
        String<unsigned> sccs;
        unsigned numSccs = stronglyConnectedComponents(tmpGraph, sccs);
        // for (unsigned i = 0; i < length(sccs); ++i)
        //     std::cerr << "SCC\t" << i << "\t==\t" << sccs[i] << "\n";
        String<unsigned> sccTarget;  // where to contract to
        resize(sccTarget, numSccs, maxValue<unsigned>());

        for (unsigned i = 0; i < numVertices(tmpGraph); ++i)
            if (sccTarget[sccs[i]] == maxValue<unsigned>())
                sccTarget[sccs[i]] = i;
        // for (unsigned i = 0; i < length(sccTarget); ++i)
        //     std::cerr << "SCC TARGET\t" << i << "\t==\t" << sccTarget[i] << "\n";

        std::set<unsigned> done;
        for (typename Iterator<TGraph, EdgeIterator>::Type it(tmpGraph); !atEnd(it); goNext(it))
        {
            unsigned destS = sccTarget[sccs[sourceVertex(it)]];
            unsigned destT = sccTarget[sccs[targetVertex(it)]];
            unsigned s = sourceVertex(it);
            unsigned t = targetVertex(it);
            if (sccs[s] == sccs[t])
            {
                if (!done.count(t))
                    append(outVertexMap[destS], tmpVertexMap[t]);
                if (!done.count(s))
                    append(outVertexMap[destS], tmpVertexMap[s]);
            }
            else
            {
                if (!done.count(s))
                    append(outVertexMap[destS], tmpVertexMap[s]);
                if (!done.count(t))
                    append(outVertexMap[destT], tmpVertexMap[t]);
                addEdge(outGraph, destS, destT, cargo(*it));
            }
            done.insert(s);
            done.insert(t);
        }

        // Copy over for singletons.
        std::set<unsigned> all;
        for (unsigned i = 0; i < numVertices(tmpGraph); ++i)
            all.insert(i);
        std::set<unsigned> singletons;
        std::set_difference(all.begin(), all.end(), done.begin(), done.end(), std::inserter(singletons, singletons.end()));
        for (std::set<unsigned>::const_iterator it = singletons.begin(); it != singletons.end(); ++it)
            outVertexMap[*it] = tmpVertexMap[*it];

        SEQAN_ASSERT_EQ(numVertices(outGraph), numVertices(inGraph));
        SEQAN_ASSERT_EQ(lengthSum(outVertexMap), lengthSum(inVertexMap));

        // { std::ofstream af("contracted.dot"); write(af, outGraph, DotDrawing()); }
    }

    void _buildEnredoGraph(TEnredoGraph & eg, TGraph const & inGraph, String<unsigned> const & vertexLengths)
    {
        // If v is a vertex in inGraph then (2 * v) is the tail vertex in eg and (2 * v + 1) is the head vertex in eg.

        // Add vertices and sequence edges.
        for (unsigned i = 0; i < numVertices(inGraph); ++i)
        {
            addVertex(eg);  // idx: 2 * i
            addVertex(eg);  // idx: 2 * i + 1
            addEdge(eg, 2 * i, 2 * i + 1, std::pair<int, unsigned>(i, vertexLengths[i]));
        }

        // Add edges.
        for (typename Iterator<TGraph, EdgeIterator>::Type it(const_cast<TGraph &>(inGraph)); !atEnd(it); goNext(it))
        {
            unsigned u = 2 * sourceVertex(it) + 1;  // head
            unsigned v = 2 * targetVertex(it);      // tail
            addEdge(eg, u, v, std::pair<int, unsigned>(-1, 0));
        }
    }

    template <typename T>
    std::pair<T, T> _sortedPair(T lhs, T rhs)
    {
        return std::make_pair(std::min(lhs, rhs), std::max(lhs, rhs));
    }

    // Enumerate small cycles in Enredo graph.
    void _enumerateSmallCycles(StringSet<String<unsigned> > & smallCycles,
                               StringSet<String<TEnredoEdge> > & smallCycleEdges,
                               TEnredoGraph const & graph)
    {
        // Compute spanning tree.
        String<unsigned> treeEdges;
        String<unsigned> edgeWeights;
        resize(edgeWeights, numEdges(graph), 1);
        kruskalsAlgorithm(graph, /*ignored*/0, edgeWeights, treeEdges);

        // Build set of edges in the MST.
        std::set<std::pair<unsigned, unsigned> > treeSet;
        for (unsigned i = 0; i < length(treeEdges); i += 2)
        {
            // std::cerr << "INSERTING TREE\t" << treeEdges[i] << "\t" << treeEdges[i + 1] << "\n";
            treeSet.insert(_sortedPair(treeEdges[i], treeEdges[i + 1]));
        }

        // Build undirected graph for iteration.
        typedef Graph<Undirected<TEnredoEdge> > TUG;
        TUG ug;
        for (unsigned i = 0; i < numVertices(graph); ++i)
            addVertex(ug);
        for (typename Iterator<TEnredoGraph, EdgeIterator>::Type itE(const_cast<TEnredoGraph &>(graph)); !atEnd(itE); goNext(itE))
        {
            addEdge(ug, sourceVertex(itE), targetVertex(itE), *itE);
            // if (treeSet.count(_sortedPair(sourceVertex(itE), targetVertex(itE))))
            //     std::cerr << "TREE EDGE    \t" << sourceVertex(itE) << "\t" << targetVertex(itE) << "\n";
            // else
            //     std::cerr << "NON-TREE EDGE\t" << sourceVertex(itE) << "\t" << targetVertex(itE) << "\n";
        }

        // Run DFS on tree for each non-tree edge to find cycle with this edge.
        for (typename Iterator<TEnredoGraph, EdgeIterator>::Type itE(const_cast<TEnredoGraph &>(graph)); !atEnd(itE); goNext(itE))
        {
            if (treeSet.count(_sortedPair(sourceVertex(itE), targetVertex(itE))))
                continue;  // Skip tree edges.
            // std::cerr << "STARTING FROM EDGE " << sourceVertex(itE) << ", " << targetVertex(itE) << "\n";
            unsigned s = sourceVertex(itE), t = targetVertex(itE);

            // We search the unique cycle from s to t using DFS.  Proceed through at most girth characters.
            std::deque<std::pair<unsigned, unsigned> > Q;  // queue (distance, vertex id)
            std::map<unsigned, unsigned> pred;  // pred[v] predecessor of v in iteration
            std::map<unsigned, TEnredoEdge> predEdge;  // edge that v is reached from pred[v]
            Q.push_back(std::make_pair(s, 0));
            pred[s] = s;

            while (!Q.empty())
            {
                unsigned v = Q.front().first, d = Q.front().second;
                // std::cerr << "REACHED\t" << v << "\t" << d << "\n";
                Q.pop_front();

                for (typename Iterator<TUG, OutEdgeIterator>::Type it(ug, v); !atEnd(it); goNext(it))
                {
                    // std::cerr << "CONSIDERING\t" << sourceVertex(it) << "\t" << targetVertex(it) << "\n";
                    if (!treeSet.count(_sortedPair(sourceVertex(it), targetVertex(it))))
                    {
                        // std::cerr << "  SKIPPING NON-TREE\n";
                        continue;  // Skip non-tree edges.
                    }

                    // int labelVertex = cargo(cargo(*it)).first;
                    unsigned labelLen = cargo(cargo(*it)).second;

                    SEQAN_ASSERT_EQ(sourceVertex(it), v);
                    unsigned x = targetVertex(it);
                    // std::cerr << "EDGE\t" << sourceVertex(it) << "\t" << targetVertex(it) << "\t" << labelLen << "\n";
                    if (pred.count(targetVertex(it)))
                        continue;  // Skip, already reached
                    if (d + labelLen > girth)
                        continue;  // Skip, too many chars.
                    pred[x] = v;
                    predEdge[x] = cargo(*it);
                    // std::cerr << "    pred[" << x << "] = " << v << "\n";

                    if (x == t)  // found our cycle
                    {
                        // std::cerr << "FOUND CYCLE FROM s == " << s << " TO t == " << t << "\n";
                        // Write out cycle.
                        String<unsigned> cycle;
                        std::vector<TEnredoEdge> cycleEdges;
                        unsigned z = x;
                        for (; pred[z] != z; z = pred[z])
                        {
                            appendValue(cycle, z);
                            SEQAN_ASSERT(predEdge.count(z));
                            appendValue(cycleEdges, predEdge[z]);
                        }
                        appendValue(cycle, z);
                        appendValue(cycleEdges, *itE);

                        // std::cerr << "FOUND CYCLE\t";
                        // for (unsigned i = 0; i < length(cycle); ++i)
                        //     std::cerr << " " << cycle[i];
                        // std::cerr << "]--\n";
                        appendValue(smallCycles, cycle);
                        appendValue(smallCycleEdges, cycleEdges);

                        // Done.  Break out.
                        Q.clear();
                        break;
                    }

                    // Continue from d.
                    // std::cerr << "PUSH(" << x << ", " << d + labelLen << ")\n";
                    Q.push_back(std::make_pair(x, d + labelLen));
                }
            }
        }
    }

    void _contractCycles(TGraph & outGraph,
                         StringSet<String<unsigned> > & outVertexMap,
                         TGraph const & inGraph,
                         StringSet<String<unsigned> > const & inVertexMap,
                         TEnredoGraph const & eg,
                         StringSet<String<unsigned> > const & smallCycles)
    {
        // -------------------------------------------------------------------
        // For each cycle: find left-/rightmost vertices in cycles
        // -------------------------------------------------------------------

        // Obtain topological sorting of vertices.
        String<unsigned> tmpOrder;
        topologicalSort(const_cast<TEnredoGraph &>(eg), tmpOrder);
        // for (unsigned i = 0; i < length(tmpOrder); ++i)
        //     std::cerr << "tmpOrder[" << i << "]\t=\t" << tmpOrder[i] << "\n";
        std::vector<std::pair<unsigned, unsigned> > order;
        for (unsigned i = 0; i < length(tmpOrder); ++i)
            order.push_back(std::make_pair(tmpOrder[i], i));
        std::sort(order.begin(), order.end());        
        // for (unsigned i = 0; i < length(order2); ++i)
        //     std::cerr << "order2[" << i << "]\t=\t" << order2[i].first << ", " << order2[i].second << "\n";

        // Compute leftmost/rightmost vertices
        std::vector<unsigned> leftmost(length(smallCycles), maxValue<unsigned>());
        std::vector<unsigned> rightmost(length(smallCycles), maxValue<unsigned>());
        for (unsigned i = 0; i < length(smallCycles); ++i)
        {
            for (unsigned j = 0; j < length(smallCycles[i]); ++j)
                if (leftmost[i] == maxValue<unsigned>() || order[smallCycles[i][j]] < order[leftmost[i]])
                    leftmost[i] = smallCycles[i][j];
            for (unsigned j = 0; j < length(smallCycles[i]); ++j)
                if (rightmost[i] == maxValue<unsigned>() || order[smallCycles[i][j]] > order[rightmost[i]])
                    rightmost[i] = smallCycles[i][j];
            // std::cerr << "leftmost[" << i << "] == " << leftmost[i] << "\n";
            // std::cerr << "rightmost[" << i << "] == " << rightmost[i] << "\n";
        }

        // -------------------------------------------------------------------
        // Prepare Contraction (build map)
        // -------------------------------------------------------------------

        // Now contract all vertices of a cycle except for the leftmost one into the rightmost of the cycle.
        UnionFind<int> uf;
        resize(uf, numVertices(eg));
        for (unsigned i = 0; i < length(smallCycles); ++i)
            for (unsigned j = 0; j < length(smallCycles[i]); ++j)
                if (smallCycles[i][j] != leftmost[i])
                    joinSets(uf, findSet(uf, smallCycles[i][j]), findSet(uf, rightmost[i]));

        // Build a partition of the vertices, by UF clusters.
        StringSet<String<unsigned> > partition;
        resize(partition, numVertices(eg));
        for (unsigned u = 0; u < numVertices(eg); ++u)
        {
            // std::cerr << "CLUSTER\t" << u << " => " << findSet(uf, u) << "\n";
            appendValue(partition[findSet(uf, u)], u);
        }

        // Build an array of the rightmost entries for each cluster.
        rightmost.clear();
        rightmost.resize(length(partition), maxValue<unsigned>());
        for (unsigned i = 0; i < length(partition); ++i)
        {
            for (unsigned j = 0; j < length(partition[i]); ++j)
                if (rightmost[i] == maxValue<unsigned>() || order[partition[i][j]] > order[rightmost[i]])
                    rightmost[i] = partition[i][j];
            // std::cerr << "rightmost[" << i << "] == " << rightmost[i] << "\n";
        }

        // Build a map from input vertex id to final id.
        String<unsigned> contractionMap;
        SEQAN_ASSERT_EQ(numVertices(eg) % 2u, 0u);
        resize(contractionMap, numVertices(eg) / 2, maxValue<unsigned>());
        for (unsigned i = 0; i < numVertices(eg); ++i)
            contractionMap[i / 2] = rightmost[findSet(uf, i)] / 2;  // (/ 2) => EG vertex to A-Bruijn vertex

        // -------------------------------------------------------------------
        // Create output from input graph
        // -------------------------------------------------------------------

        // Prepare output graph with vertices.
        clear(outGraph);
        for (unsigned i = 0; i < numVertices(inGraph); ++i)
            addVertex(outGraph);
        clear(outVertexMap);
        resize(outVertexMap, length(inVertexMap));

        // Build edge list of output graph using the contraction mapping (mapping to multiplicity).
        typedef std::map<std::pair<unsigned, unsigned>, unsigned> TMultMap;
        TMultMap edgeMult;
        for (typename Iterator<TGraph, EdgeIterator>::Type it(const_cast<TGraph &>(inGraph)); !atEnd(it); goNext(it))
        {
            std::pair<unsigned, unsigned> key(contractionMap[sourceVertex(it)],
                                              contractionMap[targetVertex(it)]);
            if (key.first != key.second)  // no self-edges
                edgeMult[key] = std::max(edgeMult[key], cargo(*it));
        }

        // Add edges.
        for (TMultMap::const_iterator it = edgeMult.begin(); it != edgeMult.end(); ++it)
            addEdge(outGraph, it->first.first, it->first.second, it->second);

        // Add vertex map.
        for (typename Iterator<TGraph, VertexIterator>::Type it(const_cast<TGraph &>(inGraph)); !atEnd(it); goNext(it))
            append(outVertexMap[contractionMap[*it]], inVertexMap[*it]);
    }
};

template <typename TCargo, typename TSpec>
class VisitingBlockSeparator_
{
public:
    // The A-Bruijn graph type to use.
    typedef Graph<Directed<TCargo, TSpec> > TGraph;
    // Lengths shorter than this thresholds are short.
    unsigned blockLen;
    // Dis-/enable logging.
    bool logging;
    
    explicit
    VisitingBlockSeparator_(unsigned blockLen) : blockLen(blockLen), logging(false)
    {}

    template <typename TAlignmentGraph>
    void run(TGraph & outGraph,
             StringSet<String<unsigned> > & outVertexMap,
             TGraph const & inGraph,
             StringSet<String<unsigned> > const & inVertexMap,
             String<unsigned> vertexLengths,  // used as buffer below, thus copy
             TAlignmentGraph const & ag)
    {
        TGraph tmpGraph = inGraph;
        StringSet<String<unsigned> > tmpVertexMap = inVertexMap;
        String<unsigned> tmpVertexLengths;

        for (unsigned round = 0; true; ++round)
        {
            (void)round;
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "ROUND\t" << round << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
            if (!_runRound(outGraph, outVertexMap, tmpVertexLengths, tmpGraph, tmpVertexMap, vertexLengths, ag))
                return;  // outGraph, outVertexMap untouched, return
#if DEBUG_INCONSISTENT_LEN 
            { std::stringstream ss; ss << "round_" << round << "_in.dot"; std::ofstream of(ss.str().c_str()); write(of, tmpGraph, DotDrawing()); }
            { std::stringstream ss; ss << "round_" << round << "_out.dot"; std::ofstream of(ss.str().c_str()); write(of, outGraph, DotDrawing()); }
#endif  // #if DEBUG_INCONSISTENT_LEN 
            tmpGraph = outGraph;
            tmpVertexMap = outVertexMap;
        }
    }

    // Remove one visiting block.  Removing more than one at once might be faster but also much more complicated with
    // the book keeping.
    //
    // Returns true if a path was found and corrected.
    template <typename TAlignmentGraph>
    bool _runRound(TGraph & outGraph,
                   StringSet<String<unsigned> > & outVertexMap,
                   String<unsigned> & outVertexLengths,
                   TGraph const & inGraph,
                   StringSet<String<unsigned> > const & inVertexMap,
                   String<unsigned> const & origVertexLengths,
                   TAlignmentGraph const & ag)
    {
        // Compute topological sorting of input graph.
        String<unsigned> order;
        topologicalSort(const_cast<TGraph &>(inGraph), order);
        std::map<unsigned, unsigned> orderMap;
        for (unsigned i = 0; i < length(order); ++i)
        {
            orderMap[order[i]] = i;
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "orderMap[" << order[i] << "] = " << i << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        }

#if DEBUG_INCONSISTENT_LEN 
        for (typename Iterator<TAlignmentGraph, VertexIterator>::Type it(const_cast<TAlignmentGraph &>(ag)); !atEnd(it); goNext(it))
        {
            std::cerr << "AG sequence id " << sequenceId(ag, *it) << " => " << *it << "\n";
        }
#endif  // #if DEBUG_INCONSISTENT_LEN 

        // Compute incoming/outgoing vertices for each vertex in graph.
        StringSet<String<unsigned> > incoming, outgoing;
        resize(incoming, numVertices(inGraph));
        resize(outgoing, numVertices(inGraph));
        for (typename Iterator<TGraph, EdgeIterator>::Type it(const_cast<TGraph &>(inGraph)); !atEnd(it); goNext(it))
        {
            appendValue(outgoing[sourceVertex(it)], targetVertex(it));
            appendValue(incoming[targetVertex(it)], sourceVertex(it));
        }
#if DEBUG_INCONSISTENT_LEN 
        for (unsigned i = 0; i < length(incoming); ++i)
        {
            std::cerr << "incoming[" << i << "] == ";
            for (unsigned j = 0; j < length(incoming[i]); ++j)
                std::cerr << "\t" << incoming[i][j];
            std::cerr << "\n";
        }
        for (unsigned i = 0; i < length(outgoing); ++i)
        {
            std::cerr << "outgoing[" << i << "] == ";
            for (unsigned j = 0; j < length(outgoing[i]); ++j)
                std::cerr << "\t" << outgoing[i][j];
            std::cerr << "\n";
        }
#endif  // #if DEBUG_INCONSISTENT_LEN 

        // Traverse inGraph in DFS fashion starting at topological smallest vertices.  During this traversal, we collect
        // the non-branching paths of length <= blockLen;
        StringSet<String<unsigned> > tmpPaths;
        String<bool> visited;
        resize(visited, numVertices(inGraph), false);
        for (unsigned i = 0; i < length(order); ++i)
        {
            unsigned pathLen = 0;
            String<unsigned> path;
            std::vector<unsigned> stack;
            stack.push_back(order[i]);

            while (!stack.empty())
            {
                // Get next in DFS order, make sure that we visit each vertex once and mark current vertex as visited.
                unsigned u = stack.back();
                stack.pop_back();
                if (visited[u])
                    continue;  // Skip, already visited.
                visited[u] = true;

                // In case the in-degree is greater than 1 and we are on a path, the path stops before u.  u is not
                // considered visited however.  If the in-degree is greater than 1 and it is the first vertex on then
                // we continue the path as in all other cases.  This means that the path is only one vertex long in case
                // that the out degree is also greater than 1.
                //
                // If all of this is not the case then we extend the path by u.
                unsigned inDegree = length(incoming[u]);
                if (inDegree > 1u && !empty(path))
                {
                    visited[u] = false;
                }
                else
                {
                    appendValue(path, u);
                    pathLen += origVertexLengths[u];
                }

                // Enqueue next vertices to visit.
                unsigned oldSize = stack.size();
                for (typename Iterator<TGraph, OutEdgeIterator>::Type it(inGraph, u); !atEnd(it); goNext(it))
                    if (!visited[targetVertex(it)])
                        stack.push_back(targetVertex(it));

                // If the out-degree is not 1 or all outgoing vertices have already been visited then we are in a
                // cul-de-sac or at a branching block.  In this case, we start a new path.  The old path is stored in
                // case it is not too long.
                unsigned outDegree = length(outgoing[u]);
                if (outDegree != 1u || oldSize == stack.size())
                {
                    // In case of being at the end of a non-branching block, add path to tmpPaths if not too long.
                    if (pathLen > 0 && pathLen <= blockLen)
                    {
                        // Do not add singleton vertices that have no correspondence to vertices in the AG any more.
                        // Such vertices are created in a previous contraction step.
                        if (length(path) > 1u || !empty(inVertexMap[path[0]]))
                            appendValue(tmpPaths, path);
                        // std::cerr << "pathLen == " << pathLen << "\n";
                    }
                    clear(path);
                    pathLen = 0;
                }
            }
        }

        // Remove paths where the first vertex in degree <= 1 and the last vertex out degree also >= 1.
        StringSet<String<unsigned> > paths;
        for (unsigned i = 0; i < length(tmpPaths); ++i)
        {
            if (length(incoming[front(tmpPaths[i])]) <= 1 || length(outgoing[back(tmpPaths[i])]) <= 1)
                continue;
            appendValue(paths, tmpPaths[i]);
        }

        if (DEBUG_INCONSISTENT_LEN || logging)
        {
            std::cerr << "PATHS\n";
            for (unsigned i = 0; i < length(paths); ++i)
            {
                std::cerr << i << "\t::\t";
                for (unsigned j = 0; j < length(paths[i]); ++j)
                    std::cerr << "\t" << paths[i][j];
                std::cerr << "\n";
            }
        }

        // For the first such short path P that is not a cul-de-sac, compute extended paths (u, P, v) where u and v are
        // vertices in the graph.  Compute extended path with maximal multiplicity.  If no such extended path exists,
        // use edge (u, P) or (P, v) with highest multiplicity (there must be such an edge).
        unsigned i = maxValue<unsigned>();

        unsigned firstVertex = 0;
        unsigned lastVertex = 0;

        for (unsigned k = 0; k < length(paths); ++k)
        {
            // Find leftmost/rightmost in path.
            SEQAN_ASSERT(orderMap.count(paths[k][0]));
            unsigned orderMin = orderMap[paths[k][0]];
            unsigned orderMax = orderMin;
            unsigned idxMin = 0, idxMax = 0;
            for (unsigned j = 1; j < length(paths[k]); ++j)
            {
                SEQAN_ASSERT(orderMap.count(paths[k][j]));
                unsigned orderJ = orderMap[paths[k][j]];
                if (orderJ < orderMin)
                {
                    idxMin = j;
                    orderMin = orderJ;
                }
                if (orderJ > orderMax)
                {
                    idxMax = j;
                    orderMax = orderJ;
                }
            }

            firstVertex = paths[k][idxMin];
            lastVertex = paths[k][idxMax];
            if (!empty(incoming[firstVertex]) && !empty(outgoing[lastVertex]))
            {
                i = k;
#if DEBUG_INCONSISTENT_LEN 
                std::cerr << "#incoming[" << firstVertex << "] == " << length(incoming[firstVertex]) << "\n";
                std::cerr << "#outgoing[" << lastVertex << "] == " << length(outgoing[lastVertex]) << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                break;
            }
        }

        if (i == maxValue<unsigned>())  // no such path found
        {
            outGraph = inGraph;
            outVertexMap = inVertexMap;
            outVertexLengths = origVertexLengths;
            return false;
        }

        // Build set of path vertices.
        std::set<unsigned> pathVertices;
        std::copy(begin(paths[i], Standard()), end(paths[i], Standard()), std::inserter(pathVertices, pathVertices.end()));

        // Build list of fragments that span over the path to the left/right.
        //
        // We collect a map from sequence id to the vertex that is left/right of the path.
        std::map<unsigned, unsigned> spanningL, spanningR;
#if DEBUG_INCONSISTENT_LEN 
        std::cerr << "firstVertex == " << firstVertex << ", lastVertex == " << lastVertex << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        for (unsigned i = 0; i < length(incoming[firstVertex]); ++i)
        {
            unsigned u = incoming[firstVertex][i];  // vertex u in inGraph
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "incoming[" << firstVertex << "][" << i << "] == " << u << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
            for (unsigned j = 0; j < length(inVertexMap[u]); ++j)
            {
                unsigned v = inVertexMap[u][j];  // vertex v in alignment graph ag, retrieve segments
#if DEBUG_INCONSISTENT_LEN 
                std::cerr << "v == " << v << " == inVertexMap[" << u << "][" << j << "]\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                // SEQAN_ASSERT_NOT(spanningL.count(sequenceId(ag, v)));
                // If the same sequence occurs in two vertices then take the topographically greatest one.
                if (!spanningL.count(sequenceId(ag, v)) || orderMap[u] > orderMap[spanningL[sequenceId(ag, v)]])
                {
                    spanningL[sequenceId(ag, v)] = u;
#if DEBUG_INCONSISTENT_LEN 
                    std::cerr << "  ---> orderMap[" << u << "] == " << orderMap[u] << ", orderMap[" << spanningL[sequenceId(ag, v)] << "] == " << orderMap[spanningL[sequenceId(ag, v)]] << "\n";
                    std::cerr << "spanningL[sequenceId(ag, " << v << ") == " << sequenceId(ag, v) << "] = " << spanningL[sequenceId(ag, v)] << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                }
            }
        }
        for (unsigned i = 0; i < length(outgoing[lastVertex]); ++i)
        {
            unsigned u = outgoing[lastVertex][i];  // vertex u in inGraph
            // std::cerr << "outgoing[" << lastVertex << "][" << i << " == " << u << "\n";
            for (unsigned j = 0; j < length(inVertexMap[u]); ++j)
            {
                unsigned v = inVertexMap[u][j];  // vertex v in alignment graph ag, retrieve segments
                // SEQAN_ASSERT_NOT(spanningR.count(sequenceId(ag, v)));
                // spanningR[sequenceId(ag, v)] = u;
#if DEBUG_INCONSISTENT_LEN 
                std::cerr << "spanningR[sequenceId(ag, " << v << ") == " << sequenceId(ag, v) << "] = " << spanningR[sequenceId(ag, v)] << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                // If the same sequence occurs in two vertices then take the topographically smallest one.
                if (!spanningR.count(sequenceId(ag, v)) || orderMap[u] < orderMap[spanningR[sequenceId(ag, v)]])
                {
                    spanningR[sequenceId(ag, v)] = u;
#if DEBUG_INCONSISTENT_LEN 
                    std::cerr << "  ---> orderMap[" << u << "] == " << orderMap[u] << ", orderMap[" << spanningR[sequenceId(ag, v)] << "] == " << orderMap[spanningR[sequenceId(ag, v)]] << "\n";
                    std::cerr << "spanningR[sequenceId(ag, " << v << ") == " << sequenceId(ag, v) << "] = " << spanningR[sequenceId(ag, v)] << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                }
            }
        }

        // Build a list of spanning fragments and get thickest edge (u, v) at the same time.
        std::map<std::pair<unsigned, unsigned>, unsigned> mult;
        unsigned maxMult = 0;
        std::pair<unsigned, unsigned> maxKey(maxValue<unsigned>(), maxValue<unsigned>());
        for (std::map<unsigned, unsigned>::const_iterator it = spanningL.begin(); it != spanningL.end(); ++it)
            if (spanningR.count(it->first))
            {
                std::pair<unsigned, unsigned> key(it->second, spanningR[it->first]);
                mult[key] += 1;
                if (maxKey.first == maxValue<unsigned>() || mult[key] > maxMult)
                {
                    maxKey = key;
                    maxMult = mult[key];
                }
            }
        // The keys of mult are now the existing spanning edges and maxKey has the highest-multiplicity edge (or is
        // maxValue<unsigned>().
        if (DEBUG_INCONSISTENT_LEN)
        {
            for (std::map<std::pair<unsigned, unsigned>, unsigned>::const_iterator it = mult.begin(); it != mult.end(); ++it)
                std::cerr << "SPANNING BOTH\t" << it->first.first << ", " << it->first.second << "\n";
            std::cerr << "MAX KEY\t" << maxKey.first << ", " << maxKey.second << "\n";
        }

        // Create a set of ids with sequences that have any spanning edge at all (left or right or both).
        std::set<unsigned> anySpanning;
        for (std::map<unsigned, unsigned>::const_iterator it = spanningL.begin(); it != spanningL.end(); ++it)
            anySpanning.insert(it->first);
        for (std::map<unsigned, unsigned>::const_iterator it = spanningR.begin(); it != spanningR.end(); ++it)
            anySpanning.insert(it->first);

        // If we have no spanning fragment then take a fragment spanning from left or right as maximum multiplicity
        // extended path.
        if (maxMult == 0)
        {
            std::map<unsigned, unsigned> multL;
            unsigned maxL = 0, maxLIdx = maxValue<unsigned>();
            for (std::map<unsigned, unsigned>::const_iterator it = spanningL.begin(); it != spanningL.end(); ++it)
            {
                multL[it->second] += 1;
                if (maxLIdx == maxValue<unsigned>() || multL[it->second] > maxL)
                {
                    maxLIdx = it->second;
                    maxL = multL[it->second];
                }
            }
            std::map<unsigned, unsigned> multR;
            unsigned maxR = 0, maxRIdx = maxValue<unsigned>();
            for (std::map<unsigned, unsigned>::const_iterator it = spanningR.begin(); it != spanningR.end(); ++it)
            {
                multR[it->second] += 1;
                if (maxRIdx == maxValue<unsigned>() || multR[it->second] > maxR)
                {
                    maxRIdx = it->second;
                    maxR = multR[it->second];
                }
            }
            SEQAN_CHECK((maxR > 0u || maxL > 0u), "Graph model wrong/bug.");
            if (maxR > maxL)
                maxKey = std::make_pair(maxValue<unsigned>(), maxRIdx);
            else
                maxKey = std::make_pair(maxLIdx, maxValue<unsigned>());
        }

#if DEBUG_INCONSISTENT_LEN 
        std::cerr << "MAX SPANNING\t" << maxKey.first << "\t" << maxKey.second << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        
        // ---------------------------------------------------------------
        // Create output graph with new structure where the path is split.
        // ---------------------------------------------------------------

        // We will create a new graph.  The graph will be based on inGraph with the following modifications:
        //
        // (1) Remove all but the incoming/outgoing edge(s) identified by maxKey, also remove path edges.
        // (2) Add edges (u, v) for all keys of mult except the heaviest one.
        // (3) The alignment graph vertices of the path will be reassigned for all but the sequences on the edge(s)
        //     identified by maxKey and the sequences that have no value in spanningL or spanningR.  Sequences with a
        //     value in spanningL will be assigned to the the value's vertex.  Likewise for sequences with a value in
        //     spanningR.

        // Add vertices to output graph.
        clear(outGraph);
        for (unsigned i = 0; i < numVertices(inGraph); ++i)
            addVertex(outGraph);

        // Build set of edges that we will not add back any more (we will skip ignoring for highest-scoring extended
        // path below).
        std::set<std::pair<unsigned, unsigned> > ignoredEdges;
        for (unsigned i = 0; i < length(incoming[firstVertex]); ++i)
        {
            ignoredEdges.insert(std::make_pair(incoming[firstVertex][i], firstVertex));
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "ignoring " << incoming[firstVertex][i] << ", " << firstVertex << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        }
        for (unsigned i = 0; i < length(outgoing[lastVertex]); ++i)
        {
            ignoredEdges.insert(std::make_pair(lastVertex, outgoing[lastVertex][i]));
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "ignoring " << lastVertex << ", " << outgoing[lastVertex][i] << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        }

        // Copy over the edges from the old graph, ignoring the to-be-removed ones.
        for (typename Iterator<TGraph const, EdgeIterator>::Type it(inGraph); !atEnd(it); goNext(it))
        {
            bool ignore = ignoredEdges.count(std::make_pair(sourceVertex(it), targetVertex(it)));
            // std::cerr << "  in ignored edges(" << sourceVertex(it) << ", " << targetVertex(it) << ") ==> " << ignore << "\n";
            if (maxKey.first != maxValue<unsigned>())
                if (sourceVertex(it) == maxKey.first && targetVertex(it) == firstVertex)
                {
#if DEBUG_INCONSISTENT_LEN 
                    std::cerr << "    highest weight incoming one.\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                    ignore = false;  // Do not ignore the highest weight incoming one.
                }
            if (maxKey.second != maxValue<unsigned>())
                if (sourceVertex(it) == lastVertex && targetVertex(it) == maxKey.second)
                {
#if DEBUG_INCONSISTENT_LEN 
                    std::cerr << "    highest weight outgoing one.\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
                    ignore = false;  // Do not ignore the hight weight outgoing one.
                }
            if (ignore)
            {
                // std::cerr << "IGNORING EDGE\t" << sourceVertex(it) << "\t" << targetVertex(it) << "\n";
                continue;  // Skip if to be ignored.
            }
            addEdge(outGraph, sourceVertex(it), targetVertex(it), cargo(*it));
#if DEBUG_INCONSISTENT_LEN 
            std::cerr << "ADDING EDGE\t" << sourceVertex(it) << "\t" << targetVertex(it) << "\n";
#endif  // #if DEBUG_INCONSISTENT_LEN 
        }

        // Add new edges to graph.
        for (std::map<std::pair<unsigned, unsigned>, unsigned>::const_iterator it = mult.begin(); it != mult.end(); ++it)
            if (it->first != maxKey)
                addEdge(outGraph, it->first.first, it->first.second);

#if DEBUG_INCONSISTENT_LEN 
        std::cerr << "IN VERTEX MAP\n";
        for (unsigned i = 0; i < length(inVertexMap); ++i)
        {
            std::cerr << i;
            for (unsigned j = 0; j < length(inVertexMap[i]); ++j)
            {
                std::cerr << "\t" << sequenceId(ag, inVertexMap[i][j]) << " : " << fragmentBegin(ag, inVertexMap[i][j])
                          << ", " << fragmentLength(ag, inVertexMap[i][j]);
            }
            std::cerr << "\n";
        }
#endif  // #if DEBUG_INCONSISTENT_LEN 

        // Build new vertex map.
        clear(outVertexMap);
        resize(outVertexMap, length(inVertexMap));
        // First, copy over the mapping for non-path vertices.
        for (unsigned i = 0; i < length(inVertexMap); ++i)
            if (!pathVertices.count(i))
                outVertexMap[i] = inVertexMap[i];
        // Then, reassign the mappings for the path vertices if necessary.
        for (std::set<unsigned>::const_iterator it = pathVertices.begin(); it != pathVertices.end(); ++it)
        {
            unsigned gV = *it;  // vertex in inGraph
            for (unsigned i = 0; i < length(inVertexMap[gV]); ++i)
            {
                unsigned agV = inVertexMap[gV][i];  // in alignment graph
                unsigned seqID = sequenceId(ag, agV);  // ID of the sequence's fragment
                if (DEBUG_INCONSISTENT_LEN)
                    std::cerr << "seqID=" << seqID << "\tspanningL.count(seqID)=" << spanningL.count(seqID) << "\tspanningR.count(seqID)=" << spanningR.count(seqID) << "\n";
                if (spanningL.count(seqID) && spanningL[seqID] != maxKey.first)
                {
                    // contract left
                    if (DEBUG_INCONSISTENT_LEN)
                        std::cerr << "CONTRACTING LEFT\t" << agV << "\tto\t" << spanningL[seqID] << "\n";
                    appendValue(outVertexMap[spanningL[seqID]], agV);
                }
                else if (spanningR.count(seqID) && spanningR[seqID] != maxKey.second)
                {
                    // contract right
                    if (DEBUG_INCONSISTENT_LEN)
                        std::cerr << "CONTRACTING RIGHT\t" << agV << "\tto\t" << spanningR[seqID] << "\n";
                    appendValue(outVertexMap[spanningR[seqID]], agV);
                }
                else
                {
                    // keep here
                    if (DEBUG_INCONSISTENT_LEN)
                        std::cerr << "KEEPING\t" << agV << "\tat\t" << gV << "\n";
                    appendValue(outVertexMap[gV], agV);
                }
            }
        }

        // Rebuild vertexLengths.
        clear(outVertexLengths);
        resize(outVertexLengths, length(origVertexLengths), 0u);
        for (unsigned i = 0; i < length(outVertexMap); ++i)
        {
            std::map<unsigned, unsigned> lengths;
            for (unsigned j = 0; j < length(outVertexMap[i]); ++j)
                lengths[sequenceId(ag, outVertexMap[i][j])] += fragmentLength(ag, outVertexMap[i][j]);
            unsigned maxLen = 0;
            for (std::map<unsigned, unsigned>::const_iterator it = lengths.begin(); it != lengths.end(); ++it)
                maxLen = std::max(maxLen, it->second);
            outVertexLengths[i] = maxLen;
        }

        return true;
    }
};

// ---------------------------------------------------------------------------
// Class NonBranchingPathMerger_
// ---------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
class NonBranchingPathMerger_
{
public:
    // The A-Bruijn graph type to use.
    typedef Graph<Directed<TCargo, TSpec> > TGraph;

    void run(TGraph & outGraph,
             StringSet<String<unsigned> > & outVertexMap,
             TGraph const & inGraph,
             StringSet<String<unsigned> > const & inVertexMap)
    {
        // Build in/out degree map for each vertex.
        String<unsigned> inDegree, outDegree;
        resizeVertexMap(inGraph, inDegree, 0);
        resizeVertexMap(inGraph, outDegree, 0);
        for (typename Iterator<TGraph const, EdgeIterator>::Type it(inGraph); !atEnd(it); goNext(it))
        {
            inDegree[targetVertex(it)] += 1;
            outDegree[sourceVertex(it)] += 1;
        }

        // Join vertices along edges if outDegree(tail) <= 1 and inDegree(head) <= 1.
        UnionFind<unsigned> uf;
        resizeVertexMap(inGraph, uf);
        for (typename Iterator<TGraph const, EdgeIterator>::Type it(inGraph); !atEnd(it); goNext(it))
            if (outDegree[sourceVertex(it)] <= 1 && inDegree[targetVertex(it)] <= 1)
                joinSets(uf, findSet(uf, sourceVertex(it)), findSet(uf, targetVertex(it)));

        // Map component representants in uf to vertex id in new graph.
        unsigned n = 0;
        std::map<unsigned, unsigned> componentToVertex;
        for (unsigned i = 0; i < numVertices(inGraph); ++i)
            if (!componentToVertex.count(findSet(uf, i)) && !empty(inVertexMap[i]))
                componentToVertex[findSet(uf, i)] = n++;

        // Create new graph.
        clear(outGraph);
        for (unsigned i = 0; i < n; ++i)
            addVertex(outGraph);

        // Build edge list to add with cargo.
        std::map<std::pair<unsigned, unsigned>, unsigned> edgeMult;
        for (typename Iterator<TGraph const, EdgeIterator>::Type it(inGraph); !atEnd(it); goNext(it))
        {
            std::pair<unsigned, unsigned> key(componentToVertex[findSet(uf, sourceVertex(it))],
                                              componentToVertex[findSet(uf, targetVertex(it))]);
            if (key.first != key.second)
                edgeMult[key] = std::max(edgeMult[key], cargo(*it));
        }

        // Add edges to output graph.
        for (std::map<std::pair<unsigned, unsigned>, unsigned>::const_iterator it = edgeMult.begin(); it != edgeMult.end(); ++it)
            addEdge(outGraph, it->first.first, it->first.second, it->second);

        // Build output vertex map.
        clear(outVertexMap);
        resize(outVertexMap, numVertices(outGraph));
        for (typename Iterator<TGraph const, VertexIterator>::Type it(inGraph); !atEnd(it); goNext(it))
            append(outVertexMap[componentToVertex[findSet(uf, *it)]], inVertexMap[*it]);
    }
};

// ----------------------------------------------------------------------------
// Class MyStoreConfig_
// ----------------------------------------------------------------------------

// Configuration for the FragmentStore we use for visualizing the multi-read alignments of the contigs.

struct MyStoreConfig_
{
	typedef seqan::String<seqan::Dna5Q>	TReadSeq;
	typedef seqan::String<seqan::Dna5Q>	TContigSeq;
	
	typedef double			TMean;
	typedef double			TStd;
	typedef signed char		TMappingQuality;
		
	typedef void					TReadStoreElementSpec;
	typedef seqan::Owner<>                 TReadSeqStoreSpec;
	typedef void					TMatePairStoreElementSpec;
	typedef void					TLibraryStoreElementSpec;
	typedef void					TContigStoreElementSpec;
	typedef void					TContigFileSpec;
	typedef void					TAlignedReadStoreElementSpec;
	typedef seqan::Owner<seqan::ConcatDirect<> >	TAlignedReadTagStoreSpec;
	typedef void					TAnnotationStoreElementSpec;
    
    typedef seqan::Alloc<>					TReadNameSpec;
	typedef seqan::Owner<>	TReadNameStoreSpec;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Function alignmentGraphToABruijnGraph()
// ---------------------------------------------------------------------------

template <typename TCargo, typename TSpec, typename TStringSet, typename TCargo2, typename TSpec2>
void alignmentGraphToABruijnGraph(Graph<Directed<TCargo, TSpec> > & out,
                                  StringSet<String<unsigned> > & vertexMap,
                                  String<unsigned> & vertexLengths,
                                  Graph<Alignment<TStringSet, TCargo2, TSpec2> > const & in)
{
    AlignmentGraphToABruijnGraphConverter_<TCargo, TSpec, TStringSet, TCargo2, TSpec2> converter;
    converter.run(out, vertexMap, vertexLengths, in);
}

// ---------------------------------------------------------------------------
// Function removeSmallBulges()
// ---------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
void removeSmallBulges(Graph<Directed<TCargo, TSpec> > & outGraph,
                       StringSet<String<unsigned> > & outVertexMap,
                       Graph<Directed<TCargo, TSpec> > const & inGraph,
                       StringSet<String<unsigned> > const & vertexMap,
                       String<unsigned> const & vertexLengths,
                       unsigned girth)
{
    SmallBulgeRemover_<TCargo, TSpec> bulgeRemover(girth);
    bulgeRemover.run(outGraph, outVertexMap, inGraph, vertexMap, vertexLengths);
}

// ---------------------------------------------------------------------------
// Function separateVisitingBlocks()
// ---------------------------------------------------------------------------

template <typename TCargo, typename TSpec, typename TAlignmentGraph>
void separateVisitingBlocks(Graph<Directed<TCargo, TSpec> > & outGraph,
                            StringSet<String<unsigned> > & outVertexMap,
                            Graph<Directed<TCargo, TSpec> > const & inGraph,
                            StringSet<String<unsigned> > const & vertexMap,
                            String<unsigned> const & vertexLengths,
                            TAlignmentGraph const & ag,
                            unsigned len)
{
    VisitingBlockSeparator_<TCargo, TSpec> blockSeparator(len);
    blockSeparator.run(outGraph, outVertexMap, inGraph, vertexMap, vertexLengths, ag);
}

// ---------------------------------------------------------------------------
// Function mergeNonBranchingPaths()
// ---------------------------------------------------------------------------

template <typename TCargo, typename TSpec>
void mergeNonBranchingPaths(Graph<Directed<TCargo, TSpec> > & outGraph,
                            StringSet<String<unsigned> > & outVertexMap,
                            Graph<Directed<TCargo, TSpec> > const & inGraph,
                            StringSet<String<unsigned> > const & vertexMap)
{
    NonBranchingPathMerger_<TCargo, TSpec> pathMerger;
    pathMerger.run(outGraph, outVertexMap, inGraph, vertexMap);
}

// ---------------------------------------------------------------------------
// Function removeEdgesBetweenPartitionEntries()
// ---------------------------------------------------------------------------

// Remove alignment edges from ag that lie between two partition entries (vertexPartition is a partition of alignment
// graph vertices).

template <typename TStringSet, typename TScore>
void removeEdgesBetweenPartitionEntries(Graph<Alignment<TStringSet, TScore> > & ag,
                                        StringSet<String<unsigned> > const & vertexPartition)
{
    typedef Graph<Alignment<TStringSet, TScore> > TAG;

    SEQAN_ASSERT_EQ(lengthSum(vertexPartition), numVertices(ag));
    String<unsigned> component;  // vertex -> component
    resize(component, numVertices(ag));
    for (unsigned i = 0; i < length(vertexPartition); ++i)
        for (unsigned j = 0; j < length(vertexPartition[i]); ++j)
            component[vertexPartition[i][j]] = i;

    std::set<std::pair<unsigned, unsigned> > toRemove;
    for (typename Iterator<TAG, EdgeIterator>::Type it(ag); !atEnd(it); goNext(it))
        if (component[sourceVertex(it)] != component[targetVertex(it)])
            toRemove.insert(std::make_pair(sourceVertex(it), targetVertex(it)));

    // std::cerr << "REMOVING\t" << toRemove.size() << " / " << numEdges(ag) << "\n";

    for (std::set<std::pair<unsigned, unsigned> >::const_iterator it = toRemove.begin(); it != toRemove.end(); ++it)
        removeEdge(ag, it->first, it->second);
}

// ---------------------------------------------------------------------------
// Function alignmentGraphToSmoothFragmentStore()
// ---------------------------------------------------------------------------

// Note that this function relies on the "all mated, adjacent reads" assumption.

template <typename TFragmentStore, typename TSequence, typename TCargo, typename TSetSpec, typename TSpec>
bool alignmentGraphToFragmentStore(TFragmentStore & store,
                                   seqan::Graph<seqan::Alignment<seqan::StringSet<TSequence, TSetSpec>, TCargo, TSpec> > const & g,
                                   seqan::Graph<seqan::Undirected<double> > const & distances,
                                   seqan::String<unsigned> const & component,
                                   seqan::String<unsigned> const & order,
                                   unsigned numComponents,
                                   bool logging)
{
    // std::cerr << ">>>>>>>>>>>>\n<<<<<<<<<<<<<<<<\n";
    // NOTE: seqToCluster is indexed by POSITION in the read set of g and not by the ID.

    // TODO(holtgrew): This function is very similar to the updateStoreFromAlignmentGraph computeProfiles functions. Maybe we can share the commonality?
    using namespace seqan;

    typedef Graph<Alignment<StringSet<TSequence, TSetSpec>, TCargo, TSpec> > TAlignmentGraph;
    // typedef Dna5 TAlphabet;

    // Allocate information for which sequence supports the profile at which position.
    resize(store.readSeqStore, length(stringSet(g)));
    resize(store.readStore, length(stringSet(g)));
    resize(store.alignedReadStore, length(stringSet(g)));

    // We can fill the mate pair store here since we know that we have an even number of reads in the read set and the
    // read with a given id is part of the pair (id / 2) and is the (id % 2)-th read in the pair.
    SEQAN_ASSERT_EQ(length(store.readStore) % 2, 0u);
    resize(store.matePairStore, length(store.readStore) / 2);
    for (unsigned i = 0; i < length(store.matePairStore); ++i)
    {
        store.matePairStore[i].readId[0] = i * 2;
        store.matePairStore[i].readId[1] = i * 2 + 1;
    }

    // -----------------------------------------------------------------------
    // Get connected components of distances / read alignment clusters.
    // -----------------------------------------------------------------------

    // Each cluster corresponds to a contig.

    // A cluster is a CC in the graph where each sequences is a vertex and two vertices are connected if they have an
    // overlap alignment.
    String<unsigned> seqToCluster;
    if (logging)
        std::cerr << "# vertices: " << numVertices(distances) << "\n"
                  << "# edges: " << numEdges(distances) << "\n";
    unsigned numClusters = connectedComponents(distances, seqToCluster);
    if (logging)
        std::cerr << "# clusters: " << numClusters << std::endl
                  << "# components: " << numComponents << std::endl;
    resize(store.contigStore, numClusters);
    String<unsigned> contigLengths;
    resize(contigLengths, numClusters, 0);

    for (unsigned i = 0; i < numClusters; ++i)
    {
        std::stringstream ss;
        ss << "contig_" << i;
        appendValue(store.contigNameStore, ss.str());
    }

    // -----------------------------------------------------------------------
    // Visit components in topological order and generate profile sequences.
    // -----------------------------------------------------------------------

    // Get mapping from component to vertices.
    String<String<unsigned> > componentVertices;
    resize(componentVertices, numComponents);
    typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV))
        appendValue(componentVertices[getProperty(component, *itV)], *itV);

    // For each cluster, the currently overlapping reads.
    std::vector<std::set<unsigned> > activeReads(numClusters);

    std::vector<unsigned> gapCount(length(stringSet(g)), 0);

    // Iterate vertices in topological order.
    for (typename Iterator<String<unsigned> const, Rooted>::Type it = begin(order, Rooted()); !atEnd(it); goNext(it))
    {
        unsigned c = *it;     // Current component.
        unsigned fLen = fragmentLength(g, front(componentVertices[c]));
        for (unsigned i = 1; i < length(componentVertices[c]); ++i)
            SEQAN_ASSERT_EQ(fragmentLength(g, front(componentVertices[c][0])),
                            fragmentLength(g, front(componentVertices[c][i])));
        unsigned cl = seqToCluster[idToPosition(stringSet(g), sequenceId(g, front(componentVertices[c])))];  // Current cluster/contig.

        // Update contig lengths.
        unsigned from = contigLengths[cl];
        contigLengths[cl] += fLen;
        if (DEBUG_INCONSISTENT_LEN)
            std::cerr << "==== c == " << c << "\n" << "     from == " << from << "\n";

        // The currently active reads that we see in this round.  Required for inserting gaps below.
        std::set<unsigned> seen;
        std::set<unsigned> done;

        // Insert gaps.
        typedef typename Iterator<String<unsigned>, Rooted>::Type TDescIt;
        for (TDescIt itV = begin(componentVertices[c], Rooted()); !atEnd(itV); goNext(itV))
        {
            unsigned idx = idToPosition(stringSet(g), sequenceId(g, *itV));
            seen.insert(idx);
            unsigned fBeg = fragmentBegin(g, *itV);
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "           fBeg of " << *itV << " is " << fBeg << " (idx == " << idx << ")\n";

            // Register sequence as supporting in profile cl starting at position from in profile.
            if (fBeg == 0u)
            {
                store.readSeqStore[idx] = getValueById(stringSet(g), sequenceId(g, *itV));
                store.readStore[idx].matePairId = idx / 2;  // TODO(holtgrew): Strong assumption :|
                SEQAN_ASSERT_NOT(empty(store.readSeqStore[idx]));
                activeReads[cl].insert(idx);
                store.alignedReadStore[idx].id = idx;
                store.alignedReadStore[idx].readId = idx;
                store.alignedReadStore[idx].contigId = cl;
                store.alignedReadStore[idx].beginPos = from;
                store.alignedReadStore[idx].endPos = from;
                store.alignedReadStore[idx].pairMatchId = idx / 2;
                if (DEBUG_INCONSISTENT_LEN)
                {
                    std::cerr << "store.alignedReadStore[" << idx << "].beginPos == " << from << " | *itV == " << *itV << "\n";
                    std::cerr << "store.alignedReadStore[" << idx << "].endPos (= endPos) == " << from << "| *itV == " << *itV << "\n";
                }
            }
            store.alignedReadStore[idx].endPos = from + fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "store.alignedReadStore[" << idx << "].endPos = " << from << " + " << fLen << " == " << (from + fLen) << "| *itV == " << *itV << "\n";

            unsigned fEnd = fBeg + fLen;

            if (fEnd == length(stringSet(g)[idx]))
                done.insert(idx);
        }

        // Get not seen reads.
        typedef std::set<unsigned>::iterator TSetIt;
        std::set<unsigned> notSeen;
        for (TSetIt it = activeReads[cl].begin(); it != activeReads[cl].end(); ++it)
            notSeen.insert(*it);
        for (TSetIt it = seen.begin(); it != seen.end(); ++it)
            notSeen.erase(*it);
        // Insert gaps into these reads.
        for (TSetIt itS = notSeen.begin(); itS != notSeen.end(); ++itS)
        {
            typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
            typedef typename Value<TAlignedReadStore>::Type TAlignedRead;
            typedef typename TAlignedRead::TGapAnchors TGapAnchors;
            typedef typename TFragmentStore::TReadSeq TReadSeq;
            SEQAN_ASSERT_NOT(empty(store.readSeqStore[*itS]));
            Gaps<TReadSeq, AnchorGaps<TGapAnchors> > gaps(store.readSeqStore[*itS], store.alignedReadStore[*itS].gaps);
            insertGaps(gaps, from - store.alignedReadStore[*itS].beginPos, fLen);
            store.alignedReadStore[*itS].endPos += fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "store.alignedReadStore[" << *itS << "].endPos += " << fLen << " == " << store.alignedReadStore[*itS].endPos << "\n";
            gapCount[*itS] += fLen;
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "gapCount[" << *itS << "] == " << gapCount[*itS] << "\n";
        }

        // Deactive done reads.
        for (TSetIt it = done.begin(); it != done.end(); ++it)
            activeReads[cl].erase(*it);
    }
 
// #if SEQAN_ENABLE_DEBUG
    {
        // Check for consistency.
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
        typedef typename TFragmentStore::TReadSeq           TReadSeq;

        TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
        for (TAlignedReadIter it2 = begin(store.alignedReadStore, Standard()); it2 != itEnd; ++it2)
        {
            typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
            TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
            SEQAN_ASSERT_EQ(length(readGaps) - length(store.readSeqStore[it2->readId]), gapCount[it2->readId]);
            if (DEBUG_INCONSISTENT_LEN)
                std::cerr << "READ GAPS\t" << (it2 - begin(store.alignedReadStore, Standard())) << "\t>>>" << readGaps << "<<< (" << length(readGaps) << ")\n"
                          << "  beginPos == " << it2->beginPos << ", endPos == " << it2->endPos << ", gapCount == " << gapCount[it2->readId] << "\n";
            if ((unsigned)abs(it2->endPos - it2->beginPos) != length(readGaps))
            {
                SEQAN_FAIL("Inconsistent begin/endPos");
            }
        }
    }
// #endif  // #if SEQAN_ENABLE_DEBUG

    return true;
}

template <typename TFragmentStore, typename TSequence, typename TCargo, typename TSetSpec, typename TSpec>
bool alignmentGraphToFragmentStore(TFragmentStore & store,
                                   seqan::Graph<seqan::Alignment<seqan::StringSet<TSequence, TSetSpec>, TCargo, TSpec> > const & g,
                                   seqan::Graph<seqan::Undirected<double> > const & distances,
                                   bool logging)
{
    using namespace seqan;
	typedef std::map<unsigned, unsigned> TComponentLength;

    // -----------------------------------------------------------------------
    // Compute connected components and get topological sorting of them.
    // -----------------------------------------------------------------------
	String<unsigned> component;
	String<unsigned> order;
	TComponentLength componentLength;
    if (empty(g))
        return true;  // Nothing to do for empty graphs.
	if (!convertAlignment(g, component, order, componentLength))
        return false;
    unsigned numComponents = length(order);

    return alignmentGraphToFragmentStore(store, g, distances, component, order, numComponents, logging);
}

// ---------------------------------------------------------------------------
// Function alignmentGraphToSmoothFragmentStore()
// ---------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TStringSet, typename TScore>
void alignmentGraphToSmoothFragmentStore(FragmentStore<TSpec, TConfig> & store,
                                         StringSet<String<unsigned> > & windowBorders,  // window i has borders windowBorders[i] and windowBorders[i + 1].
                                         Graph<Alignment<TStringSet, TScore> > & ag,
                                         seqan::Graph<seqan::Undirected<double> > const & distances,
                                         bool logging = false)
{
    {
        seqan::FragmentStore<void, MyStoreConfig_> store;
        bool b = alignmentGraphToFragmentStore(store, ag, distances, logging);
        (void)b;
        SEQAN_ASSERT(b);

        for (unsigned idx = 0; idx < length(store.contigStore); ++idx)
        {
            // Print FragmentStore.
            seqan::AlignedReadLayout layout;
            layoutAlignment(layout, store);
            __int64 l = 0;
            __int64 r = l;
            for (unsigned j = 0; j < length(layout.contigRows[idx]); ++j)
            {
                unsigned id = back(layout.contigRows[idx][j]);
                if (r < store.alignedReadStore[id].beginPos)
                    r = store.alignedReadStore[id].beginPos;
                if (r < store.alignedReadStore[id].endPos)
                    r = store.alignedReadStore[id].endPos;
            }
            if (logging)
            {
                std::cerr << "FRAGMENT STORE\n";
                printAlignment(std::cout, seqan::Raw(), layout, store, idx, l, r, 0, 1000);
            }
        }
    }

    // Compute (DAG) ABA graph from Alignment Graph.
    Graph<Directed<unsigned> > aba;
    StringSet<String<unsigned> > vertexMap;
    String<unsigned> vertexLengths;
    alignmentGraphToABruijnGraph(aba, vertexMap, vertexLengths, ag);
    SEQAN_ASSERT_EQ(lengthSum(vertexMap), numVertices(ag));
    if (DEBUG_INCONSISTENT_LEN)
    {
        { std::ofstream f3("aba.2.dot"); write(f3, aba, DotDrawing()); }
        for (unsigned i = 0; i < length(vertexMap); ++i)
            for (unsigned j = 0; j < length(vertexMap[i]); ++j)
                if (sequenceId(ag, vertexMap[i][j]) == 8)
                    std::cerr << "vertexMap[" << i << "][" << j << "] == " << vertexMap[i][j] << "\n";
    }

    // Remove small bulges from ABA graph.  Note that the vertex ids remain constant.
    unsigned smallGirth = 3;  // threshold for small bulges to remove
    Graph<Directed<unsigned> > aba2;
    StringSet<String<unsigned> > vertexMap2;
    removeSmallBulges(aba2, vertexMap2, aba, vertexMap, vertexLengths, smallGirth);
    SEQAN_ASSERT_EQ(lengthSum(vertexMap2), numVertices(ag));
    if (DEBUG_INCONSISTENT_LEN)
    {
        { std::ofstream f4("aba2.2.dot"); write(f4, aba2, DotDrawing()); }
        for (unsigned i = 0; i < length(vertexMap2); ++i)
            for (unsigned j = 0; j < length(vertexMap2[i]); ++j)
                if (sequenceId(ag, vertexMap2[i][j]) == 8)
                    std::cerr << "vertexMap2[" << i << "][" << j << "] == " << vertexMap2[i][j] << "\n";
    }

    // Separate short visiting blocks.
    /*
    unsigned smallLength = 10;
    Graph<Directed<unsigned> > aba3;
    StringSet<String<unsigned> > vertexMap3;
    separateVisitingBlocks(aba3, vertexMap3, aba2, vertexMap2, vertexLengths, ag, smallLength);
    SEQAN_ASSERT_EQ(lengthSum(vertexMap3), numVertices(ag));
    if (DEBUG_INCONSISTENT_LEN)
    {
        { std::ofstream f4("aba3.2.dot"); write(f4, aba3, DotDrawing()); }
        for (unsigned i = 0; i < length(vertexMap3); ++i)
            for (unsigned j = 0; j < length(vertexMap3[i]); ++j)
                if (sequenceId(ag, vertexMap3[i][j]) == 8)
                    std::cerr << "vertexMap3[" << i << "][" << j << "] == " << vertexMap3[i][j] << "\n";
    }
    */

    // Merge non-branching paths.
    Graph<Directed<unsigned> > aba4;
    StringSet<String<unsigned> > vertexMap4;
    mergeNonBranchingPaths(aba4, vertexMap4, aba2, vertexMap2);
    SEQAN_ASSERT_EQ(lengthSum(vertexMap4), numVertices(ag));
    if (DEBUG_INCONSISTENT_LEN)
    {
        { std::ofstream f4("aba4.2.dot"); write(f4, aba4, DotDrawing()); }
        for (unsigned i = 0; i < length(vertexMap4); ++i)
            for (unsigned j = 0; j < length(vertexMap4[i]); ++j)
                if (sequenceId(ag, vertexMap4[i][j]) == 8)
                    std::cerr << "vertexMap4[" << i << "][" << j << "] == " << vertexMap4[i][j] << "\n";
    }

    // for (typename Iterator<Graph<Directed<unsigned> >, EdgeIterator>::Type it(aba4); !atEnd(it); goNext(it))
    //     std::cerr << "EDGE\t" << sourceVertex(it) << "\t->\t" << targetVertex(it) << "\n";

    // Remove edges from ag that are in different vertices in the aba4.
    removeEdgesBetweenPartitionEntries(ag, vertexMap4);
    // write(std::cout, ag, Raw());

    // -----------------------------------------------------------------------
    // Build custom topological order of AG vertices.
    // -----------------------------------------------------------------------
    String<unsigned> abaMap;  // pointing from AG vertex to the ABA vertex
    resizeVertexMap(ag, abaMap, maxValue<unsigned>());
    for (unsigned i = 0; i < length(vertexMap4); ++i)
        for (unsigned j = 0; j < length(vertexMap4[i]); ++j)
            abaMap[vertexMap4[i][j]] = i;
    // Obtain topological order of ABA vertices.
    String<unsigned> abaOrder;
    topologicalSort(aba4, abaOrder);
    String<unsigned> abaOrderMap;
    resize(abaOrderMap, length(abaOrder));
    for (unsigned i = 0; i < length(abaOrder); ++i)
        abaOrderMap[abaOrder[i]] = i;
    // Now update abaMap for assigning an order to AG vertices.
    String<unsigned> ccMap;  // AG vertex to connected component in aba4
    reserve(ccMap, length(numVertices(ag)));
    for (unsigned i = 0; i < numVertices(ag); ++i)
    {
        abaMap[i] = abaOrderMap[abaMap[i]];
        // std::cerr << "(1) AG vertex\t" << i << "\t=>\t" << abaMap[i] << "\t";
        // std::cerr << infixWithLength(stringSet(ag)[sequenceId(ag, i)], fragmentBegin(ag, i), fragmentLength(ag, i));
        // std::cerr << "\n";
    }

    // Obtain ordinary topological order using convertAlignment().
	String<unsigned> component;
	String<unsigned> order;
    std::map<unsigned, unsigned> componentLength;
	SEQAN_CHECK(convertAlignment(ag, component, order, componentLength), "Invalid alignment!");
    unsigned numComponents = length(order);
    // Build order map from that.
    String<unsigned> orderMap;
    resize(orderMap, length(order), maxValue<unsigned>());
    // for (unsigned i = 0; i < length(order); ++i)
    //     std::cerr << "order[" << i << "] == " << order[i] << "\n";
    for (unsigned i = 0; i < length(order); ++i)
        orderMap[order[i]] = i;
    // for (unsigned i = 0; i < numVertices(ag); ++i)
    // {
    //     std::cerr << "(2) AG vertex\t" << i << "\t=>\t" << orderMap[component[i]] << "\t";
    //     std::cerr << infixWithLength(stringSet(ag)[sequenceId(ag, i)], fragmentBegin(ag, i), fragmentLength(ag, i));
    //     std::cerr << "\n";
    // }

    // Combine orders.  Primary order comes from abaOrderMap, order within the vertices comes from order.
    typedef std::map<std::pair<unsigned, unsigned>, std::vector<unsigned> > TFinalMap;
    TFinalMap finalOrderMap;
    for (unsigned i = 0; i < numVertices(ag); ++i)
        finalOrderMap[std::make_pair(abaMap[i], orderMap[component[i]])].push_back(i);
    unsigned c = 0;
    for (TFinalMap::const_iterator it = finalOrderMap.begin(); it != finalOrderMap.end(); ++it, ++c)
        for (std::vector<unsigned>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
        {            
            component[*it2] = c;
            // std::cerr << "component[" << *it2 << "] == " << c << "\t";
            // std::cerr << infixWithLength(stringSet(ag)[sequenceId(ag, *it2)], fragmentBegin(ag, *it2), fragmentLength(ag, *it2));
            // std::cerr << "\n";
        }
    for (unsigned i = 0; i < length(order); ++i)
        order[i] = i;
    numComponents = finalOrderMap.size();

    bool b = alignmentGraphToFragmentStore(store, ag, distances, component, order, numComponents, false);
    (void)b;
    SEQAN_ASSERT(b);

    // Get mapping from sequence to contig id.
    String<unsigned> seqToRef;
    resize(seqToRef, length(store.readStore), 0);
    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
        seqToRef[store.alignedReadStore[i].readId] = store.alignedReadStore[i].contigId;

    // Compute window borders 
    resize(windowBorders, length(store.contigStore));
    unsigned x = maxValue<unsigned>();
    for (TFinalMap::const_iterator it = finalOrderMap.begin(); it != finalOrderMap.end(); ++it, ++c)
    {
        unsigned c = seqToRef[sequenceId(ag, it->second.front())];
        if (empty(windowBorders[c]))
            appendValue(windowBorders[c], 0);

        if (x != it->first.first)
        {
            appendValue(windowBorders[c], fragmentLength(ag, it->second.front()));
        }
        else
        {
            back(windowBorders[c]) += fragmentLength(ag, it->second.front());
        }
        x = it->first.first;
    }
    for (unsigned j = 0; j < length(windowBorders); ++j)
    {
        std::partial_sum(begin(windowBorders[j], Standard()), end(windowBorders[j], Standard()),
                         begin(windowBorders[j], Standard()));
        if (logging)
            for (unsigned i = 0; i < length(windowBorders[j]); ++i)
                std::cerr << "  " << windowBorders[j][i] << "\n";
    }
}

// Edges between segments in different clusters in ag are removed.

template <typename TSpec, typename TConfig, typename TStringSet, typename TScore>
void alignmentGraphToSmoothFragmentStore(FragmentStore<TSpec, TConfig> & store,
                                         StringSet<String<unsigned> > & windowBorders,  // window i has borders windowBorders[i] and windowBorders[i + 1].
                                         Graph<Alignment<TStringSet, TScore> > & ag)
{
    typedef Graph<Alignment<TStringSet, TScore> > TAlignmentGraph;

    // Get distance tree between reads, used for finding connected components in AG.
    seqan::Graph<seqan::Undirected<double> > distances;
    for (unsigned i = 0; i < length(stringSet(ag)); ++i)
        addVertex(distances);
    std::set<std::pair<unsigned, unsigned> > edges;
    for (typename Iterator<TAlignmentGraph, EdgeIterator>::Type it(const_cast<TAlignmentGraph &>(ag)); !atEnd(it); goNext(it))
        edges.insert(std::make_pair(sequenceId(ag, sourceVertex(it)), sequenceId(ag, targetVertex(it))));
    for (std::set<std::pair<unsigned, unsigned> >::const_iterator it = edges.begin(); it != edges.end(); ++it)
        addEdge(distances, it->first, it->second);

    alignmentGraphToSmoothFragmentStore(store, windowBorders, ag, distances);
}

// ---------------------------------------------------------------------------
// Function mergeWindows()
// ---------------------------------------------------------------------------

// Merge windows of length < len into the larger neighbour.

inline void mergeWindows(String<unsigned> & windowBorders,
                         unsigned len)
{
    if (empty(windowBorders))
        return;

    String<unsigned> sizes;
    for (unsigned i = 0; i + 1 < length(windowBorders); ++i)
        appendValue(sizes, windowBorders[i + 1] - windowBorders[i]);

    bool added = true;
    String<unsigned> result;
    appendValue(result, sizes[0]);
    for (unsigned i = 1; i < length(sizes); ++i)
    {
        added = true;

        unsigned leftSize = 0, rightSize = 0;
        leftSize = back(result);
        if (i + 1 < length(sizes))
            rightSize = sizes[i + 1];

        if (sizes[i] >= len)
        {
            appendValue(result, sizes[i]);
            continue;
        }
        // Always true here: sizes[i] < len

        if (leftSize >= rightSize)
        {
            back(result) += sizes[i];
        }
        else
        {
            sizes[i + 1] += sizes[i];
            added = false;
        }
    }
    if (!added)
        appendValue(result, back(sizes));
    insertValue(result, 0, 0);

    std::partial_sum(begin(result, Standard()), end(result, Standard()), begin(result, Standard()));
    SEQAN_ASSERT_EQ(front(result), front(windowBorders));
    SEQAN_ASSERT_EQ(back(result), back(windowBorders));
    swap(result, windowBorders);
}

inline void mergeWindows(StringSet<String<unsigned> > & windowBorders,
                         unsigned len)
{
    for (unsigned i = 0; i < length(windowBorders); ++i)
        mergeWindows(windowBorders[i], len);
}

}  // namespace seqan

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_AG_TO_FRAG_STORE_H_
