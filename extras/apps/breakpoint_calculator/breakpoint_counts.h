// ==========================================================================
//                           Breakpoint Calculator
// ==========================================================================
// Copyright (C) 2012 by Birte Kehr
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_
#define SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_

#include <seqan/align.h>
#include <seqan/random.h>

#include <lemon/lgf_writer.h>
#include <lemon/smart_graph.h>
#include <lemon/matching.h>


using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BlockInMatchinGraph
// ----------------------------------------------------------------------------

struct BlockInMatchingGraph
{
    typedef lemon::SmartGraph::Node TNode;
    TNode head;
    TNode tail;
    TNode headTelomere;
    TNode tailTelomere;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function sequencesOfBlocks()
// ----------------------------------------------------------------------------

// Thread each genome through the alignment blocks.  Choose only one random copy of repeats.

template <typename TSeqId, typename TNumber, typename TSize>
int
sequencesOfBlocks(std::map<TSeqId, StringSet<String<int> > > & blockSeqSets,
                  TNumber const seed,
                  String<std::map<CharString, String<AlignmentBlockRow<TSize, TSize> > > > & idToRowsMaps)
{
    typedef Pair<CharString, TSize> TPair;
    typedef AlignmentBlockRow<TSize, TSize> TRow;
    typedef std::map<CharString, String<TRow> >  TIdRowsMap;

    // set up random number generator
    Rng<MersenneTwister> rng(seed);

    // map of seq ids to maps of <chromosome, start position> pairs to block numbers
    typedef std::map<CharString, std::map<TPair, int> > TSeqMaps;
    TSeqMaps sequenceMaps;

    // maps of chromosomeIds per seq id
    typedef std::map<CharString, std::map<CharString, unsigned> > TChromosomeMaps;
    TChromosomeMaps chromosomeMaps;

    // fill sequenceMaps
    for (TSize i = 0; i < length(idToRowsMaps); ++i)
    {
        TIdRowsMap map = value(idToRowsMaps, i);

        typename TIdRowsMap::const_iterator it = map.begin();
        while (it != map.end())
        {
            // pick a random copy if block is repeat
            typename Value<Rng<MersenneTwister> >::Type randomIndex = 0;
            if (length(it->second) > 1)
                randomIndex = pickRandomNumber(rng) % length(it->second);

            // insert chromosomeId to chromosomeMap if not yet present
            TRow row = (*it).second[randomIndex];
            std::map<CharString, unsigned> & chrMap = chromosomeMaps[it->first];
            if (chrMap.count(row.chromosomeId) == 0)
            {
                unsigned mapSize = chrMap.size();
                chrMap[row.chromosomeId] = mapSize;
            }

            // append block number to sequenceMap
            if (row.orientation)
                sequenceMaps[it->first][TPair(row.chromosomeId, row.startPos)] = i;
            else
                sequenceMaps[it->first][TPair(row.chromosomeId, row.startPos)] = -static_cast<int>(i);

            ++it;
        }
    }

    // initialize StringSets of genomes for corresponding number of chromosomes
    typename TChromosomeMaps::const_iterator chrIt = chromosomeMaps.begin();
    while (chrIt != chromosomeMaps.end())
    {
        resize(blockSeqSets[chrIt->first], chrIt->second.size());
        ++chrIt;
    }

    // thread block seqs using sequenceMaps
    TSize i = 0;
    typename TSeqMaps::const_iterator mapIt = sequenceMaps.begin();
    while (mapIt != sequenceMaps.end())
    {
        typename std::map<TPair, int>::const_iterator seqIt = mapIt->second.begin();
        while (seqIt != (*mapIt).second.end())
        {
            unsigned chrIndex = chromosomeMaps[mapIt->first][seqIt->first.i1];
            appendValue(blockSeqSets[mapIt->first][chrIndex], seqIt->second);
            ++seqIt;
        }
        ++mapIt;
        ++i;
    }

   return 0;
}

// ----------------------------------------------------------------------------
// Function commonBlocks()                                      [2-way variant]
// ----------------------------------------------------------------------------

// Given 2 sequences of blocks, store all blocks in set that are common to both.

template <typename TBlockId, typename TSize>
void
commonBlocks(std::map<TBlockId, TSize> & set,
             StringSet<String<TBlockId> > const & seq1,
             StringSet<String<TBlockId> > const & seq2)
{
    typedef typename Iterator<ConcatenatorManyToOne<StringSet<String<TBlockId> > > >::Type TIterator;

    std::set<TBlockId> seq1Set;

    // insert all blocks of seq1 to seq1Set
    TIterator end1 = end(concat(seq1));
    for (TIterator it = begin(concat(seq1)); it != end1; ++it)
        seq1Set.insert(abs(*it));

    // insert all blocks of seq2 that are in seq1Set to set
    TIterator end2 = end(concat(seq2));
    TSize pos = 1;
    for (TIterator it = begin(concat(seq2)); it != end2; ++it)
    {
        if (seq1Set.count(abs(*it)) != 0)
        {
            set[abs(*it)] = pos;
            ++pos;
        }
    }
}

// ----------------------------------------------------------------------------
// Function commonBlocks()                                      [3-way variant]
// ----------------------------------------------------------------------------

// Given 3 sequences of blocks, store all blocks in set that are common to all three.

template <typename TBlockId, typename TSize>
void
commonBlocks(std::map<TBlockId, TSize> & set,
             StringSet<String<TBlockId> > const & seq1,
             StringSet<String<TBlockId> > const & seq2,
             StringSet<String<TBlockId> > const & seq3)
{
    typedef typename Iterator<ConcatenatorManyToOne<StringSet<String<TBlockId> > > >::Type TIterator;

    std::set<TBlockId> seq1Set, seq2Set;

    // insert all blocks of seq1 to seq1Set
    TIterator end1 = end(concat(seq1));
    for (TIterator it = begin(concat(seq1)); it != end1; ++it)
        seq1Set.insert(abs(*it));

    // insert all blocks of seq2 that are in seq1Set to set
    TIterator end2 = end(concat(seq2));
    TSize pos = 1;
    for (TIterator it = begin(concat(seq2)); it != end2; ++it)
    {
        if (seq1Set.count(abs(*it)) != 0)
        {
            set[abs(*it)] = pos;
            ++pos;
        }
        else
            seq2Set.insert(abs(*it));
    }

    // insert all blocks of seq3 that are in seq1Set or seq2Set to set
    TIterator end3 = end(concat(seq3));
    for (TIterator it = begin(concat(seq3)); it != end3; ++it)
    {
        if (set.count(abs(*it)) == 0 && (seq1Set.count(abs(*it)) != 0 || seq2Set.count(abs(*it)) != 0))
        {
            set[abs(*it)] = pos;
            ++pos;
        }
    }
}

// ----------------------------------------------------------------------------
// Function uniqueBlocks()
// ----------------------------------------------------------------------------

// Given 2 sequences of blocks, store those in set that are unique to one.

template <typename TBlockId>
void
uniqueBlocks(std::set<TBlockId> & set, String<TBlockId> const & seq1, String<TBlockId> const & seq2)
{
    typedef typename Iterator<String<TBlockId> >::Type TIterator;

    // insert all blocks of seq1 to set
    TIterator end1 = end(seq1);
    for (TIterator it1 = begin(seq1); it1 != end1; ++it1)
        set.insert(abs(*it1));

    // erase all blocks of seq2 if present from set, insert otherewise
    TIterator end2 = end(seq2);
    for (TIterator it2 = begin(seq2); it2 != end2; ++it2)
    {
        if (set.count(abs(*it2)) == 0)
            set.insert(abs(*it2));
        else
            set.erase(abs(*it2));
    }
}

// ----------------------------------------------------------------------------
// Function createMatchingGraph()
// ----------------------------------------------------------------------------

template <typename TGraph, typename TSize>
void
createMatchingGraph(TGraph & graph, String<BlockInMatchingGraph> & nodes, TSize numBlocks)
{
    resize(nodes, numBlocks);

    for (TSize i = 0; i < numBlocks; ++i)
    {
        // add nodes for one block
        value(nodes, i).head = graph.addNode();
        value(nodes, i).tail = graph.addNode();
        value(nodes, i).headTelomere = graph.addNode();
        value(nodes, i).tailTelomere = graph.addNode();

        // add telomere edges for block
        graph.addEdge(value(nodes, i).head, value(nodes, i).headTelomere);
        graph.addEdge(value(nodes, i).tail, value(nodes, i).tailTelomere);
    }
}

// ----------------------------------------------------------------------------
// Function addSequenceToMatchingGraph()
// ----------------------------------------------------------------------------

// Increase edge weights, given a sequence of blocks.

template <typename TGraph, typename TEdgeMap, typename TBlockId, typename TSize>
TSize
addSequenceToMatchingGraph(TGraph & graph,
                           TEdgeMap & eMap,
                           String<TBlockId> const & seq,
                           String<BlockInMatchingGraph> const & nMap,
                           std::map<TBlockId, TSize> & blocks)
{
    typedef typename Iterator<String<TBlockId> const>::Type TIterator;
    typedef typename TGraph::Node TNode;

    if (blocks.size() == 0)
        return 0;

    TSize weight = 0;
    TIterator it = begin(seq);
    TIterator itEnd = end(seq);

    TNode u, v, uTelo, vTelo;

    while (it != itEnd && blocks.count(abs(*it)) == 0)
        ++it;
    if (it == itEnd)
        return 0;

    // add head telomere weight
    if (*it > 0)
    {
        ++eMap[lemon::findEdge(graph, value(nMap, blocks[abs(*it)] - 1).head, value(nMap, blocks[abs(*it)] - 1).headTelomere)];
        ++weight;

        u = value(nMap, blocks[abs(*it)] - 1).tail;
        uTelo = value(nMap, blocks[abs(*it)] - 1).tailTelomere;
    }
    else
    {
        ++eMap[lemon::findEdge(graph, value(nMap, blocks[abs(*it)] - 1).tail, value(nMap, blocks[abs(*it)] - 1).tailTelomere)];
        ++weight;

        u = value(nMap, blocks[abs(*it)] - 1).head;
        uTelo = value(nMap, blocks[abs(*it)] - 1).headTelomere;
    }

    // add adjacency weights
    for (++it; it != itEnd; ++it)
    {
        if (blocks.count(abs(*it)) == 0)
            continue;

        if (*it > 0)
            v = value(nMap, blocks[abs(*it)] - 1).head;
        else
            v = value(nMap, blocks[abs(*it)] - 1).tail;

        typename TGraph::Edge edge = lemon::findEdge(graph, u, v);
        if (edge == lemon::Invalid())
        {
            edge = graph.addEdge(u, v);

            if (*it > 0)
                vTelo = value(nMap, blocks[abs(*it)] - 1).headTelomere;
            else
                vTelo = value(nMap, blocks[abs(*it)] - 1).tailTelomere;
            graph.addEdge(uTelo, vTelo);
        }
        ++eMap[edge];
        ++weight;

        if (*it > 0)
        {
            u = value(nMap, blocks[abs(*it)] - 1).tail;
            uTelo = value(nMap, blocks[abs(*it)] - 1).tailTelomere;
        }
        else
        {
            u = value(nMap, blocks[abs(*it)] - 1).head;
            uTelo = value(nMap, blocks[abs(*it)] - 1).headTelomere;
        }
    }

    // add tail telomere weight
    --it;
    while (blocks.count(abs(*it)) == 0)
        --it;
    if (*it > 0)
        vTelo = value(nMap, blocks[abs(*it)] - 1).tailTelomere;
    else
        vTelo = value(nMap, blocks[abs(*it)] - 1).headTelomere;
    ++eMap[lemon::findEdge(graph, u, vTelo)];
    ++weight;

    return weight;
}

template <typename TGraph, typename TEdgeMap, typename TBlockId, typename TSize>
TSize
addSequenceToMatchingGraph(TGraph & graph,
                           TEdgeMap & eMap,
                           StringSet<String<TBlockId> > const & seq,
                           String<BlockInMatchingGraph> & nMap,
                           std::map<TBlockId, TSize> & blocks)
{
    typedef typename Iterator<StringSet<String<TBlockId> > const>::Type TIterator;

    TSize weight = 0;
    for (TIterator it = begin(seq); it != end(seq); ++it)
        weight += addSequenceToMatchingGraph(graph, eMap, *it, nMap, blocks);

    return weight;
}

// ----------------------------------------------------------------------------
// Function pwCount()
// ----------------------------------------------------------------------------

// DEBUG code: Compute pairwise breakpoints by sorting of one sequence and counting of unique adjacencies.

template <typename TBlockId, typename TSize>
typename Size<String<TBlockId> >::Type
pwCount(String<TBlockId> const & seq1, String<TBlockId> const & seq2, std::map<TBlockId, TSize> & blocks)
{
    typedef typename Position<String<TBlockId> >::Type TPosition;
    typedef Pair<TPosition, bool> TPair;
    typename Size<String<TBlockId> >::Type count = 0;

    // store positions of seq1 blocks in map
    std::map<TBlockId, TPair> seq1Map;
    TPosition pos = 0;
    for (TPosition i = 0; i < length(seq1); ++i)
    {
        if (blocks.count(abs(seq1[i])) == 0)
            continue;
        if (seq1[i] >= 0)
            seq1Map[seq1[i]] = TPair(pos, true);
        else
            seq1Map[-seq1[i]] = TPair(pos, false);
        ++pos;
    }

    // sort seq2 according to seq1
    String<TPair> seq2Sorted;
    for (TPosition i = 0; i < length(seq2); ++i)
    {
        if (blocks.count(abs(seq2[i])) == 0)
            continue;

        bool orientation = true;
        if ((seq1Map[abs(seq2[i])].i2 && seq2[i] < 0) || (!seq1Map[abs(seq2[i])].i2 && seq2[i] > 0))
            orientation = false;

        appendValue(seq2Sorted, TPair(seq1Map[abs(seq2[i])].i1, orientation));
    }

    for (TPosition i = 1; i < length(seq2Sorted); ++i)
    {
        if (seq2Sorted[i - 1].i2)
        {
            if (seq2Sorted[i].i2)
            {
                if (seq2Sorted[i - 1].i1 + 1 != seq2Sorted[i].i1)
                    ++count;
            }
            else
            {
                ++count;
            }
        }
        else
        {
            if (!seq2Sorted[i].i2)
            {
                if (seq2Sorted[i - 1].i1 != seq2Sorted[i].i1 + 1)
                    ++count;
            }
            else
                ++count;
        }
    }

    return count;
}

// ----------------------------------------------------------------------------
// Function pairwiseCounts()
// ----------------------------------------------------------------------------

// Create matching graph for all pairs of sequences of blocks and count pairwise breakpoints using perfect matchings.

template <typename TSeqId, typename TBlockId>
typename Size<String<TBlockId> >::Type
pairwiseCounts(std::map<TSeqId, StringSet<String<TBlockId> > > & blockSeqs, bool detailed)
{
    typedef typename Size<String<TBlockId> >::Type TSize;
    typedef typename std::map<TSeqId, StringSet<String<TBlockId> > >::const_iterator TSeqIterator;
    typedef typename lemon::SmartGraph TGraph;
    typedef typename TGraph::EdgeMap<int> TEdgeMap;

    if (blockSeqs.size() < 2)
        return 0;

    if (detailed)
        std::cout << "# 2-way counts: seq1, seq2, count\n";

    TSize pairwiseBreakpoints = 0;

    // iterate over all pairs of sequences
    TSeqIterator endSeqs = --blockSeqs.end();
    for (TSeqIterator itSeq1 = blockSeqs.begin(); itSeq1 != endSeqs; ++itSeq1)
    {
        ++endSeqs;

        TSeqIterator itSeq2 = itSeq1;
        for (++itSeq2; itSeq2 != endSeqs; ++itSeq2)
        {
            std::map<TBlockId, TBlockId> blocks;
            commonBlocks(blocks, itSeq1->second, itSeq2->second);

            TGraph graph;
            TEdgeMap eMap(graph);
            String<BlockInMatchingGraph> nMap;
            TSize graphWeight = 0;

            // setup graph
            createMatchingGraph(graph, nMap, length(blocks));

            // set edge weights
            graphWeight += addSequenceToMatchingGraph(graph, eMap, itSeq1->second, nMap, blocks);
            graphWeight += addSequenceToMatchingGraph(graph, eMap, itSeq2->second, nMap, blocks);

            lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
            if (!matching.run())
            {
                std::cerr << "ERROR: Pairwise perfect Matching failed!";
            }

            SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
            TSize weightRemoved = graphWeight - matching.matchingWeight();
            if (detailed)
            {
                std::cout << itSeq1->first << "\t" << itSeq2->first << "\t" << weightRemoved;
                std::cout << std::endl;
            }

            pairwiseBreakpoints += weightRemoved;
        }
        --endSeqs;
    }

    return pairwiseBreakpoints;
}

// ----------------------------------------------------------------------------
// Function pairwiseWeightRemoved()
// ----------------------------------------------------------------------------

// Compute d_{abc}(g_a, M_{ab}) + d_{abc}(g_b, M_{ab}) (with g_a == seq1, g_b == seq2).

template <typename TBlockId, typename TSize>
TSize
pairwiseWeightRemoved(StringSet<String<TBlockId> > const & seq1,
                      StringSet<String<TBlockId> > const & seq2,
                      std::map<TBlockId, TSize> & blocks)
{
    typedef typename lemon::SmartGraph TGraph;
    typedef TGraph::EdgeMap<int> TEdgeMap;

    TGraph graph;
    TEdgeMap eMap(graph);
    String<BlockInMatchingGraph> nMap;
    TSize graphWeight = 0;

    // setup graph
    createMatchingGraph(graph, nMap, length(blocks));

    // set edge weights
    graphWeight += addSequenceToMatchingGraph(graph, eMap, seq1, nMap, blocks);
    graphWeight += addSequenceToMatchingGraph(graph, eMap, seq2, nMap, blocks);

    lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
    if (!matching.run())
    {
        std::cerr << "ERROR: Perfect matching failed while computing triplet breakpoint counts!\n";
        return maxValue<TSize>();
    }

    SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
    TSize weightRemoved = graphWeight - matching.matchingWeight();

    return weightRemoved;
}

// ----------------------------------------------------------------------------
// Function tripletWeightRemoved()
// ----------------------------------------------------------------------------

// Compute d_M.

template <typename TBlockId, typename TSize>
TSize
tripletWeightRemoved(StringSet<String<TBlockId> > const & seq1,
                     StringSet<String<TBlockId> > const & seq2,
                     StringSet<String<TBlockId> > const & seq3,
                     std::map<TBlockId, TSize> & blocks)
{
    typedef typename lemon::SmartGraph TGraph;
    typedef TGraph::EdgeMap<int> TEdgeMap;

    TGraph graph;
    TEdgeMap eMap(graph);
    String<BlockInMatchingGraph> nMap;
    TSize graphWeight = 0;

    // setup graph
    createMatchingGraph(graph, nMap, length(blocks));

    // set edge weights
    graphWeight += addSequenceToMatchingGraph(graph, eMap, seq1, nMap, blocks);
    graphWeight += addSequenceToMatchingGraph(graph, eMap, seq2, nMap, blocks);
    graphWeight += addSequenceToMatchingGraph(graph, eMap, seq3, nMap, blocks);

    lemon::MaxWeightedPerfectMatching<TGraph, TEdgeMap> matching(graph, eMap);
    if (!matching.run())
    {
        std::cerr << "ERROR: Perfect matching failed while computing triplet breakpoint counts!\n";
        return maxValue<TSize>();
    }

    SEQAN_ASSERT_LEQ(matching.matchingWeight(), static_cast<__int64>(graphWeight));
    TSize weightRemoved = graphWeight - matching.matchingWeight();

    return weightRemoved;
}

// ----------------------------------------------------------------------------
// Function tripletCounts()
// ----------------------------------------------------------------------------

// Create matching graph for all triplets of sequences of blocks and count three-way hidden breakpoints using perfect
// matchings.

template <typename TSeqId, typename TBlockId>
double
tripletCounts(std::map<TSeqId, StringSet<String<TBlockId> > > & blockSeqs, bool detailed)
{
    typedef typename std::map<TSeqId, StringSet<String<TBlockId> > >::const_iterator TSeqIterator;
    typedef typename Size<String<TBlockId> >::Type TSize;

    if (blockSeqs.size() < 3)
        return 0;

    if (detailed)
        std::cout << "# 3-way counts: seq1, seq2, seq3, count\n";

    double tripletCounts = 0;

    // iterate over all triplets of sequences
    TSeqIterator endSeqs = ---- blockSeqs.end();
    for (TSeqIterator itSeq1 = blockSeqs.begin(); itSeq1 != endSeqs; ++itSeq1)
    {
        ++endSeqs;
        TSeqIterator itSeq2 = itSeq1;
        for (++itSeq2; itSeq2 != endSeqs; ++itSeq2)
        {
            ++endSeqs;
            TSeqIterator itSeq3 = itSeq2;
            for (++itSeq3; itSeq3 != endSeqs; ++itSeq3)
            {
                std::map<TBlockId, TSize> blocks;
                commonBlocks(blocks, itSeq1->second, itSeq2->second, itSeq3->second);

                TSize count12 = pairwiseWeightRemoved(itSeq1->second, itSeq2->second, blocks);
                TSize count13 = pairwiseWeightRemoved(itSeq1->second, itSeq3->second, blocks);
                TSize count23 = pairwiseWeightRemoved(itSeq2->second, itSeq3->second, blocks);
                TSize count123 = tripletWeightRemoved(itSeq1->second, itSeq2->second, itSeq3->second, blocks);

                SEQAN_ASSERT_GEQ(count123, (count12 + count13 + count23) / 2);

                double count = count123 - (count12 + count13 + count23) / 2.0;
                if (detailed)
                    std::cout << itSeq1->first << "\t" << itSeq2->first << "\t" << itSeq3->first << "\t" << count << "\n";
                tripletCounts += count;
            }
            --endSeqs;
        }
        --endSeqs;
    }

    return tripletCounts;
}

#endif  // #ifndef SANDBOX_BKEHR_APPS_BREAKPOINT_CALCULATOR_BREAKPOINT_COUNTS_H_
