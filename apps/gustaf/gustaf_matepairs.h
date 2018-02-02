// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2018, Kathrin Trappe, FU Berlin
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
// Author: Kathrin Trappe <kathrin.trappe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_APPS_GUSTAF_GUSTAF_MATEPAIRS_H_
#define SEQAN_APPS_GUSTAF_GUSTAF_MATEPAIRS_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include "stellar_routines.h"

using namespace seqan;

// Imports mate pairs from two files, joins them, stores the joining position in
// String readJoinPositions, and stores the sequences in the StringSet seqs and
// their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileNameL,
                 CharString const & fileNameR,
                 bool revCompl,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids,
                 StringSet<CharString> & quals,
                 String<unsigned> & readJoinPositions
                 )
{
    typedef typename Iterator<StringSet<TId>, Standard>::Type TIdSetIterator;

    seqan::SeqFileIn l, r;
    if (!open(l, toCString(fileNameL)) || !open(r, toCString(fileNameR)))
    {
        std::cerr << "Failed to open file." << std::endl;
        return false;
    }

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId id;
    TId sId;
    CharString qual;
    CharString qualL;
    CharString qualR;
    while (!atEnd(l) || !atEnd(r))
    {
        if (readRecord(id, seqL, qualL, l) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }
        // Note: overwrites id with ID of right mate, which is done also by gustaf_join_mates
        if (readRecord(id, seqR, qualR, r) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }

        appendValue(readJoinPositions, length(seqL));
        if (revCompl)
        {
            reverseComplement(seqR);
            reverse(qualR);
        }
        append(seq, seqL);
        append(seq, seqR);
        append(qual, qualL);
        append(qual, qualR);
        appendValue(seqs, seq, Generous());
        appendValue(quals, qual, Generous());

        _getShortId(sId, id);
        appendValue(ids, sId, Generous());
        clear(seq);
        clear(qual);
    }
    
    // Check for dupliacte id entries.
    StringSet<TId> uniqueIds = ids;
    std::sort(begin(uniqueIds, Standard()), end(uniqueIds, Standard()), IdComparator<TId>());  // O(n*log(n))
    TIdSetIterator itOldEnd = end(uniqueIds, Standard());
    TIdSetIterator itNewEnd = std::unique(begin(uniqueIds, Standard()), end(uniqueIds, Standard()), IdComparator<TId>());  // O(n)

    --itNewEnd;
    unsigned diff = itOldEnd - itNewEnd;
    if (length(ids) - diff > 0)
    {
        std::cout << "Found nonunique sequence IDs" << std::endl;
        return false;
    }
    std::cout << "Loaded " << length(ids) << " paired sequences" << std::endl;
    return true;
}

// /////////////////////////////////////////////////////////////////////////////
// Imports mate pairs from two files, joins them, stores the joining position in
// String readJoinPositions, and stores the sequences in the StringSet seqs and
// their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileNameL,
                 CharString const & fileNameR,
                 bool revCompl,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids,
                 String<unsigned> & readJoinPositions
                 )
{
    typedef typename Iterator<StringSet<TId>, Standard>::Type TIdSetIterator;

    seqan::SeqFileIn leftMates, rightMates;
    if (!open(leftMates, toCString(fileNameL)))
    {
        std::cerr << "Failed to open " << fileNameL << " file." << std::endl;
        return false;
    }
    if (!open(rightMates, toCString(fileNameR)))
    {
        std::cerr << "Failed to open " << fileNameR << " file." << std::endl;
        return false;
    }

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId idL;
    TId idR;
    TId sId;
    unsigned seqCount = 0;
    for (unsigned i = 0; !atEnd(leftMates) && !atEnd(rightMates); ++i, ++seqCount)
    {
        try
        {
            readRecord(idL, seqL, leftMates);
            readRecord(idR, seqR, rightMates);
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "Problem reading mates (" << ioErr.what() << ")\n";
            return false;
        }

        if (revCompl)
            reverseComplement(seqR);
        appendValue(readJoinPositions, length(seqL));
        append(seq, seqL);
        append(seq, seqR);
        // Note: saving ID of right(!) mate per default
        appendValue(seqs, seq, Generous());

        _getShortId(sId, idR);
        appendValue(ids, sId, Generous());
        clear(seq);
    }
    if (!atEnd(leftMates) || !atEnd(rightMates))
    {
        std::cerr << "Unequal number of left and right mates!\n";
        return false;
    }

    std::cout << "Loaded " << seqCount << " mate pair sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;

    StringSet<TId> uniqueIds = ids;
    // Check for dupliacte id entries.
    std::sort(begin(uniqueIds, Standard()), end(uniqueIds, Standard()), IdComparator<TId>());  // O(n*log(n))
    TIdSetIterator itOldEnd = end(uniqueIds, Standard());
    TIdSetIterator itNewEnd = std::unique(begin(uniqueIds, Standard()), end(uniqueIds, Standard()), IdComparator<TId>());  // O(n)

    --itNewEnd;
    unsigned diff = itOldEnd - itNewEnd;
    if (length(ids) - diff > 0)
    {
        std::cout << "Found nonunique sequence IDs" << std::endl;
        return false;
    }
    std::cout << "Loaded " << length(ids) << " paired sequences" << std::endl;
    return true;
}
template <typename TSequence, typename TId, typename TPos>
inline bool
_isLeftMate(StellarMatch<TSequence, TId> const & sMatch, TPos const & joiningPos)
{
    if (sMatch.begin2 > joiningPos)
        return false;
    if (sMatch.end2 < joiningPos)
        return true;
    if ((joiningPos - sMatch.begin2) > (sMatch.end2 - joiningPos))
        return true;
    return false;
}

template <typename TString, typename TId, typename TMatch>
inline bool
_countMateMatches(TString const & confMateMatches, String<TMatch> const & queryMatches, TId & dbID, unsigned support)
{
    unsigned count = 0;
    // StellarMatches in confMateMatches have the correct mate distance, but the IntervalTree does not check for correct
    // database IDs, this is done here (i.e. only matches on the same reference are counted)
    // Note: matches on the opposite strand in the correct distance are considered valid
    for (unsigned i = 0; i < length(confMateMatches); ++i)
    {
        unsigned index = confMateMatches[i];
        if (queryMatches[index].id == dbID)
            ++count;
    }
    // Support must/should dependend on the size of mateIntervalTrees!
    return count >= support;
}

// sMatch is part of the left mate, valid (right) mate matches can only be upstream but cannot exceed database length
template <typename TMatch, typename TMSplazerChain>
inline bool
_checkRightMateMatches(TMatch const & sMatch, String<TMatch> const & queryMatches, TMSplazerChain const & gustafChain, MSplazerOptions const & options)
{
    // typedef typename TMatch::TPos TPos;
    typedef unsigned TPos;
    String<unsigned> confMateMatches;
    // Begin and end position of valid region for mate matches
    unsigned dbLength = length(source(sMatch.row1));
    TPos mateIntervalBegin = _min(dbLength, sMatch.begin1 + options.libSize - gustafChain.mateJoinPosition - options.libError);
    TPos mateIntervalEnd = _min(dbLength, sMatch.end1 + options.libSize - gustafChain.mateJoinPosition + options.libError);
    SEQAN_ASSERT_GEQ_MSG(dbLength, mateIntervalBegin, "Calculated valid region for right mate exceeds database length");
    SEQAN_ASSERT_GEQ_MSG(dbLength, mateIntervalEnd, "Calculated valid region for right mate exceeds database length");
    findIntervals(confMateMatches, gustafChain.rightMateTree, mateIntervalBegin, mateIntervalEnd);
    /*
    std::cout << "checkRightMatches " << mateIntervalBegin << "," << mateIntervalEnd << ";";
    for (unsigned i = 0; i < length(confMateMatches); ++i)
        std::cout << confMateMatches[i] << ",";
    std::cout << std::endl;
    */
    return _countMateMatches(confMateMatches, queryMatches, sMatch.id, options.mateSupport);
}

// sMatch is part of the right mate, valid (left) mate matches can only be downstream but cannot be below 0
template <typename TMatch, typename TMSplazerChain>
inline bool
_checkLeftMateMatches(TMatch const & sMatch, String<TMatch> const & queryMatches, TMSplazerChain const & gustafChain, MSplazerOptions const & options)
{
    // typedef typename TMatch::TPos TPos;
    typedef unsigned TPos;
    String<unsigned> confMateMatches;
    // Begin and end position of valid region for mate matches
    TPos mateIntervalBegin = _max(0u, sMatch.begin1 - options.libSize + gustafChain.mateJoinPosition - options.libError);
    TPos mateIntervalEnd = _max(0u, sMatch.end1 - options.libSize + gustafChain.mateJoinPosition + options.libError);
    SEQAN_ASSERT_GEQ_MSG(mateIntervalBegin, 0u, "Calculated valid region for left mate is < 0");
    SEQAN_ASSERT_GEQ_MSG(mateIntervalEnd, 0u, "Calculated valid region for left mate is < 0");
    findIntervals(confMateMatches, gustafChain.leftMateTree, mateIntervalBegin, mateIntervalEnd);
    /*
    std::cout << "checkLeftMatches " << mateIntervalBegin << "," << mateIntervalEnd << ";";
    for (unsigned i = 0; i < length(confMateMatches); ++i)
        std::cout << confMateMatches[i] << ",";
    std::cout << std::endl;
    */
    return _countMateMatches(confMateMatches, queryMatches, sMatch.id, options.mateSupport);
}

// Check for one match if it is confirmed by ANY of the other (mate) matches
template <typename TMatch, typename TMSplazerChain>
inline bool
_checkMateMatches(TMatch const & sMatch, String<TMatch> const & queryMatches, TMSplazerChain const & gustafChain, MSplazerOptions const & options)
{
    if (_isLeftMate(sMatch, gustafChain.mateJoinPosition))
            return _checkRightMateMatches(sMatch, queryMatches, gustafChain, options);
    return _checkLeftMateMatches(sMatch, queryMatches, gustafChain, options);
}


// Check for two matches, each from a different mate, if they both apply to the libSize+sd, i.e. the BP is artificial
// Assumptions: sMatch1 < sMatch2, valid order and valid gap between matches regarding read sequence
template <typename TMatch, typename TMSplazerChain>
inline bool
_artificialBP(TMatch const & sMatch1, TMatch const & sMatch2, TMSplazerChain const & gustafChain)//, MSplazerOptions const & options)
{
    // Check if both matches are from different mates
    if (_isLeftMate(sMatch1, gustafChain.mateJoinPosition) && _isLeftMate(sMatch2, gustafChain.mateJoinPosition))
        return false;
    if (!_isLeftMate(sMatch1, gustafChain.mateJoinPosition) && !_isLeftMate(sMatch2, gustafChain.mateJoinPosition))
        return false;
    // Check the distance between inner match position on the reference
    /*
    typedef typename TMatch::TPos TPos;
    TPos dist = sMatch2.begin1 - sMatch1.end1;
    if (dist < (options.libSize - options.libError))
        return false;
    if (dist > (options.libSize + options.libError))
        return false;
        */
    return true;
}


// Intitialisation of graph structure for combinable StellarMatches of a read
template <typename TSequence, typename TId, typename TMSplazerChain>
void _initialiseGraphMatePairsNoBreakends(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TMSplazerChain & chain,
                      MSplazerOptions const & options)
{
    // std::cerr << " Initialising graph structure " << std::endl;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<String<StellarMatch<TSequence, TId> > >::Type TIterator;
    typedef typename TMSplazerChain::TInterval TInterval;

    String<TInterval> leftMateMatches;
    String<TInterval> rightMateMatches;
    TIterator itStellarMatches = begin(queryMatches.matches);
    TIterator itEndStellarMatches = end(queryMatches.matches);
    // The default vertex descriptor is an integer. So if inserted in the same order, the vertex descriptor value is
    // the same as the position of the corresponding vertex within the QueryMatches --> since we can easily iterate
    // through the QueryMatches and use the iterator, we wont keep track of the vertex descriptors.
    unsigned index = 0;
    for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
    {
        addVertex(chain.graph);
        // Collect intervals for mate IntervalTrees
        if (_isLeftMate(*itStellarMatches, chain.mateJoinPosition))
            appendValue(leftMateMatches, TInterval((*itStellarMatches).begin1, (*itStellarMatches).end1, index));
        else
            appendValue(rightMateMatches, TInterval((*itStellarMatches).begin1, (*itStellarMatches).end1, index));
        ++index;
    }

    /*
    std::cout << "Left " << std::endl;
    for (unsigned i = 0; i < length(leftMateMatches); ++i)
        std::cout << leftMateMatches[i].i1 << "," << leftMateMatches[i].i2 << "," << leftMateMatches[i].cargo << " ";
    std::cout << std::endl;
    std::cout << "Right " << std::endl;
    for (unsigned i = 0; i < length(rightMateMatches); ++i)
        std::cout << rightMateMatches[i].i1 << "," << rightMateMatches[i].i2 << "," << rightMateMatches[i].cargo << " ";
    std::cout << std::endl;
    */

    // Create IntervalTree for mates
    createIntervalTree(chain.rightMateTree, rightMateMatches);
    createIntervalTree(chain.leftMateTree, leftMateMatches);

    // std::cerr << " Created graph " << std::endl;
    // Add start and end to graph and property map
    chain.startVertex = addVertex(chain.graph);
    chain.endVertex = addVertex(chain.graph);

    int cargo = 0;
    resize(chain.breakpoints.slotLookupTable, 2 * length(queryMatches.matches));
    // Adding edges to start and end vertices
    for (unsigned i = 0; i < length(queryMatches.matches); ++i)
    {
        cargo = static_cast<int>(queryMatches.matches[i].begin2);
        if (cargo < (options.initGapThresh + 1))
        {
            cargo += chain.matchDistanceScores[i];
            // Note that the penalty for no mate matches is also always on the ingoing edge,
            // so for the graph initialsation, we only check this for the start vertex.
            if (!_checkMateMatches(queryMatches.matches[i], queryMatches.matches, chain, options))
                cargo += options.noMateMatchesPen;
            TEdgeDescriptor edge = addEdge(chain.graph, chain.startVertex, i, cargo);
            resizeEdgeMap(chain.breakpoints.slotLookupTable, chain.graph);
            assignProperty(chain.breakpoints, edge);
        }
        cargo = static_cast<int>(length(source(queryMatches.matches[i].row2))) -
                static_cast<int>(queryMatches.matches[i].end2);
        if (cargo < (options.initGapThresh + 1))
        {
            TEdgeDescriptor edge = addEdge(chain.graph, i, chain.endVertex, cargo);
            resizeEdgeMap(chain.breakpoints.slotLookupTable, chain.graph);
            assignProperty(chain.breakpoints, edge);
        }
    }
}

// Intitialisation of graph structure for combinable StellarMatches of a read
// including breakends
template <typename TSequence, typename TId, typename TMSplazerChain>
void _initialiseGraphMatePairs(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TId & queryId,
                      TMSplazerChain & chain,
                      MSplazerOptions const & options)
{

    // std::cerr << " Initialising graph structure " << std::endl;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<String<StellarMatch<TSequence, TId> > >::Type TIterator;
    typedef typename TMSplazerChain::TInterval TInterval;
    typedef Breakpoint<TSequence, TId> TBreakpoint;

    String<TInterval> leftMateMatches;
    String<TInterval> rightMateMatches;
    TIterator itStellarMatches = begin(queryMatches.matches);
    TIterator itEndStellarMatches = end(queryMatches.matches);
    // The default vertex descriptor is an integer. So if inserted in the same order the vertex descriptor value is
    // the same as the position of the corresponding vertex within the QueryMatches --> since we can easily iterate
    // through the QueryMatches and use the iterator we wont keep track of the vertex descriptors
    unsigned index = 0;
    for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
    {
        addVertex(chain.graph);
        // Collect intervals for mate IntervalTrees
        if (_isLeftMate(*itStellarMatches, chain.mateJoinPosition))
            appendValue(leftMateMatches, TInterval((*itStellarMatches).begin1, (*itStellarMatches).end1, index));
        else
            appendValue(rightMateMatches, TInterval((*itStellarMatches).begin1, (*itStellarMatches).end1, index));
        ++index;
    }

    // Create IntervalTree for mates
    createIntervalTree(chain.rightMateTree, rightMateMatches);
    createIntervalTree(chain.leftMateTree, leftMateMatches);

    // std::cerr << " Created graph " << std::endl;
    // Add start and end to graph and property map
    chain.startVertex = addVertex(chain.graph);
    chain.endVertex = addVertex(chain.graph);

    int cargo = 0;
    resize(chain.breakpoints.slotLookupTable, 2 * length(queryMatches.matches));
    // Adding edges to start and end vertices
    for (unsigned i = 0; i < length(queryMatches.matches); ++i)
    {
        // start vertex
        cargo = static_cast<int>(queryMatches.matches[i].begin2);
        if (cargo < (options.breakendThresh + 1))
        {
            cargo += chain.matchDistanceScores[i];
            // Note that the penalty for no mate matches is also always on the ingoing edge,
            // so for the graph initialsation, we only check this for the start vertex.
            if (!_checkMateMatches(queryMatches.matches[i], queryMatches.matches, chain, options))
                cargo += options.noMateMatchesPen;
            TEdgeDescriptor edge = addEdge(chain.graph, chain.startVertex, i, cargo);
            resizeEdgeMap(chain.breakpoints.slotLookupTable, chain.graph);
            if (cargo < (options.initGapThresh + 1))
                assignProperty(chain.breakpoints, edge);
            else
            {
                // TODO(ktrappe): needs trimming of x-drop/sloppy end
                TBreakpoint bp(queryMatches.matches[i].id,
                               queryMatches.matches[i].id,
                               queryMatches.matches[i].orientation,
                               queryMatches.matches[i].orientation,
                               queryMatches.matches[i].begin1,
                               queryMatches.matches[i].begin1,
                               queryMatches.matches[i].begin2,
                               queryMatches.matches[i].begin2,
                               queryId);
                bp.svtype = TBreakpoint::BREAKEND;
                // bp.imprecise = true;
                bp.breakend = 0; // left breakend
                assignProperty(chain.breakpoints, edge, bp);
            }
        }
        // end vertex
        cargo = static_cast<int>(length(source(queryMatches.matches[i].row2))) -
                static_cast<int>(queryMatches.matches[i].end2);
        if (cargo < (options.breakendThresh + 1))
        {
            TEdgeDescriptor edge = addEdge(chain.graph, i, chain.endVertex, cargo);
            resizeEdgeMap(chain.breakpoints.slotLookupTable, chain.graph);
            if (cargo < (options.initGapThresh + 1))
                assignProperty(chain.breakpoints, edge);
            else
            {
                TBreakpoint bp(queryMatches.matches[i].id,
                               queryMatches.matches[i].id,
                               queryMatches.matches[i].orientation,
                               queryMatches.matches[i].orientation,
                               queryMatches.matches[i].end1,
                               queryMatches.matches[i].end1,
                               queryMatches.matches[i].end2,
                               queryMatches.matches[i].end2,
                               queryId);
                bp.svtype = TBreakpoint::BREAKEND;
                // bp.imprecise = true;
                bp.breakend = 1; // right breakend
                assignProperty(chain.breakpoints, edge, bp);
            }
        }
    }
}

#endif  // #ifndef _APPS_GUSTAF_GUSTAF_MATEPAIRS_H_
