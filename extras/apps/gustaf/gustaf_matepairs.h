// ==========================================================================
//                                  Gustaf
// ==========================================================================
// Copyright (c) 2011-2013, Kathrin Trappe, FU Berlin
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

#ifndef SEQAN_EXTRAS_APPS_GUSTAF_GUSTAF_MATEPAIRS_H_
#define SEQAN_EXTRAS_APPS_GUSTAF_GUSTAF_MATEPAIRS_H_

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
                 StringSet<TId> & sIds,
                 StringSet<CharString> & quals,
                 String<unsigned> & readJoinPositions
                 )
{
    SequenceStream l(toCString(fileNameL));
    SequenceStream r(toCString(fileNameR));
    if (!isGood(l) || !isGood(r))
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
    unsigned counter = 0;
    while (!atEnd(l) || !atEnd(r))
    {
        if (readRecord(id, seqL, qualL, l) != 0)
        {
            std::cerr << "Problem reading from first input file." << std::endl;
            return false;
        }
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
        appendValue(ids, id, Generous());

        _getShortId(sId, id);
        if (!_checkUniqueId(sId, id, ids, sIds))
            ++counter;
        appendValue(sIds, sId);
        clear(seq);
        clear(qual);
    }
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
                 StringSet<TId> & sIds,
                 String<unsigned> & readJoinPositions
                 )
{
    MultiSeqFile leftMates;
    MultiSeqFile rightMates;
    if (!open(leftMates.concat, toCString(fileNameL), OPEN_RDONLY))
    {
        std::cerr << "Failed to open " << fileNameL << " file." << std::endl;
        return false;
    }
    if (!open(rightMates.concat, toCString(fileNameR), OPEN_RDONLY))
    {
        std::cerr << "Failed to open " << fileNameR << " file." << std::endl;
        return false;
    }

    AutoSeqFormat formatL;
    guessFormat(leftMates.concat, formatL);
    split(leftMates, formatL);

    AutoSeqFormat formatR;
    guessFormat(rightMates.concat, formatR);
    split(rightMates, formatR);

    unsigned seqCount = length(leftMates);
    SEQAN_ASSERT_EQ_MSG(seqCount, length(rightMates), "Unequal number of left and right mates!");
    if (seqCount != length(leftMates))
    {
        std::cerr << "Unequal number of left and right mates!" << std::endl;
        return false;
    }

    resize(readJoinPositions, seqCount);
    reserve(seqs, seqCount, Exact());
    reserve(ids, seqCount, Exact());
    reserve(sIds, seqCount, Exact());

    TSequence seq;
    TSequence seqL;
    TSequence seqR;
    TId id;
    TId sId;
    unsigned counter = 0;
    for (unsigned i = 0; i < seqCount; ++i)
    {
        assignSeq(seqL, leftMates[i], formatL);
        assignSeq(seqR, rightMates[i], formatR);
        if (revCompl)
            reverseComplement(seqR);
        readJoinPositions[i] = length(seqL);
        // std::cout << "Read " << i+1 << " has joining position " << readJoinPositions[i] << std::endl;
        append(seq, seqL);
        append(seq, seqR);
        assignSeqId(id, leftMates[i], formatL);
        appendValue(seqs, seq, Generous());
        appendValue(ids, id, Generous());

        _getShortId(sId, id);
        if (!_checkUniqueId(sId, id, ids, sIds))
            ++counter;
        appendValue(sIds, sId);
        clear(seq);
    }

    std::cout << "Loaded " << seqCount << " mate pair sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;
    if (counter > 0)
        std::cout << "Found " << counter << " nonunique sequence IDs" << std::endl;
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
    // StellarMatches in confMateMatches have to correct mate distance, but the IntervalTree does not check for correct
    // database IDs, this is done here (i.e. only matches on the same reference are counted)
    // Note: matches on the opposite strand in the correct distance are considered valid
    for (unsigned i = 0; i < length(confMateMatches); ++i)
    {
        unsigned index = confMateMatches[i];
        if (queryMatches[index].id == dbID)
            ++count;
    }
    // Support must dependent in the size of mateIntervalTrees!
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
    findIntervals(gustafChain.rightMateTree, mateIntervalBegin, mateIntervalEnd, confMateMatches);
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
    findIntervals(gustafChain.leftMateTree, mateIntervalBegin, mateIntervalEnd, confMateMatches);
    /*
    std::cout << "checkLeftMatches " << mateIntervalBegin << "," << mateIntervalEnd << ";";
    for (unsigned i = 0; i < length(confMateMatches); ++i)
        std::cout << confMateMatches[i] << ",";
    std::cout << std::endl;
    */
    return _countMateMatches(confMateMatches, queryMatches, sMatch.id, options.mateSupport);
}

template <typename TMatch, typename TMSplazerChain>
inline bool
_checkMateMatches(TMatch const & sMatch, String<TMatch> const & queryMatches, TMSplazerChain const & gustafChain, MSplazerOptions const & options)
{
    if (_isLeftMate(sMatch, gustafChain.mateJoinPosition))
            return _checkRightMateMatches(sMatch, queryMatches, gustafChain, options);
    return _checkLeftMateMatches(sMatch, queryMatches, gustafChain, options);
}



// Intitialisation of graph structure for combinable StellarMatches of a read
template <typename TSequence, typename TId, typename TMSplazerChain>
void _initialiseGraphMatePairs(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TMSplazerChain & chain,
                      MSplazerOptions const & options)
{
    // std::cerr << " Initialising graph structure " << std::endl;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<String<StellarMatch<TSequence, TId> > >::Type TIterator;
    typedef typename TMSplazerChain::TInterval TInterval;
    typedef typename TMSplazerChain::TIntervalTree TIntervalTree;

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
            if (!_checkMateMatches(queryMatches.matches[i], queryMatches.matches, chain, options))
                cargo += options.noMateMatchesPen;
            TEdgeDescriptor edge = addEdge(chain.graph, chain.startVertex, i, cargo);
            resizeEdgeMap(chain.graph, chain.breakpoints.slotLookupTable);
            assignProperty(chain.breakpoints, edge);
        }
        cargo = static_cast<int>(length(source(queryMatches.matches[i].row2))) -
                static_cast<int>(queryMatches.matches[i].end2);
        if (cargo < (options.initGapThresh + 1))
        {
            TEdgeDescriptor edge = addEdge(chain.graph, i, chain.endVertex, cargo);
            resizeEdgeMap(chain.graph, chain.breakpoints.slotLookupTable);
            assignProperty(chain.breakpoints, edge);
        }
    }
}

#endif  // #ifndef _APPS_GUSTAF_GUSTAF_MATEPAIRS_H_
