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

#ifndef SEQAN_APPS_GUSTAF_MSPLAZER_ALGORITHMS_H_
#define SEQAN_APPS_GUSTAF_MSPLAZER_ALGORITHMS_H_

#include "msplazer.h"
// #include <seqan/PathEnumeration.h>
// #include "../../../andreotti/include/seqan/PathEnumeration.h"
#include <seqan/align_split.h>
#include "gustaf_matepairs.h"

using namespace seqan;

// Check for match overlap
template <typename TPos>
inline bool _checkMatchOverlap(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin, TPos const & m2End)
{
    // if overlap: begin position of snd match is smaller than end position of fst match
    if (m2Begin < m1End && m1End < m2End && m1Begin < m2Begin)
        return true;

    return false;
}

// Check for match similarity: The overlapping part has to be smaller than a specified percentage of each match length
template <typename TPos>
inline bool _checkMatchSim(TPos const & m1Begin,
                           TPos const & m1End,
                           TPos const & m2Begin,
                           TPos const & m2End,
                           MSplazerOptions const & msplazerOptions)
{
    // TPos overlapPartLength = m1End - m2Begin;
    double overlapPartLength = static_cast<double>(m1End - m2Begin);
    // Catch special case of required overlap=0
    if (msplazerOptions.simThresh == static_cast<double>(0.0))
        return overlapPartLength == static_cast<double>(0.0);

    // TPos match1Length = m1End - m1Begin;
    // TPos match2Length = m2End - m2Begin;
    double match1Length = static_cast<double>(m1End - m1Begin);
    double match2Length = static_cast<double>(m2End - m2Begin);
    // Check if overlapping percent of each match is lower than the allowed percent threshold
    double const EPSILON = 0.00001;

    if (((1.0 * overlapPartLength / match1Length) < (msplazerOptions.simThresh + EPSILON))
       && ((1.0 * overlapPartLength / match2Length) < (msplazerOptions.simThresh + EPSILON)))
        return true;
    return false;
}

// Check distance between matches:
// The maximal allowed distance between two matches has to be smaller than a specified threshold value
// template < typename TSequence, typename TId >
template <typename TPos>
inline bool _checkMatchDist(TPos const & begin, TPos const & end, MSplazerOptions const & msplazerOptions)
{
    // Check if distance between matches is smaller than the distance threshold, assumes that begin > end
    if ((int) (begin - end) < (msplazerOptions.gapThresh + 1))
        return true;

    return false;
}

// Check matches for same database
template <typename TSequence, typename TId>
inline bool _checkDBIds(StellarMatch<TSequence, TId> const & match1, StellarMatch<TSequence, TId> const & match2)
{
    return match1.id == match2.id;
}

// Check matches for same strand
template <typename TSequence, typename TId>
inline bool _checkMatchStrands(StellarMatch<TSequence, TId> const & match1, StellarMatch<TSequence, TId> const & match2)
{
    return match1.orientation == match2.orientation;
}

// Check order in matching strand: Function assumes that both matches are on the same strand and in the same genome
// and that match1.begin2 < match2.begin2 (ordered matches according to read)
template <typename TSequence, typename TId>
inline bool _checkMatchOrderInDB(StellarMatch<TSequence, TId> const & match1,
                                 StellarMatch<TSequence, TId> const & match2)
{
    return match1.begin1 < match2.begin1;
}

// Check match compatibility: Checks compatibility property of two stellar matches: Do they overlap (in a specified way)?
// If not, are they still close enough? Returns false, if not.
template <typename TPos>
inline bool _checkMatchComp(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin, TPos const & m2End,
                            bool & doBP, bool & insertEdge, MSplazerOptions const & msplazerOptions)
{
    doBP = false;
    insertEdge = false;
    // Check for true overlap
    if (_checkMatchOverlap(m1Begin, m1End, m2Begin, m2End))
    {
        // Check for similarity, i.e. the percentage of the overlapping part to the match length
        doBP = _checkMatchSim(m1Begin, m1End, m2Begin, m2End, msplazerOptions);
        insertEdge = doBP;
        return true;
    }
    // If not overlapping correctly, check if there is a gap and then the gap length (matchDist)
    if (m1End <= m2Begin)
    {
        insertEdge = _checkMatchDist(m2Begin, m1End, msplazerOptions);
        // If gap length is small enough, next match might still be ok. If gap is too big, the next one will be as well.
        return insertEdge;
    }
    return true;
}

// Check for a valid overlap of m1 and m2 given their begin and end positions. Valid are
// m1begin <= m2begin < m2end <= m1end
// m1begin <= m2begin < m1end <= m2end
// with minimal length constraint of 50bp for the tandem repeat length
template <typename TPos>
inline bool _isTandemOverlap(TPos const & m1Begin, TPos const & m1End, TPos const & m2Begin, int tandemThresh)
{
    if (m2Begin < m1Begin)
        return false;
    if (m1End < m2Begin)
        return false;
    if (m1End-m2Begin < static_cast<unsigned>(tandemThresh)) // length/size of tandem repeat (default 50)
        return false;

    return true;
}

// Intitialisation of graph structure for combinable StellarMatches of a read
/*
template <typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TVertexDescriptor,
          typename TBreakpointMap>
void _initialiseGraph(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TGraph & graph,
                      TScoreAlloc & matchDistanceScores,
                      TVertexDescriptor & startVertex,
                      TVertexDescriptor & endVertex,
                      TBreakpointMap & queryBreakpoints,
                      MSplazerOptions const & msplazerOptions)
*/
template <typename TSequence, typename TId, typename TMSplazerChain>
void _initialiseGraphNoBreakend(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TMSplazerChain & chain,
                      MSplazerOptions const & options)
{
    // std::cerr << " Initialising graph structure " << std::endl;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<String<StellarMatch<TSequence, TId> > >::Type TIterator;

    TIterator itStellarMatches = begin(queryMatches.matches);
    TIterator itEndStellarMatches = end(queryMatches.matches);
    // The default vertex descriptor is an integer. So if inserted in the same order the vertex descriptor value is
    // the same as the position of the corresponding vertex within the QueryMatches --> since we can easily iterate
    // through the QueryMatches and use the iterator we wont keep track of the vertex descriptors
    for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
        addVertex(chain.graph);

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
template <typename TSequence, typename TId, typename TMSplazerChain>
void _initialiseGraph(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                      TId & queryId,
                      TMSplazerChain & chain,
                      MSplazerOptions const & options)
{
    // std::cerr << " Initialising graph structure " << std::endl;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename Iterator<String<StellarMatch<TSequence, TId> > >::Type TIterator;
    typedef Breakpoint<TSequence, TId> TBreakpoint;

    TIterator itStellarMatches = begin(queryMatches.matches);
    TIterator itEndStellarMatches = end(queryMatches.matches);
    // The default vertex descriptor is an integer. So if inserted in the same order the vertex descriptor value is
    // the same as the position of the corresponding vertex within the QueryMatches --> since we can easily iterate
    // through the QueryMatches and use the iterator we wont keep track of the vertex descriptors
    for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
        addVertex(chain.graph);

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

// Match Chaining for one query: Inserts edges between compatible matches and determines their breakpoint
template <typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TBreakpointMap, typename TMSplazerChain>
void _chainMatches(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                   TId const & queryId,
                   TSequence & query,
                   TGraph & graph,
                   TScoreAlloc & matchDistanceScores,
                   TBreakpointMap & queryBreakpoints,
                   TMSplazerChain & chain,
                   MSplazerOptions const & msplazerOptions)
{
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Align<TSequence, ArrayGaps> TAlign;
    typedef typename Position<TSequence>::Type  TPos;
    typedef Breakpoint<TSequence, TId> TBreakpoint;
    typedef typename Infix<TSequence>::Type TInfix;

    // Output values for compatibility check: do breakpoint evaluation, insert edge into graph,
    // found gap --> stop iterating if gap is too big (overlap in read)
    bool doBP = false;
    bool insertEdge = false;
    // refOrder = true meaning matches have same reference order (duplication indicator)
    bool  refOrder = true;
    // Penalties
    int diffDBPen = 0, diffStrandPen = 0, diffOrderPen = 0, noMateMatchesPen = 0;
    // Terminating condition for taking the next snd match for comparison:
    // takeNextMatch == false meaning this and all the next matches are too far away
    bool takeNextMatch = true;
    // Breakpoint parameters
    // edit distance
    Score<int> scoreType(0, -1, -1, -1);
    // edit distance (affine)
    // Score<int> scoreType(0, -1, -4, -1);
    int splitPos = 0;

    // Cargo on edges
    int cargo = 0;
    // Breakpoint score
    int score = 0;

    // std::cout << "Query Id: " << queryId << std::endl;

    // loop over all query matches
    // std::cerr << "In chainQueryMatches length(queryMatches.matches): " << length(queryMatches.matches) << std::endl;
    for (unsigned m1 = 0; m1 < (length(queryMatches.matches) - 1); ++m1)
    {
        StellarMatch<TSequence, TId> & stMatch1 = queryMatches.matches[m1];
        TPos m1Begin = stMatch1.begin2;
        TPos m1End = stMatch1.end2;
        takeNextMatch = true;
        // loop over all query matches that are supposed to be compatible
        for (unsigned m2 = m1 + 1; takeNextMatch && (m2 < length(queryMatches.matches)); ++m2)
        {
            // /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Compatibility check
            StellarMatch<TSequence, TId> & stMatch2 = queryMatches.matches[m2];
            TPos m2Begin = stMatch2.begin2;
            TPos m2End = stMatch2.end2;
            // std::cout << "chain matches " << m1Begin << " " << m1End << " " << m2Begin << " " << m2End << std::endl;
            // Returns false, if next match is def. not compatible anymore
            takeNextMatch = _checkMatchComp(m1Begin, m1End, m2Begin, m2End, doBP, insertEdge, msplazerOptions);
	    refOrder = true;

            //std::cout << "insertEdge: " << insertEdge << " doBP: " << doBP << std::endl;
            // match is compatible
            if (insertEdge)
            {
                // /////////////////////////////////////////////////////////////////////////////////////////////////////
                // Penalties
                // Different reference sequence penalty
                diffDBPen = _checkDBIds(stMatch1, stMatch2) ? 0 : msplazerOptions.diffDBPen;
                // Different orientation penalty
                diffStrandPen =
                    (diffDBPen > 0 || _checkMatchStrands(stMatch1, stMatch2)) ? 0 : msplazerOptions.diffStrandPen;
                // Different order in reference than in read penalty
		refOrder = _checkMatchOrderInDB(stMatch1, stMatch2);
                diffOrderPen =
                    (diffDBPen > 0 || diffStrandPen > 0 || refOrder) ? 0 : msplazerOptions.diffOrderPen;
                // Penalty if there are no confirming mate matches
                noMateMatchesPen =
                    _checkMateMatches(stMatch2, queryMatches.matches, chain, msplazerOptions) ? 0 : msplazerOptions.noMateMatchesPen;
                // Compute edge cargo (edge weight)
                cargo = static_cast<int>(matchDistanceScores[m2])
                        + diffDBPen
                        + diffStrandPen
                        + diffOrderPen
                        + noMateMatchesPen;

                // /////////////////////////////////////////////////////////////////////////////////////////////////////
                // Breakpoint computation, Graph input
                // Do breakpoint computation if flag doBP was set 'true' during match compatibility check

                // Compute breakpoint
                TPos startSeqPos, endSeqPos, readStartPos, readEndPos;
                // if match1 and match2 are from different mates and in right distance, doBP=false
                if (doBP)
                {
                    // Create alignments from Stellar rows as input for breakpoint function
                    TAlign match1, match2;
                    resize(rows(match1), 2);
                    resize(rows(match2), 2);

                    assignSource(row(match1, 0), infix(source(stMatch1.row2), m1Begin, m2End));  // Read infix match1
                    assignSource(row(match1, 1), infix(source(stMatch1.row1), stMatch1.begin1, stMatch1.end1)); // Ref infix
                    assignSource(row(match2, 0), infix(source(stMatch2.row2), m1Begin, m2End));  // Read infix match2
                    assignSource(row(match2, 1), infix(source(stMatch2.row1), stMatch2.begin1, stMatch2.end1)); // Ref infix

                    // Reverse complement matches if they are on the reverse strand to get sequence content right
                    if (!stMatch1.orientation)
                        reverseComplement(source(row(match1, 1)));
                    if (!stMatch2.orientation)
                        reverseComplement(source(row(match2, 1)));

                    // Compute breakpoint score
                    // The resulting score is the sum of the scores of both alignments. --> substract old match scores
                    // Note: old match scores are already distances, new score is a negative score bc. we use scoring sceme (0, -1, -1, -1)
                    auto tmp = _splitAlignmentImpl(row(match1, 0), row(match1, 1), row(match2, 0), row(match2, 1),
                                                   scoreType, AlignConfig<false, false, true, true>());
                    SEQAN_ASSERT_NEQ(std::get<0>(tmp) + std::get<1>(tmp), std::numeric_limits<int>::max());
                    score = std::get<0>(tmp) + std::get<1>(tmp) + (static_cast<int>(matchDistanceScores[m1]) + static_cast<int>(matchDistanceScores[m2]));
                    splitPos = endPosition(row(match1, 0)) + stMatch1.begin2;

                    // Refine cargo by reducing distance by score, score is the number of edit errors avoided by the breakpoint/trimming
                    cargo -= score;
                    readStartPos = splitPos;
                    readEndPos = splitPos;
                    startSeqPos = endPosition(row(match1, 1)) + stMatch1.begin1;
                    endSeqPos = beginPosition(row(match2, 1)) + stMatch2.begin1;

                    // Reverse complement matches back again
                    if (!stMatch1.orientation)
                        reverseComplement(source(row(match1, 1)));
                    if (!stMatch2.orientation)
                        reverseComplement(source(row(match2, 1)));

                }
                else // No overlap but a valid gap, bp corresponds to start and end positions of the matches
                {
                    // Refine score by adding gap length to distance
                    cargo += static_cast<int>(m2Begin - m1End);
                    startSeqPos = stMatch1.end1;
                    endSeqPos = stMatch2.begin1;
                    readStartPos = stMatch1.end2;
                    readEndPos = stMatch2.begin2;
                }

                // Adjust bp position if they arise from matches on the reverse strand
                if (!stMatch1.orientation)
                    startSeqPos = stMatch1.end1 - startSeqPos + stMatch1.begin1;
                if (!stMatch2.orientation)
                    endSeqPos = stMatch2.end1 - endSeqPos + stMatch2.begin1;

                // Create breakpoint with calculated positions and match information
                TBreakpoint bp(stMatch1.id,
                               stMatch2.id,
                               stMatch1.orientation,
                               stMatch2.orientation,
                               startSeqPos,
                               endSeqPos,
                               readStartPos,
                               readEndPos,
                               queryId);

                // Imprecise breakpoint?
                // if (!doBP) bp.imprecise = true;

		//std::cout << bp << std::endl;
		//std::cout << stMatch1 << stMatch2 << std::endl;
                // Returns true for insertion type, get insertion infix then
                if (setSVType(bp, refOrder))
                {
                    if (stMatch1.end2 < stMatch2.begin2)
                    {
                        // TSequence inSeq;
                        TInfix inSeq;
                        // get insertion sequence from matches and read sequence --> infix(endPos(match1),startPos(match2))
                        inSeq = infix(query, stMatch1.end2, stMatch2.begin2);
                        setInsertionSeq(bp, inSeq);
                    }
                }
		//std::cout << bp << std::endl;
                // TODO(ktrappe): needs adjustment of positions?
                if (bp.svtype == TBreakpoint::DISPDUPLICATION && _isTandemOverlap(stMatch1.begin1, startSeqPos, endSeqPos, msplazerOptions.tandemThresh))
                {
                    // Double overlap check (not handled jet)
                    // std::cerr << "double overlap in reference and read called from read overlap" << std::endl;
                    // std::cout << "Translocation double overlap" << std::endl;
                    bp.svtype = TBreakpoint::SEQAN_TANDEM;
                }

                //std::cout << "Breakpoint " << bp << std::endl;
                // Insert breakpoint
                TEdgeDescriptor edge = addEdge(graph, m1, m2, cargo);
                resizeEdgeMap(queryBreakpoints.slotLookupTable, graph);
                // if match1 and match2 are from different mates and BP is indel check for artificial bp
                //if (bp.svtype == TBreakpoint::DELETION && _artificialBP(stMatch1, stMatch2, chain, msplazerOptions))
                if (_artificialBP(stMatch1, stMatch2, chain))
                    assignProperty(queryBreakpoints, edge);
                else
                    assignProperty(queryBreakpoints, edge, bp);
            }
            doBP = false;
            insertEdge = false;
        }
    }
}

// Match Chaining for one query: Inserts edges between compatible matches and determines their breakpoint
template <typename TSequence, typename TId, typename TGraph, typename TScoreAlloc, typename TBreakpointMap, typename TMSplazerChain>
void _chainMatchesReference(QueryMatches<StellarMatch<TSequence, TId> > & queryMatches,
                            TId const & queryId,
                            TSequence & query,
                            TGraph & graph,
                            TScoreAlloc & matchDistanceScores,
                            TBreakpointMap & queryBreakpoints,
                            TMSplazerChain & chain,
                            MSplazerOptions const & msplazerOptions)
{
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef Align<TSequence, ArrayGaps> TAlign;
    typedef typename Position<TSequence>::Type  TPos;
    typedef Breakpoint<TSequence, TId> TBreakpoint;
    typedef StellarMatch<TSequence, TId> TMatch;
    typedef typename Infix<TSequence>::Type TInfix;

    // Penalties
    int diffStrandPen = 0;
    int diffOrderPen = 0;
    int noMateMatchesPen = 0;
    // Output values for compatibility check: do breakpoint evaluation, insert edge into graph,
    // found gap --> stop iterating if gap is too big (overlap in read)
    bool doBP = false;
    bool insertEdge = false;
    // swap=true meaning matches have different order in read
    bool swap = false;
    // refOrder=true meaning matches have different order in reference
    bool refOrder = false;
    // Breakpoint parameters
    // edit distance
    // edit distance with affine gap cost (can improve alignment around split points in present of gaps)
    // Score<int> scoreType(0, -1, -2, -1);
    Score<int> scoreType(0, -1, -1, -1);
    int splitPos = 0;
    // Cargo on edges
    int cargo = 0;

    // Breakpoint score (gain of edit distance)
    int score = 0;

    //std::cout << "Query Id: " << queryId << std::endl;

    // loop over all query matches
    for (unsigned m1 = 0; m1 < (length(queryMatches.matches) - 1); ++m1)
    {
        TMatch * stMatch1 = &queryMatches.matches[m1];
        TPos m1Begin = (*stMatch1).begin1;
        TPos m1End = (*stMatch1).end1;
        // loop over all query matches that are supposed to be compatible
        for (unsigned m2 = m1 + 1; m2 < length(queryMatches.matches); ++m2)
        {
            swap = false;
	    refOrder = true;
            // /////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Compatibility check
            TMatch * stMatch2 = &queryMatches.matches[m2];

            // std::cerr << "stMatch1 and stMatch2: " << *stMatch1 << *stMatch2 << std::endl;
            // Check if on same chromosome
            if (_checkDBIds(*stMatch1, *stMatch2))
            {
                // Matches have other order in reference than in read
                // Exchange matches once, via temp reference, then change back at the end!
                if (m1Begin > (*stMatch2).begin1)
                {
                    std::swap(stMatch1, stMatch2);
                    m1Begin = (*stMatch1).begin1;
                    m1End = (*stMatch1).end1;
                    swap = true;
                }

                TPos m2Begin = (*stMatch2).begin1;
                TPos m2End = (*stMatch2).end1;
                _checkMatchComp(m1Begin, m1End, m2Begin, m2End, doBP, insertEdge, msplazerOptions);

                if (insertEdge)
                {
                    // /////////////////////////////////////////////////////////////////////////////////////////////////
                    // Penalties
                    // Insertion penalty
                    // insertionPen = msplazerOptions.insertionPen + insertionLength;
                    // Different orientation of matches
                    diffStrandPen = (_checkMatchStrands(*stMatch1, *stMatch2)) ? 0 : msplazerOptions.diffStrandPen;
                    // Different order in reference than in read penalty
                    diffOrderPen = 0;
		    refOrder = _checkMatchOrderInDB(*stMatch1, *stMatch2);
                    if (swap)
                        diffOrderPen =
                            (diffStrandPen > 0 || refOrder) ? 0 : msplazerOptions.diffOrderPen;
                    // Penalty if there are no confirming mate matches
                    noMateMatchesPen =
                        _checkMateMatches(*stMatch2, queryMatches.matches, chain, msplazerOptions) ? 0 : msplazerOptions.noMateMatchesPen;

                    // /////////////////////////////////////////////////////////////////////////////////////////////////
                    // Breakpoint computation, Graph input
                    // Do breakpoint computation if flag doBP was set 'true' during match compatibility check

                    // Compute breakpoint
                    TPos startSeqPos, endSeqPos, readStartPos, readEndPos;
                    if (doBP)
                    {
                        // Create alignments from Stellar rows as input for breakpoint function
                        TAlign match1, match2;
                        resize(rows(match1), 2);
                        resize(rows(match2), 2);

                        assignSource(row(match1, 0), infix(source((*stMatch1).row1), m1Begin, m2End)); // Ref infix match1
                        assignSource(row(match1, 1),
                                     infix(source((*stMatch1).row2), (*stMatch1).begin2, (*stMatch1).end2)); // Read infix
                        assignSource(row(match2, 0), infix(source((*stMatch2).row1), m1Begin, m2End)); // Ref infix match2
                        assignSource(row(match2, 1),
                                     infix(source((*stMatch2).row2), (*stMatch2).begin2, (*stMatch2).end2)); // Read infix

                        // Reverse complement sequence content from matches on reverse strand
                        if (!(*stMatch1).orientation)
                            reverseComplement(source(row(match1, 1)));
                        if (!(*stMatch2).orientation)
                            reverseComplement(source(row(match2, 1)));

                        // Compute breakpoint and score
                        // int lDiag = -10, uDiag = 10;
                        // score = splitAlignment(match1, match2, scoreType, lDiag, uDiag);
                        // score = splitAlignment(match1, match2, scoreType);
                        auto tmp = _splitAlignmentImpl(row(match1, 0), row(match1, 1), row(match2, 0), row(match2, 1),
                                                    scoreType, AlignConfig<false, false, true, true>());
                        score = std::get<0>(tmp) + std::get<1>(tmp);
                        SEQAN_ASSERT_EQ(endPosition(row(match1, 0)), beginPosition(row(match2, 0)));

                        // Compute cargo, reduce distance by score
                        // cargo = static_cast<int>(matchDistanceScores[m2]) - score + diffStrandPen + diffOrderPen;
                        cargo = static_cast<int>(matchDistanceScores[m1]) + static_cast<int>(matchDistanceScores[m2]) +
                                score + diffStrandPen + diffOrderPen + noMateMatchesPen;

                        // If matches have been swapped, i.e. their order has been switched, bp positions have to be computed from the other match and vice versa
                        if (!swap)
                        {
                            // Get view position in stMatch1.row2 of first source character after split
                            splitPos = endPosition(row(match1, 0)) + (*stMatch1).begin1;
                            readStartPos = endPosition(row(match1, 1)) + (*stMatch1).begin2;
                            readEndPos = beginPosition(row(match2, 1)) + (*stMatch2).begin2;
                        }
                        else
                        {
                            // Get view position in stMatch2.row2 of first source character after split
                            splitPos = endPosition(row(match2, 0)) + (*stMatch2).begin1;
                            readStartPos = endPosition(row(match2, 1)) + (*stMatch2).begin2;
                            readEndPos = beginPosition(row(match1, 1)) + (*stMatch1).begin2;
                        }
                        startSeqPos = splitPos;
                        endSeqPos = splitPos;

                        if (!swap)
                        {
                            // Adjust bp position if they arise from matches on the reverse strand
                            if (!(*stMatch1).orientation)
                                startSeqPos = (*stMatch1).end1 - startSeqPos + (*stMatch1).begin1;
                            if (!(*stMatch2).orientation)
                                endSeqPos = (*stMatch2).end1 - endSeqPos + (*stMatch2).begin1;
                        }
                        else
                        {
                            // Adjust bp position if they arise from matches on the reverse strand
                            if (!(*stMatch2).orientation)
                                startSeqPos = (*stMatch2).end1 - startSeqPos + (*stMatch2).begin1;
                            if (!(*stMatch1).orientation)
                                endSeqPos = (*stMatch1).end1 - endSeqPos + (*stMatch1).begin1;
                        }

                        // Reverse complement sequence content back again
                        if (!(*stMatch1).orientation)
                            reverseComplement(source(row(match1, 1)));
                        if (!(*stMatch2).orientation)
                            reverseComplement(source(row(match2, 1)));
                    }
                    else // No overlap but a valid gap, bp corresponds to start and end positions of the matches
                    {
                        // Compute score, add gap length to distance
                        cargo = static_cast<int>(matchDistanceScores[m2]) +
                                static_cast<int>(m2Begin - m1End) + diffStrandPen + diffOrderPen + noMateMatchesPen;
                        if (!swap)
                        {
                            startSeqPos = (*stMatch1).end1;         // m1End
                            endSeqPos = (*stMatch2).begin1;
                            readStartPos = (*stMatch1).end2;
                            readEndPos = (*stMatch2).begin2;
                            // Adjust bp position if they arise from matches on the reverse strand
                            if (!(*stMatch1).orientation)
                                startSeqPos = (*stMatch1).end1 - startSeqPos + (*stMatch1).begin1;
                            if (!(*stMatch2).orientation)
                                endSeqPos = (*stMatch2).end1 - endSeqPos + (*stMatch2).begin1;
                        }
                        else
                        {
                            startSeqPos = (*stMatch2).end1;         // m1End
                            endSeqPos = (*stMatch1).begin1;
                            readStartPos = (*stMatch2).end2;
                            readEndPos = (*stMatch1).begin2;
                            // Adjust bp position if they arise from matches on the reverse strand
                            if (!(*stMatch2).orientation)
                                startSeqPos = (*stMatch2).end1 - startSeqPos + (*stMatch2).begin1;
                            if (!(*stMatch1).orientation)
                                endSeqPos = (*stMatch1).end1 - endSeqPos + (*stMatch1).begin1;
                        }

                    }

                    // Create breakpoint with calculated positions and match information
                    TBreakpoint bp((*stMatch1).id,
                                   (*stMatch2).id,
                                   (*stMatch1).orientation,
                                   (*stMatch2).orientation,
                                   startSeqPos,
                                   endSeqPos,
                                   readStartPos,
                                   readEndPos,
                                   queryId);

                    // Imprecise breakpoint?
                    // if (!doBP) bp.imprecise = true;

                    // Set SV type of breakpoint, returns true if SV type is "insertion", if so, compute inserted sequence and assign to bp
                    if (setSVType(bp, refOrder))
                    {
                        // TSequence inSeq;
                        TInfix inSeq;
                        if (readStartPos < readEndPos)
                            inSeq = infix(query, readStartPos, readEndPos);
                        else
                            inSeq = infix(query, readEndPos, readStartPos);
                        if (length(inSeq) == 0)
                            setSVType(bp, TBreakpoint::INVALID);
                        else
                            setInsertionSeq(bp, inSeq);
                    }

                    // Put breakpoint on corresponding edge in breakpoint graph, overwrite existing bp if new one better (smaller cargo)
                    bool artificialBP = _artificialBP(*stMatch1, *stMatch2, chain);
                    TEdgeDescriptor edge;
                    edge = findEdge(graph, m1, m2);
                    // Returns 0 if edge does not exist
                    if (edge != 0 && !artificialBP)
                    {
                        TBreakpoint bp_dummy;
                        bool foundBP = getProperty(queryBreakpoints, edge, bp_dummy);
                        if (foundBP)
                        {
                            TBreakpoint & oldBp = property(queryBreakpoints, edge);
                            if (cargo < getCargo(edge))// && oldBp.svtype != TBreakpoint::DISPDUPLICATION)
                            {
                                // Check for tandem deletion or insertion
                                // if olBP is deletion and new bp is insertion, and
                                // insertion but (m2.e1 - m1.b1) > (m2.e2 - m1.b1) (equivalent to
                                // deletion but (m2.e1 - m1.b1) < (m2.e2 - m1.b1) ) keep deletion bp, prob. DEL:ME?
                                if (! ( (bp.svtype == 1 && oldBp.svtype == 2) &&
                                        ((*stMatch2).end1 - (*stMatch1).begin1) > ((*stMatch2).end2 - (*stMatch1).begin2)
                                      )
                                   ) // 1=Insertion, 2=Deletion, if both condictions true then tandem del
                                {
                                    // Replace cargo and bp
                                    assignCargo(edge, cargo);
                                    oldBp = bp;
                                }
                            }
                        }
                        else
                        {
                            assignProperty(queryBreakpoints, edge, bp);
                        }

                    }
                    else if (edge == 0)
                    {
                        edge = addEdge(graph, m1, m2, cargo);
                        resizeEdgeMap(queryBreakpoints.slotLookupTable, graph);

                        if (artificialBP)
                            assignProperty(queryBreakpoints, edge);
                        else
                            assignProperty(queryBreakpoints, edge, bp);
                    }
                }
                if (swap)
                {
                    std::swap(stMatch1, stMatch2);
                    m1Begin = (*stMatch1).begin1;
                    m1End = (*stMatch1).end1;
                }
                doBP = false;
                insertEdge = false;
            }
        }
    }
}

// Chain all matches of each query
template <typename TSequence, typename TId, typename TScoreAlloc, typename TMSplazerChain>
void _chainQueryMatches(StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & stellarMatches,
                        String<TScoreAlloc> & distanceScores,
                        String<TMSplazerChain> & queryChains,
                        StringSet<TId> & queryIds,
                        StringSet<TSequence> & queries,
                        String<unsigned> & readJoinPositions,
                        MSplazerOptions const & msplazerOptions)
{

    for (unsigned i = 0; i < length(stellarMatches); ++i)
    {
        TScoreAlloc matchDistanceScores = distanceScores[i];
        TMSplazerChain chain(matchDistanceScores);
        if (msplazerOptions.pairedEndMode)
            chain.mateJoinPosition = readJoinPositions[i];
        // TMSplazerChain chain(matchDistanceScores, readJoinPositions[i]);

        if (length(stellarMatches[i].matches) == 0)
        {
            // Insert single match into graph, no extra edges beside from start and to end
            chain.isEmpty = true;
        }
        else
        {
            // Graph init
            if (msplazerOptions.pairedEndMode)
                _initialiseGraphMatePairs(stellarMatches[i], queryIds[i], chain, msplazerOptions);
            else
                //_initialiseGraphNoBreakend(stellarMatches[i], chain, msplazerOptions);
                _initialiseGraph(stellarMatches[i], queryIds[i], chain, msplazerOptions);

            // Chain compatible matches
            _chainMatches(stellarMatches[i],
                          queryIds[i],
                          queries[i],
                          chain.graph,
                          chain.matchDistanceScores,
                          chain.breakpoints,
                          chain,
                          msplazerOptions);
            // Chain matches comptable according to reference
            _chainMatchesReference(stellarMatches[i],
                                   queryIds[i],
                                   queries[i],
                                   chain.graph,
                                   chain.matchDistanceScores,
                                   chain.breakpoints,
                                   chain,
                                   msplazerOptions);

        }
        // Reevaluate chains with translocations/duplications
        appendValue(queryChains, chain);
    }
}

// Analyze chains in read graph by calling DAG shortest path algorithm
template <typename TMSplazerChain>
void _analyzeChains(String<TMSplazerChain> & queryChains)
{
    InternalPropertyMap<int> weightMap;
    // typedef typename TMSplazerChain::TGraph TGraph;
    // typedef typename Size<TGraph>::Type TGraphSize;

    for (unsigned i = 0; i < length(queryChains); ++i)
    {
        if (!queryChains[i].isEmpty)
        {
            resizeVertexMap(queryChains[i].distMap, queryChains[i].graph);
            dagShortestPath(queryChains[i].predMap,
                            queryChains[i].distMap,
                            queryChains[i].graph,
                            queryChains[i].startVertex,
                            weightMap);
        }
        // else
        // std::cerr << " Chain for query " << i << ", is empty!" << std::endl;//", queryID: " << queryIDs[i] <<
    }
}

template <typename TBreakpoint>
inline unsigned _getSVSize(TBreakpoint const & bp)
{
    return (bp.endSeqPos - bp.startSeqPos);
}

// Function deletionSupport
// tempBP is also deletion, are both part of a transl? extract and save as transl
// if (del1.endSeqPos == del2.startSeqPos) // add position variance of at least 3bp
//      transl(del1.startSeqPos, del2.endSeqPos); transl.dupMidPos = del1.endSeqPos;
// if (del2.endSeqPos == del1.startSeqPos) // add position variance of at least 3bp
//      transl(del2.startSeqPos, del1.endSeqPos); transl.dupMidPos = del2.endSeqPos;
// tempBP is insertion, are both part of a transl/dup? (Only possible when regarding more than 1 chain per
// read
// tempBP is translocation
// Returns true if tempBP is obsolete and will then be erased
template <typename TBreakpoint>
bool _translDelSupport(TBreakpoint & bp, TBreakpoint & tempBP, bool & foundNewBP, unsigned const & bpPosRange)
{
    if (bp.startSeqId != tempBP.startSeqId)
        return false;
    /*
    if ( (bp.endSeqPos - bp.startSeqPos) < 31)
        return false;
    if ( (tempBP.endSeqPos - tempBP.startSeqPos) < 31)
        return false;
        */

    // Note: If there is a deletion supporting a duplication, the position that is equal in both deletion and duplication
    // determines the target position of a duplication.
    // If the deletion is the one that distinguishes the translocation from a duplication, then the deletion does not
    // cover the target position, but we assume that in this case we observe both deletions and if not, it is more likely
    // to observe the deletion that covers the target position.
    // Function returns true if the old bp (tempBP) can be deleted from the old set bc it became obsolete.


    // Case 1: new bp is deletion, old bp (tempBP) is duplication or translocation
    if (bp.svtype == TBreakpoint::DELETION
            && (tempBP.svtype == TBreakpoint::DISPDUPLICATION || tempBP.svtype == TBreakpoint::TRANSLOCATION) )
    {
        bool startPosSameRange = _posInSameRange(bp.startSeqPos, tempBP.startSeqPos, bpPosRange);
        bool endPosSameRange = _posInSameRange(bp.endSeqPos, tempBP.endSeqPos, bpPosRange);
        // Sometimes a deletion is called that has same start and end coordinates as duplication
        // Then we just want to add this deletion as support
        if (startPosSameRange && endPosSameRange)
        {
            appendSupportId(tempBP, bp.supportIds);
            foundNewBP = false;
            return false;
        }
        // deletion start position equals duplication start position
        // If startPos of tempBP is equal to startPos of bp, bp is a pseudodeletion of tempBP with
        // tempBP.middlePos=bp.endPos
        // If there is no second pseudodeletion marking the variant as a translocation, then the equal start
        // position is also the target position of the duplication
        if (startPosSameRange)
        {
            appendSupportId(tempBP, bp.supportIds);
            tempBP.dupTargetPos = tempBP.startSeqPos;
            tempBP.dupMiddlePos = bp.endSeqPos;
            tempBP.midPosId = bp.endSeqId;
            tempBP.translSuppStartPos = true;
            //if (tempBP.translSuppEndPos)
                //tempBP.svtype = TBreakpoint::TRANSLOCATION;
            foundNewBP = false;
            return false;
        }
        // deletion end position equals duplication end position
        // If endPos of tempBP is equal to endPos of bp, bp is a pseudodeletion of tempBP with
        // tempBP.middlePos=bp.startPos
        // If there is no second pseudodeletion marking the variant as a translocation, then the equal end
        // position is also the target position of the duplication
        if (endPosSameRange)
        {
            appendSupportId(tempBP, bp.supportIds);
            tempBP.dupTargetPos = tempBP.endSeqPos;
            tempBP.dupMiddlePos = bp.startSeqPos;
            tempBP.midPosId = bp.startSeqId;
            tempBP.translSuppEndPos = true;
            //if (tempBP.translSuppStartPos)
               // tempBP.svtype = TBreakpoint::TRANSLOCATION;
            foundNewBP = false;
            return false;
        }

        return false;
    }

    // Case 2: new bp is duplication or translocation, old bp (tempBP) is deletion
    if ( (bp.svtype == TBreakpoint::DISPDUPLICATION || bp.svtype == TBreakpoint::TRANSLOCATION)
            && tempBP.svtype == TBreakpoint::DELETION)
    {
        bool startPosSameRange = _posInSameRange(bp.startSeqPos, tempBP.startSeqPos, bpPosRange);
        bool endPosSameRange = _posInSameRange(bp.endSeqPos, tempBP.endSeqPos, bpPosRange);
        // Sometimes a deletion is called that has same start and end coordinates as duplication
        // Then we just want to add this deletion as support
        if (startPosSameRange && endPosSameRange)
        {
            appendSupportId(bp, tempBP.supportIds);
            foundNewBP = true;
            return true;
        }

        // If startPos of bp is equal to startPos of tempBP, tempBP is a pseudodeletion of bp with
        // bp.middlePos=tempBP.endPos
        if (startPosSameRange)
        {
            appendSupportId(bp, tempBP.supportIds);
            bp.dupTargetPos = bp.startSeqPos;
            bp.dupMiddlePos = tempBP.endSeqPos;
            bp.midPosId = tempBP.endSeqId;
            bp.translSuppStartPos = true;
            // tempBP = bp;
            foundNewBP = true;
            return true;
        }
        // If endPos of bp is equal to endPos of tempBP, tempBP is a pseudodeletion of bp with
        // bp.middlePos=tempBP.startPos
        if (endPosSameRange)
        {
            appendSupportId(bp, tempBP.supportIds);
            bp.dupTargetPos = bp.endSeqPos;
            bp.dupMiddlePos = tempBP.startSeqPos;
            bp.midPosId = tempBP.startSeqId;
            bp.translSuppEndPos = true;
            // tempBP = bp;
            foundNewBP = true;
            return true;
        }
        return false;
    }

    // Case 3: Two possible pseudodeletions of a translocation events
    if (bp.svtype == TBreakpoint::DELETION && tempBP.svtype == TBreakpoint::DELETION)
    {
        // If startPos of bp is equal to endPos of tempBP, bp is upstream of tempBP, they form the translocation
        // startPos=tempBP.startPos, middlePos=min(tempBP.endPos,bp.startPos), endPos=bp.endPos
        if (_posInSameRange(bp.startSeqPos, tempBP.endSeqPos, bpPosRange))
        {
            TBreakpoint translBP(tempBP.startSeqId,
                               bp.endSeqId,
                               tempBP.startSeqStrand,
                               bp.endSeqStrand,
                               tempBP.startSeqPos,
                               bp.endSeqPos,
                               tempBP.readStartPos,
                               bp.readEndPos
                              );
            translBP.svtype = TBreakpoint::TRANSLOCATION;
            translBP.dupMiddlePos = std::min(tempBP.endSeqPos,bp.startSeqPos);
            translBP.midPosId = tempBP.endSeqId;
            translBP.translSuppStartPos = true;
            translBP.translSuppEndPos = true;
            appendSupportId(translBP, bp.supportIds);
            appendSupportId(translBP, tempBP.supportIds);
            bp = translBP;
            foundNewBP = true;
            return true;
        }
        // If startPos of tempBP is equal to endPos of bp, bp is downstream of tempBP, they form the translocation
        // startPos=bp.startPos, middlePos=min(bp.endPos,tempBP.startPos), endPos=tempBP.endPos
        if (_posInSameRange(bp.endSeqPos, tempBP.startSeqPos, bpPosRange))
        {
            TBreakpoint translBP(bp.startSeqId,
                               tempBP.endSeqId,
                               bp.startSeqStrand,
                               tempBP.endSeqStrand,
                               bp.startSeqPos,
                               tempBP.endSeqPos,
                               bp.readStartPos,
                               tempBP.readEndPos
                              );
            translBP.svtype = TBreakpoint::TRANSLOCATION;
            translBP.dupMiddlePos = std::min(bp.endSeqPos, tempBP.startSeqPos);
            translBP.midPosId = bp.endSeqId;
            translBP.translSuppStartPos = true;
            translBP.translSuppEndPos = true;
            appendSupportId(translBP, bp.supportIds);
            appendSupportId(translBP, tempBP.supportIds);
            bp = translBP;
            foundNewBP = true;
            return true;
        }
        return false;
    }
    return false;
}

// Returns true if a 3rd breakpoint (newDelBP) is inferred.
template <typename TBreakpoint>
inline bool _invDelClassification(TBreakpoint & bp, TBreakpoint & tempBP, TBreakpoint & newInvBP, unsigned const & bpPosRange)
{
    bool newInv = false;
    if (bp.svtype != TBreakpoint::INVERSION || tempBP.svtype != TBreakpoint::INVERSION)
        return false;
    if (bp.startSeqId != tempBP.startSeqId)
        return false;
    if (tempBP.inferredBP)
        return false;

    bool startPosSameRange = _posInSameRange(bp.startSeqPos, tempBP.startSeqPos, bpPosRange);
    bool endPosSameRange = _posInSameRange(bp.endSeqPos, tempBP.endSeqPos, bpPosRange);

    // Case 1: new bp is outer inv inv(v,z), old (tempBP) is inner inv inv(w,z)||inv(v,w)
    // inv(v,w)
    if (startPosSameRange && (tempBP.endSeqPos < bp.endSeqPos))
    {
        // Change new bp from inv(v,z) to deletion del(w,z)
        bp.svtype = TBreakpoint::DELETION;
        bp.startSeqPos = tempBP.endSeqPos;
        bp.startSeqStrand = bp.endSeqStrand;
        bp.inferredBP = true;
        appendSupportId(bp, tempBP.supportIds);
        return false;
    }
    // inv(w,z)
    if (endPosSameRange && (bp.startSeqPos < tempBP.startSeqPos))
    {
        // Change new bp from inv(v,z) to deletion del(v,w)
        bp.svtype = TBreakpoint::DELETION;
        bp.endSeqPos = tempBP.startSeqPos;
        bp.endSeqStrand = bp.startSeqStrand;
        bp.inferredBP = true;
        appendSupportId(bp, tempBP.supportIds);
        return false;
    }
    // Case 2: new bp is inner inv inv(w,z)||inv(v,w), old (tempBP) is outer inv inv(v,z)
    // inv(v,w)
    if (startPosSameRange && (bp.endSeqPos < tempBP.endSeqPos))
    {
        // Change old bp from inv(v,z) to deletion del(w,z)
        tempBP.svtype = TBreakpoint::DELETION;
        tempBP.startSeqPos = bp.endSeqPos;
        tempBP.startSeqStrand = tempBP.endSeqStrand;
        tempBP.inferredBP = true;
        appendSupportId(tempBP, bp.supportIds);
        return false;
    }
    // inv(w,z)
    if (endPosSameRange && (bp.startSeqPos < tempBP.startSeqPos))
    {
        // Change old bp from inv(v,z) to deletion del(w,z)
        tempBP.svtype = TBreakpoint::DELETION;
        tempBP.endSeqPos = bp.startSeqPos;
        tempBP.endSeqStrand = tempBP.startSeqStrand;
        tempBP.inferredBP = true;
        appendSupportId(tempBP, bp.supportIds);
        return false;
    }
    // Case 3: new bp inv(v,w) and old bp inv(u,z) are overlapping, i.e. there are 2 adjacent dels
    // Either i) v < u < w < z, then bp becomes del(v,u), tempBP del(w,z), newInvBP inv(u,w)
    // or     ii)u < v < z < w, then bp becomes del(z,w), tempBP del(u,v), newInvBP inv(v,z)
    // i)
    if (_checkMatchOverlap(bp.startSeqPos, bp.endSeqPos, tempBP.startSeqPos, tempBP.endSeqPos))
    {
        newInvBP.svtype = TBreakpoint::INVERSION;
        newInvBP.startSeqPos = tempBP.startSeqPos;
        newInvBP.endSeqPos = bp.endSeqPos;
        newInvBP.startSeqStrand = tempBP.startSeqStrand;
        newInvBP.endSeqStrand = bp.endSeqStrand;
        newInvBP.inferredBP = true;
        bp.svtype = TBreakpoint::DELETION;
        bp.endSeqPos = newInvBP.startSeqPos;
        bp.endSeqStrand = bp.startSeqStrand;
        bp.inferredBP = true;
        tempBP.svtype = TBreakpoint::DELETION;
        tempBP.startSeqPos = newInvBP.endSeqPos;
        tempBP.startSeqStrand = tempBP.endSeqStrand;
        tempBP.inferredBP = true;
        appendSupportId(tempBP, bp.supportIds);
        newInvBP.supportIds = tempBP.supportIds;
        newInvBP.support = tempBP.support;
        bp.supportIds = tempBP.supportIds;
        bp.support = tempBP.support;
        return true;
    }
    // ii)
    if (_checkMatchOverlap(tempBP.startSeqPos, tempBP.endSeqPos, bp.startSeqPos, bp.endSeqPos))
    {
        newInvBP.svtype = TBreakpoint::INVERSION;
        newInvBP.startSeqPos = bp.startSeqPos;
        newInvBP.endSeqPos = tempBP.endSeqPos;
        newInvBP.startSeqStrand = bp.startSeqStrand;
        newInvBP.endSeqStrand = tempBP.endSeqStrand;
        newInvBP.inferredBP = true;
        bp.svtype = TBreakpoint::DELETION;
        bp.startSeqPos = newInvBP.endSeqPos;
        bp.startSeqStrand = bp.endSeqStrand;
        bp.inferredBP = true;
        tempBP.svtype = TBreakpoint::DELETION;
        tempBP.endSeqPos = newInvBP.startSeqPos;
        tempBP.endSeqStrand = tempBP.startSeqStrand;
        tempBP.inferredBP = true;
        appendSupportId(tempBP, bp.supportIds);
        newInvBP.supportIds = tempBP.supportIds;
        newInvBP.support = tempBP.support;
        bp.supportIds = tempBP.supportIds;
        bp.support = tempBP.support;
        return true;
    }

    return newInv;
}

template <typename TBreakpoint>
inline bool _insertBreakend(String<TBreakpoint> & countedBE, TBreakpoint & be, unsigned const & bpPosRange)
{
    bool newBE = true;
    // Compare new breakend to breakends already in the set and add to set if new or supportID if not
    for (unsigned i = 0; i < length(countedBE); ++i)
    {
        TBreakpoint & tempBE = countedBE[i];
        // if .breakend is equal, both breakends are both either left or right breakends
        if (_similarBreakends(be, tempBE, bpPosRange) && tempBE.breakend == be.breakend)
        {
            appendSupportId(tempBE, be.supportIds);
            newBE = false;
            // Take leftmost (left BE,0) or rightmost (right BE,1) position
            if ((tempBE.breakend && tempBE.startSeqPos > be.startSeqPos)
                    || (!tempBE.breakend && tempBE.startSeqPos < be.startSeqPos))
            {
                tempBE.startSeqPos = be.startSeqPos;
                tempBE.endSeqPos = be.endSeqPos;
            }
        }
    }
    if (newBE)
        appendValue(countedBE, be);
    return newBE;
}

template <typename TBreakpoint>
inline bool _addBeSupport(String<TBreakpoint> & countedBP, TBreakpoint & be, unsigned const & bpPosRange)
{
    bool newBE = true;
    // Compare new breakend to breakpoints already in the set and add supportIDs if similar
    for (unsigned i = 0; i < length(countedBP); ++i)
    {
        TBreakpoint & tempBP = countedBP[i];
        if (_breakendSupport(be, tempBP, bpPosRange))
        {
            appendSupportId(tempBP, be.supportIds);
            newBE = false;
            // return false;
        }
    }
    return newBE;
    // return true;
}

template <typename TBreakpoint>
inline void _getBeSupport(String<TBreakpoint> & countedBE, TBreakpoint & bp, unsigned const & bpPosRange)
{
    String<unsigned> toBeErased;
    // Compare new breakpoint to breakends in the set, add supportIDs and erase breakend if similar
    for (unsigned i = 0; i < length(countedBE); ++i)
    {
        TBreakpoint & tempBE = countedBE[i];
        if (_breakendSupport(tempBE, bp, bpPosRange))
        {
            // Append support Ids of old bp (tempBP) to new bp
            appendSupportId(bp, tempBE.supportIds);
            appendValue(toBeErased, i);
        }
    }
    for (unsigned i = 0; i < length(toBeErased); ++i)
    {
        unsigned index = length(toBeErased) - 1 - i;
        erase(countedBE, toBeErased[index]);
    }
}

// Combines simple BPs (deletions, inversions) w/ or to complex SVs (translocations, duplications, combined inv-del)
template <typename TBreakpoint>
inline void _inferComplexBP(String<TBreakpoint> & globalBreakpoints, TBreakpoint & bp, unsigned const & bpPosRange, unsigned & similarBPId)
{
    String<unsigned> toBeErased;
    bool keepBP = true;
    bool newInv = false;
    TBreakpoint newInvBP = bp;
    // Compare new breakpoint to breakpoints in the set, add supportIDs if similar or create new, more complex
    // breakpoint if supporting breakpoint
    for (unsigned i = 0; i < length(globalBreakpoints); ++i)
    {
        TBreakpoint & tempBP = globalBreakpoints[i];
        // Special case: one of the breakpoints is a deletion, the other a duplication or translocation, then the del
        // can either be part of the duplication and should not be in the output list or it can distinguish the dup
        // from a translocation in which case the svtype is updated from duplication to translocation
        // Note: both deletions cannot be distinguished for which it is, so only if both have been witnessed we can
        // update from duplication to translocation
        // Note: if both are oberserved at the same time without a duplication/translocation, we infer the presence of
        // a translocation and insert a new translocation breakpoint instead of both deletions

        if (_translDelSupport(bp, tempBP, keepBP, bpPosRange))
            appendValue(toBeErased, i);

        // Special case: both breakpoints are nested inversions. Then, it is more likely that there is actually a
        // combination of an inversion with 1 (or 2) adjacent deletion(s).
        // If one includes the other, we keep the inner inversion and infer an deletion of the remaining region.
        // If both overlap, then the overlapping region is the inversion, the remaining regions are inferred dels.
        // In this case, we need to create a third bp
        else if (_invDelClassification(bp, tempBP, newInvBP, bpPosRange))
            newInv = true;

    }
    // Erase obsolete bps in descending order
    for (unsigned i = 0; i < length(toBeErased); ++i)
    {
        unsigned index = length(toBeErased) - 1 - i;
        erase(globalBreakpoints, toBeErased[index]);
    }
    // Append breakpoint if new
    if (keepBP)
    {
        if (bp.similar == std::numeric_limits<unsigned>::max())
        {
            bp.similar = similarBPId;
            ++similarBPId;
        }
        appendValue(globalBreakpoints, bp);
    }
    // Append new deletion breakpoint if appl./insert new deletion bc. there could be one already in the set?
    if (newInv)
        appendValue(globalBreakpoints, newInvBP);
        // _insertBreakpoint(countedBP, newDelbp, bpPosRange, similarBPId++?);
}

// Insert Breakpoint into string of breakpoints if it is not already in the set. Returns true if breakpoint was new and
// has been inserted or false if breakpoint was already in the set (and just has been counted).
// For insertions, also the insertion length has to be the same
template <typename TBreakpoint>
inline void _insertBreakpoint(String<TBreakpoint> & countedBP, TBreakpoint & bp, unsigned const & bpPosRange, unsigned & similarBPId)
{
    bool newBP = true;
    // Compare new breakpoint to breakpoints in the set, add supportIDs if similar or create new, more complex
    // breakpoint if supporting breakpoint
    for (unsigned i = 0; i < length(countedBP); ++i)
    {
        TBreakpoint & tempBP = countedBP[i];
        if (bp == tempBP)
        {
            appendSupportId(tempBP, bp.supportIds);
            // Very special case, should be no other support likely
            newBP = false;
        }
        else if (_similarBreakpoints(bp, tempBP, bpPosRange))
        {
            appendSupportId(tempBP, bp.supportIds);
	    if (bp.dupMiddlePos != std::numeric_limits<unsigned>::max() && tempBP.dupMiddlePos != std::numeric_limits<unsigned>::max()
		&& bp.dupMiddlePos < tempBP.dupMiddlePos)
            {
                tempBP.dupMiddlePos = bp.dupMiddlePos; 
	    }
	    if (bp.startSeqPos < tempBP.startSeqPos)
		tempBP.startSeqPos = bp.startSeqPos;
	    if (bp.endSeqPos < tempBP.endSeqPos)
		tempBP.endSeqPos = bp.endSeqPos;
            // Very special case, should be no other support likely
            newBP = false;
            //bp.similar = tempBP.similar;
            // Mark breakpoint to be similar via an ID
        }

    }
    // Append breakpoint if new
    if (newBP)
    {
        if (bp.similar == std::numeric_limits<unsigned>::max())
        {
            bp.similar = similarBPId;
            ++similarBPId;
        }
        appendValue(countedBP, bp);
    }
}

// Wrapper insert breakpoint function with separated breakends, calls insert breakpoint fct below
template <typename TBreakpoint>
void _insertBreakpoints(String<TBreakpoint> & countedBP,
                        String<TBreakpoint> & countedBE,
                        String<TBreakpoint> & newBP,
                        MSplazerOptions const & options,
                        unsigned & similarBPId)
{
    for (unsigned i = 0; i < length(newBP); ++i)
    {
        TBreakpoint & bp = newBP[i];
        // If new bp is a breakend, see if it is supporting any breakpoint already in the set and if not, insert into set of
        // breakends
        // If new bp is a breakpoint, see if there is a supporting breakend in the breakend set and insert into set of
        // breakpoints
        if (bp.svtype == TBreakpoint::BREAKEND)
        {
            // If breakend is not supporting any real breakpoint, insert into breakend set
            if (_addBeSupport(countedBP, bp, options.breakpointPosRange))
            {
                _insertBreakend(countedBE, bp, options.breakpointPosRange);
            }
        }
        else
        {
            _getBeSupport(countedBE, bp, options.breakpointPosRange);
            _insertBreakpoint(countedBP, bp, options.breakpointPosRange, similarBPId);
        }
    }
}

// Insert breakpoint function w/o breakends
template <typename TBreakpoint>
void _insertBreakpoints(String<TBreakpoint> & countedBP, String<TBreakpoint> & newBP, MSplazerOptions const & options)
{
    for (unsigned i = 0; i < length(newBP); ++i)
        _insertBreakpoint(countedBP, newBP[i], options.breakpointPosRange);
}

//TODO(ktrappe): Trimming functionality is disabled atm and needs adaption to new alignment properties and moduls
template <typename TMatch, typename TPos>
void _trimMatches(TMatch & fstMatch, TMatch & sndMatch, String<TPos> & splitPos)
{
    // Note: The "end" of the first match has to be trimmed with the "begin" position of the breakpoint,
    // and the "begin" of the second match with the "end" position of the breakpoint
    // 0:referenceStart, 1:referenceEnd, 2:readStart, 3:readEnd
    TPos & refEndSplitPos = splitPos[0];
    TPos & refBeginSplitPos = splitPos[1];
    TPos & readEndSplitPos = splitPos[2];
    TPos & readBeginSplitPos = splitPos[3];

    SEQAN_ASSERT_EQ(length(fstMatch.row1), length(fstMatch.row2));
    _trimMatchEnd(fstMatch, refEndSplitPos, readEndSplitPos);
    SEQAN_ASSERT_EQ(length(fstMatch.row1), length(fstMatch.row2));
    SEQAN_ASSERT_EQ(length(sndMatch.row1), length(sndMatch.row2));
    _trimMatchBegin(sndMatch, refBeginSplitPos, readBeginSplitPos);
    SEQAN_ASSERT_EQ(length(sndMatch.row1), length(sndMatch.row2));
}

template <typename TMatch, typename TPos>
void _trimMatches(String<TMatch> & matchChain, StringSet<String<TPos> > & splitPos)
{
    for (unsigned i = 0; i < length(splitPos); ++i)
    {
        // Only trim matches that overlap in either reference or read sequence
        if (_checkMatchOverlap(matchChain[i + 1].begin1, matchChain[i + 1].end1, matchChain[i].begin1,
                               matchChain[i].end1)
           || _checkMatchOverlap(matchChain[i + 1].begin2, matchChain[i + 1].end2, matchChain[i].begin2,
                                 matchChain[i].end2))
            _trimMatches(matchChain[i + 1], matchChain[i], splitPos[i]);
    }
}

// TODO(ktrappe): Indel extraction is disabled atm cause it depends on match trimming
template <typename TBreakpoint, typename TSequence, typename TId, typename TMatch>
void _getChainIndels(String<String<TMatch> > & bestChains,
                     String<TBreakpoint> & globalStellarIndels,
                     TId const & queryId,
                     TSequence & query)
{
    typedef typename Iterator<String<TMatch> >::Type TIterator;
    for (unsigned i = 0; i < length(bestChains); ++i)
    {
        TIterator itStellarMatches = begin(bestChains[i]);
        TIterator itEndStellarMatches = end(bestChains[i]);
        for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
            _getStellarIndel(*itStellarMatches, globalStellarIndels, queryId, query);
    }
}

template <typename TMSplazerChain>
void _findPartialChains(TMSplazerChain & queryChain)
{
    if (outDegree(queryChain.graph, queryChain.startVertex) > 0 && inDegree(queryChain.graph, queryChain.endVertex) > 0)
        return;

    if (inDegree(queryChain.graph, queryChain.endVertex) > 0)
    {
        queryChain.isPartial = true;
        return;
    }

    if (outDegree(queryChain.graph, queryChain.startVertex) > 0)
    {
        queryChain.isPartial = true;
        return;
    }
    return;
}

// Finding the best chain (belonging to the shortest path) and reporting the breakpoints, if any.
template <typename TBreakpoint, typename TMSplazerChain, typename TMatch>
bool _findBestChain(TMSplazerChain & queryChain, String<TMatch> & stellarMatches,
                    String<TBreakpoint> & tmpGlobalBreakpoints,
                    // MSplazerOptions const & msplazerOptions,
                    unsigned & bcc)
{
    typedef typename TMSplazerChain::TVertexDescriptor TVertexDescriptor;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    typedef typename TBreakpoint::TPos TPos;
    bool foundBP = false;
    TVertexDescriptor spVertex1, spVertex2;
    TEdgeDescriptor edge;

    String<TMatch> tmpBestChain;
    StringSet<String<TPos> > tmpSplitPos;

    if (queryChain.isEmpty)
        return 0;

    spVertex1 = queryChain.endVertex;
    while (spVertex1 != queryChain.startVertex)
    {
        // Getting vertex descriptor of anchestor
        spVertex2 = getProperty(queryChain.predMap, spVertex1);
        // if vertex descriptor has max value then there is no (shortest) path from start to end.
        if (spVertex2 == std::numeric_limits<TVertexDescriptor>::max())
        {
            // queryChain.isPartial = true;
            ++bcc;
            return 0;
        }
        // Append tmpStellarMatch
        // std::cerr << "ID: " << _getId(spVertex2) << std::endl;
        if (_getId(spVertex2) != queryChain.startVertex)
        {
            TMatch stMatch = stellarMatches[_getId(spVertex2)];
            appendValue(tmpBestChain, stMatch);
        }
        // Getting edge between both vertices; findEdge returns 0 if there is no such edge
        edge = findEdge(queryChain.graph, spVertex2, spVertex1);
        // If you ever want just a copy of the breakpoints in the global set, redefine getProperty to give back
        // getValue() instead of value(), keep in mind that bp has to be defined somewhere else or it gets lost
        TBreakpoint bp;

        // Get breakpoint on edge
        if (edge != 0)
            foundBP = getProperty(queryChain.breakpoints, edge, bp);

        // Append breakpoint to tmpGlobalBreakpoints
        if (foundBP)
        {
            // insert BP into tmpGlobalBreakpoints
            appendValue(tmpGlobalBreakpoints, bp);
            foundBP = false;
            // Append splitPos, which are the source position of the breakoint in read and reference
            // referenceStart, referenceEnd, readStart, readEnd
            String<TPos> splitPos;
            resize(splitPos, 4);

            // Store breakpoint pos., mind matches of different order for translocations and reverse strand deletions
            // if ((bp.svtype == 6 && bp.startSeqId == bp.endSeqId) || bp.revStrandDel) // 6=translocation
            // TODO(ktrappe): Adjust SVType
            if ((bp.svtype == TBreakpoint::TRANSLOCATION) || bp.revStrandDel)
            {
                splitPos[0] = bp.endSeqPos;
                splitPos[1] = bp.startSeqPos;
                splitPos[2] = bp.readStartPos;
                splitPos[3] = bp.readEndPos;
            }
            else
            {
                splitPos[0] = bp.startSeqPos;
                splitPos[1] = bp.endSeqPos;
                splitPos[2] = bp.readStartPos;
                splitPos[3] = bp.readEndPos;
            }
            appendValue(tmpSplitPos, splitPos);
        }
        spVertex1 = spVertex2;
    }

    // only if chain was complete
    // trim matches with split pos
    // if (msplazerOptions.simThresh > 0)
    // _trimMatches(tmpBestChain, tmpSplitPos);
    // Adjust bp positions after match trimming for

    // insert bestChain into queryChain object
    insertBestChain(queryChain, tmpBestChain);
    return 1;
}

// Finding the best chain (belonging to the shortest path) and reporting the breakpoints, if any.
template <typename TMSplazerChain, typename TBreakpoint, typename TQueryMatches>
// , typename TSequence, typename TId>
void _findAllBestChains(String<TMSplazerChain> & queryChains,
                        StringSet<TQueryMatches> & queryMatches,
                        // StringSet<TSequence> & queries,
                        // StringSet<TId> const & queryIds,
                        String<TBreakpoint> & globalBreakpoints,
                        // String<TBreakpoint> & globalStellarIndels,
                        MSplazerOptions const & msplazerOptions
                        )
{
    /*
    String<unsigned> chainSizeCount;
    resize(chainSizeCount, 10);
    for(unsigned i = 0; i < length(chainSizeCount); ++i)
        chainSizeCount[i] = 0;
    */
    String<TBreakpoint> globalBreakends;
    String<TBreakpoint> tmpGlobalBreakpoints;
    unsigned brokenChainCount = 0;
    unsigned similarBPId = 0;
    // 1st pass scanning:
    // Extract 
    for (unsigned i = 0; i < length(queryChains); ++i)
    {
        String<TBreakpoint> localBreakpoints;
        // String<TBreakpoint> tmpStellarIndels;
        if (_findBestChain(queryChains[i], queryMatches[i].matches, localBreakpoints, // msplazerOptions,
                           brokenChainCount))
        {
            _insertBreakpoints(tmpGlobalBreakpoints, globalBreakends, localBreakpoints, msplazerOptions, similarBPId);
            // get small indels from matches
            // _getChainIndels(queryChains[i].bestChains, globalStellarIndels, queryIds[i], queries[i]);
            /*
            for(unsigned j = 0; j < length(queryChains[i].bestChains); ++j)
                ++chainSizeCount[length(queryChains[i].bestChains[j])];
            */
        }
    }
    // 2nd pass scanning:
    // Only keep BPs w/ enough support and check if these can be combined to more complex ones
    for (unsigned i = 0; i < length(tmpGlobalBreakpoints); ++i)
    {
        if (tmpGlobalBreakpoints[i].support >= msplazerOptions.support)
        {
            if (msplazerOptions.inferComplexBP)
                _inferComplexBP(globalBreakpoints, tmpGlobalBreakpoints[i], msplazerOptions.breakpointPosRange, similarBPId);
            else
                appendValue(globalBreakpoints, tmpGlobalBreakpoints[i]);
        }
    }
    append(globalBreakpoints, globalBreakends);
}

// TODO(ktrappe): Finding suboptimal chains using Sandros enumerateAcyclicPaths function (disabled atm)
// template <typename TBreakpoint, typename TMSplazerChain, typename TMatch>
template <typename TMSplazerChain>
bool _findSuboptimalChains(TMSplazerChain & queryChain // ,
                           // String<TMatch> & stellarMatches,
                           // String<TBreakpoint> & tmpGlobalBreakpoints,
                           // MSplazerOptions const & msplazerOptions
                           )
{
    typedef typename TMSplazerChain::TVertexDescriptor TVertexDescriptor;
    typedef typename TMSplazerChain::TGraph TGraph;
    typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
    // typedef typename TBreakpoint::TPos TPos;
    // typedef typename TBreakpoint::TId TId;

    typedef::std::multimap<int, unsigned> TChainScoreMap;
    TChainScoreMap chainScoreMap;

    // Get all paths from start to end vertex of the graph
    StringSet<String<TVertexDescriptor> > all_paths;
    // enumerateAcyclicPaths(queryChain.graph, queryChain.startVertex, queryChain.endVertex, all_paths);

    // Get edit distance of all paths (cargo on edges) and store path id and edit distance in a sorted list/map
    for (unsigned i = 0; i < length(all_paths); ++i)
    {
        String<TVertexDescriptor> & path = all_paths[i];
        int pathScore = 0;
        // Sum up edge weights along the path
        TVertexDescriptor spVertex1, spVertex2;
        TEdgeDescriptor edge;
        for (unsigned j = 0; j < length(path) - 1; ++j)
        {
            spVertex1 = path[j];
            spVertex2 = path[j + 1];
            edge = findEdge(queryChain.graph, spVertex1, spVertex2);
            pathScore += getCargo(edge);

        }
        // Insert tuple of score and position into chainScoreMap
        chainScoreMap.insert(std::make_pair(pathScore, i));
    }
    std::cout << "Number of paths: " << length(all_paths) << std::endl;
    std::multimap<int, unsigned>::iterator it = chainScoreMap.begin();
    // Print best * chains to std out
    std::cout << "paths" << std::endl;
    unsigned count = 0;
    while (it != chainScoreMap.end() && count < 20)
    {
        std::cout << "score id" << std::endl;
        std::cout << (*it).first << " " << (*it).second << std::endl;
        std::cout << "path: " << std::endl;
        String<TVertexDescriptor> & path = all_paths[(*it).second];
        for (unsigned i = 0; i < length(path); ++i)
        {
            std::cout << path[i] << " ";
        }
        std::cout << std::endl;
        ++it;
        ++count;
    }


    // Get breakpoints on the edges and put them into the tmpBreakpoint set

    /*
    TVertexDescriptor spVertex1, spVertex2;
    TEdgeDescriptor edge;

    String<TMatch> tmpBestChain;
    StringSet<String<TPos> > tmpSplitPos;

    if (queryChain.isEmpty)
        return 0;

    spVertex1 = queryChain.endVertex;
    //std::cerr << "Length of stellarMatches: " << length(stellarMatches) << std::endl;
    while (spVertex1 != queryChain.startVertex)
    {
        //Getting vertex descriptor of anchestor
        spVertex2 = getProperty(queryChain.predMap, spVertex1);
        //if vertex descriptor has max value then there is no (shortest) path from start to end.
        if (spVertex2 == maxValue<TVertexDescriptor>())
        {
            //queryChain.isPartial = true;
            ++bcc;
            return 0;
        }
        // Append tmpStellarMatch
        //std::cerr << "ID: " << _getId(spVertex2) << std::endl;
        if (_getId(spVertex2) != queryChain.startVertex)
        {
            TMatch stMatch = stellarMatches[_getId(spVertex2)];
            appendValue(tmpBestChain, stMatch);
        }
        // Getting edge between both vertices
        edge = findEdge(queryChain.graph, spVertex2, spVertex1);
        // If you ever want just a copy of the breakpoints in the global set, redefine getProperty to give back
        // getValue() instead of value(), keep in mind that bp has to be defined somewhere else or it gets lost
        TBreakpoint bp;

        // Get breakpoint on edge
        if (edge != 0)
            foundBP = getProperty(queryChain.breakpoints, edge, bp);
        //else
        //std::cerr << " Edge does not exist, which makes no sense, since you are on the shortest path!" << std::endl;

        // Append breakpoint to tmpGlobalBreakpoints
        if (foundBP)
        {
            // insert BP into tmpGlobalBreakpoints
            appendValue(tmpGlobalBreakpoints, bp);
            foundBP = false;
            // Append splitPos
            String<TPos> splitPos;
            resize(splitPos, 4);

            //std::cerr << bp << std::endl;

            // Store breakpoint pos., mind matches of different order for translocations and reverse strand deletions
            if ((getSVType(bp) == static_cast<TId>("translocation") && bp.startSeqId == bp.endSeqId) || bp.revStrandDel)
            {
                splitPos[0] = bp.endSeqPos;
                splitPos[1] = bp.startSeqPos;
                splitPos[2] = bp.readStartPos;
                splitPos[3] = bp.readEndPos;
            }
            else
            {
                splitPos[0] = bp.startSeqPos;
                splitPos[1] = bp.endSeqPos;
                splitPos[2] = bp.readStartPos;
                splitPos[3] = bp.readEndPos;
            }
            appendValue(tmpSplitPos, splitPos);
        }
        spVertex1 = spVertex2;
    }

    // only if chain was complete
    // trim matches with split pos
    if (msplazerOptions.simThresh > 0)
        _trimMatches(tmpBestChain, tmpSplitPos);

    // insert bestChain into queryChain object
    insertBestChain(queryChain, tmpBestChain);
    */
    return 1;
}

// Finding the best chain (belonging to the shortest path) and reporting the breakpoints, if any.
// template <typename TMSplazerChain, typename TBreakpoint, typename TQueryMatches, typename TSequence, typename TId>
// Note: This is just a place holder atm for future extension of suboptimal path extraction.
template <typename TMSplazerChain>
void _findAllChains(String<TMSplazerChain> & queryChains // ,
                    //    StringSet<TQueryMatches> & queryMatches,
                    //    StringSet<TSequence> & queries,
                    //    StringSet<TId> const & queryIds,
                    //    String<TBreakpoint> & globalBreakpoints,
                    //    String<TBreakpoint> & globalStellarIndels,
                    //    MSplazerOptions const & msplazerOptions
                    )
{
    for (unsigned i = 0; i < length(queryChains); ++i)
    {
        // String<TBreakpoint> tmpGlobalBreakpoints;
        // String<TBreakpoint> tmpStellarIndels;
        // if (_findSuboptimalChains(queryChains[i], queryMatches[i].matches, tmpGlobalBreakpoints, msplazerOptions))
        if (_findSuboptimalChains(queryChains[i]))
        {}
    }
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_MSPLAZER_ALGORITHMS_H_
