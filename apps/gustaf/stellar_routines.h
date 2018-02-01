// ==========================================================================
//                                   Gustaf
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

#ifndef SEQAN_APPS_GUSTAF_STELLAR_ROUTINES_H_
#define SEQAN_APPS_GUSTAF_STELLAR_ROUTINES_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "../stellar/stellar.h"
#include "create_stellarmatches_from_file.h"
using namespace seqan;

// Compare StellarMatches using in order of priority (1) read begin, (2) read end position,
// (3) chromosome begin, (4) chromosome end pos
// Note: Assumes matches derive from the same read/contig
template <typename TSequence, typename TId>
struct CompareStellarMatches
{
    bool operator()(StellarMatch<TSequence, TId> const & match1, StellarMatch<TSequence, TId> const & match2) const
    {
        if (match1.begin2 != match2.begin2)
            return match1.begin2 < match2.begin2;

        if (match1.end2 != match2.end2)
            return match1.end2 < match2.end2;

        if (match1.begin1 != match2.begin1)
            return match1.begin1 < match2.begin1;

        return match1.end1 < match2.end1;
    }

};

// ----------------------------------------------------------------------------
// Function copyGaps()
// ----------------------------------------------------------------------------

template <typename TSeqTo, typename TSeqFrom, typename TSpec>
inline void
copyGaps(Gaps<TSeqTo, TSpec> & toGaps, Gaps<TSeqFrom, TSpec> & fromGaps)
{
    typedef Gaps<TSeqFrom, TSpec>     TGapsFrom;
    typedef typename Iterator<TGapsFrom>::Type TIterFrom;
    typedef Gaps<TSeqTo, TSpec>     TGapsTo;
    typedef typename Iterator<TGapsTo>::Type TIterTo;

    TIterFrom itF = begin(fromGaps);
    TIterFrom itFE = end(fromGaps);
    TIterTo itT = begin(toGaps);

    for (; itF != itFE; ++itF, ++itT)
        if (isGap(itF))
            insertGap(itT);
}

// /////////////////////////////////////////////////////////////////////////////
// Wrapper functions to get Stellar matches

// Function to convert format of begin and end of Stellar matches on the reverse to Seqan convention.
// In Stellar matches: begin < end, additional orientation information needed.
// Seqan: begin > end
template <typename TSequence, typename TID, typename TMatch>
void _convertStellarReverseMatches(TSequence & database, TID & databaseID,
                                   StringSet<QueryMatches<TMatch> > & stellarMatches,
                                   StringSet<QueryMatches<TMatch> > & revStellarMatches)
{
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename StellarMatch<TSequence, TID>::TAlign TAlign;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Size<StringSet<TSequence> >::Type TSize;
    typedef typename Iterator<String<TMatch> >::Type TIter;
    typedef typename Position<TMatch>::Type TPos;

    TSize db_length = length(database);

    // Iterating over matches with Stellar reverse format
    // StringSet is sorted by query
    for (unsigned i = 0; i < length(revStellarMatches); ++i)
    {
        TIter stM = begin(revStellarMatches[i].matches);
        TIter stME = end(revStellarMatches[i].matches);

        for (; stM != stME; goNext(stM))
        {
            // Convert begin and end position and extract the right infix
            TPos begin1 = db_length - (*stM).end1;
            TPos end1 = db_length - (*stM).begin1;
            TInfix dbInf(database, begin1, end1);

            // Create new align object for as input for new Stellar object
            TAlign localAlign;
            resize(rows(localAlign), 1);
            setSource(row(localAlign, 0), host(dbInf));
            appendValue(rows(localAlign), (*stM).row2);

            TRow & row1 = row(localAlign, 0);

            setBeginPosition(row1, begin1);
            setEndPosition(row1, end1);

            copyGaps(row1, (*stM).row1);

            TMatch match(localAlign, databaseID, (*stM).orientation);
            appendValue(stellarMatches[i].matches, match);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _getStellarMatches()
// ----------------------------------------------------------------------------

// Compute Stellar matches using stellar function
// Note: requires to convert matches on reverse strand: StellarMatches of the reverse Strand are being modified,
//  in the sense that they correspond to the right positions within the forward strand
template <typename TSequence, typename TMatches>
void _getStellarMatches(StringSet<TSequence> & queries, StringSet<TSequence> & databases,
                        StringSet<CharString> & databaseIDs, StellarOptions & stellarOptions, TMatches & stellarMatches)
{
    // Finder
    typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;
    // typedef Finder<TSequence, Swift<SwiftLocal> > TFinder;

    // Using Stellars structure for queries
    typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndex;
    // typedef Index<StringSet<TSequence, Dependent<> >, IndexQGram<SimpleShape> > TQGramIndex;
    TQGramIndex qgramIndex(queries);
    resize(indexShape(qgramIndex), stellarOptions.qGram);

    // pattern
    /*
    //Using FragmenStore
    typedef Index<StringSet<TSequence, Owner<ConcatDirect<> > >, IndexQGram<SimpleShape, OpenAddressing> > TQGramIndexFrSt;
    TQGramIndexFrSt qgramIndexFrSt(fragments.readSeqStore);
    stellarOptions.qGram = 4;
    resize(indexShape(qgramIndexFrSt), stellarOptions.qGram);
    typedef Pattern<TQGramIndexFrSt, Swift<SwiftLocal> > TPatternFrSt;
    TPatternFrSt swiftPatternFrSt(qgramIndexFrSt);
    std::cout << "FrStPattern: " << value(host(needle(swiftPatternFrSt)),0) << std::endl;
    // Construct index
    std::cout << "Constructing index..." << std::endl;
    indexRequire(qgramIndexFrSt, QGramSADir());
    std::cout << std::endl;
    */

    // Using Stellars structure for queries
    typedef Pattern<TQGramIndex, Swift<SwiftLocal> > TPattern;
    TPattern swiftPattern(qgramIndex);

    // Construct index
    std::cout << "Constructing index..." << std::endl;
    indexRequire(qgramIndex, QGramSADir());
    std::cout << std::endl;


    // Call Stellar for each database sequence
    double start = sysTime();
    for (unsigned i = 0; i < length(databases); ++i)
    {
        // Using long stellar() to calculate stellarMatches on + strand
        if (stellarOptions.forward)
        {
            TFinder swiftFinder(databases[i], stellarOptions.minRepeatLength, stellarOptions.maxRepeatPeriod);
            stellar(swiftFinder, swiftPattern, stellarOptions.epsilon, stellarOptions.minLength, stellarOptions.xDrop,
                    stellarOptions.disableThresh, stellarOptions.compactThresh, stellarOptions.numMatches,
                    stellarOptions.verbose, databaseIDs[i], true, stellarMatches, AllLocal());
        }

        // Store reverse matches in temporary container and transform before appending to stellarMatches
        TMatches stRevMatches;
        resize(stRevMatches, length(queries));

        // - strand
        if (stellarOptions.reverse)
        {
            reverseComplement(databases[i]);
            TFinder revSwiftFinder(databases[i], stellarOptions.minRepeatLength, stellarOptions.maxRepeatPeriod);
            stellar(revSwiftFinder,
                    swiftPattern,
                    stellarOptions.epsilon,
                    stellarOptions.minLength,
                    stellarOptions.xDrop,
                    stellarOptions.disableThresh,
                    stellarOptions.compactThresh,
                    stellarOptions.numMatches,
                    stellarOptions.verbose,
                    databaseIDs[i],
                    false,
                    stRevMatches, // stellarMatches, // Replace with stRevMatches
                    AllLocal());
            // Convert stRevMatches and append to stellarMatches
            _convertStellarReverseMatches(databases[i], databaseIDs[i], stellarMatches, stRevMatches);
            reverseComplement(databases[i]);
        }
    }
    std::cout << "TIME stellar " << (sysTime() - start) << "s" << std::endl;
}

// ----------------------------------------------------------------------------
// Function _getScore()
// ----------------------------------------------------------------------------

// Get edit distance score of two alignment rows
template <typename TSequence, typename TId, typename TValue>
inline void _getScore(StellarMatch<TSequence, TId> & match, TValue & alignDistance)
{
    SEQAN_ASSERT_EQ(length(match.row1), length(match.row2));
    typedef typename StellarMatch<TSequence, TId>::TRow TRow;
    typedef typename Iterator<TRow>::Type TIter;

    if (!match.orientation)
        reverseComplement(infix(source(match.row1), match.begin1, match.end1));

    TIter itRow1 = begin(match.row1);
    TIter itRow2 = begin(match.row2);
    TIter itRow1End = end(match.row1);

    alignDistance = 0;
    for (; itRow1 != itRow1End; ++itRow1, ++itRow2)
    {
        if (isGap(itRow1) || isGap(itRow2))
            ++alignDistance;
        else if (*itRow1 != *itRow2)
            ++alignDistance;
    }
    // std::cerr << match << std::endl;
    // std::cerr << match.row1 << '\n' << match.row2 << std::endl;

    if (!match.orientation)
        reverseComplement(infix(source(match.row1), match.begin1, match.end1));
}

// ----------------------------------------------------------------------------
// Function _getMatchDistance()
// ----------------------------------------------------------------------------

// Computes distance score for each Stellar match and stores i in distanceScores
// Note: Matches are being sorted within this function
template <typename TScoreAlloc, typename TSequence, typename TId>
void _getMatchDistanceScore(
    StringSet<QueryMatches<StellarMatch<TSequence, TId> > > & queryMatchesSet,
    String<TScoreAlloc> & distanceScores,  // Note that this is a string of strings corresponding to the queryMatchesSet
    unsigned & numThreads)
{
    typedef StellarMatch<TSequence, TId> TMatch;
    typedef typename Size<typename TMatch::TAlign>::Type TSize;
    typedef typename Iterator<String<TMatch>, Standard>::Type TIterator;

    typedef StringSet<QueryMatches<TMatch> > TQueryMatchSet;
    typedef typename Iterator<TQueryMatchSet, Standard>::Type TQueryMatchSetIterator;

    // TODO(rmaerker): Add threads to options and then remove line below or adapt to parallelization scheme.
    omp_set_num_threads(numThreads);
    Splitter<TQueryMatchSetIterator> setSplitter(begin(queryMatchesSet, Standard()), end(queryMatchesSet, Standard()));

    SEQAN_OMP_PRAGMA(parallel for)
    for(int jobId = 0; jobId < static_cast<int>(length(setSplitter)); ++jobId)
    {
        for (TQueryMatchSetIterator it = setSplitter[jobId]; it != setSplitter[jobId + 1]; ++it)
        {
            TScoreAlloc & matchDistanceScores = distanceScores[it - begin(queryMatchesSet, Standard())];
            resize(matchDistanceScores, length((*it).matches));
            // Sorting matches according to query begin position (begin2)
            std::sort(begin((*it).matches), end((*it).matches), CompareStellarMatches<TSequence, TId>());
            TIterator itStellarMatches = begin((*it).matches, Standard());
            TIterator itEndStellarMatches = end((*it).matches, Standard());
            TSize matchIndex = 0;
            for (; itStellarMatches < itEndStellarMatches; goNext(itStellarMatches))
            {
                int alignDistance = 0;  // Calculating distance score of each match using Stellars _analayzeAlignment function
                // Compute edit distance score
                _getScore(*itStellarMatches, alignDistance);
                matchDistanceScores[matchIndex] = alignDistance;
                
                ++matchIndex;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function trimMatchBegin()
// ----------------------------------------------------------------------------

// TODO(ktrappe): Rewrite according to new align modul (disabled atm)
template <typename TMatch, typename TPos>
void _trimMatchBegin(TMatch & stMatch, TPos const & splitPos, TPos const & projSplitPos)
{
    setBeginPosition(stMatch.row2, projSplitPos);
    if (stMatch.orientation)
        setBeginPosition(stMatch.row1, splitPos);
    else
    {
        TPos diff = splitPos - beginPosition(stMatch.row1);
        _transformCoordinates(stMatch);
        setBeginPosition(stMatch.row1, stMatch.begin1 + diff);
        stMatch.begin1 = beginPosition(stMatch.row1);
        _transformCoordinates(stMatch);
    }
    stMatch.begin1 = beginPosition(stMatch.row1);
    stMatch.begin2 = beginPosition(stMatch.row2);
}

// ----------------------------------------------------------------------------
// Function trimMatchEnd()
// ----------------------------------------------------------------------------

// TODO(ktrappe): Rewrite according to new align modul (disabled atm)
template <typename TMatch, typename TPos>
void _trimMatchEnd(TMatch & stMatch, TPos & splitPos, TPos & projSplitPos)
{
    std::cerr << stMatch << '\n';
    std::cerr << stMatch.row1 << '\n' << stMatch.row2 << '\n' << splitPos << '\n' << projSplitPos << '\n';

    setEndPosition(stMatch.row2, projSplitPos);
    if (stMatch.orientation)
    {
        setEndPosition(stMatch.row1, splitPos);

        /*
        SEQAN_ASSERT_GT(length(stMatch.row1), 0u);
        // Clip away trailing gaps if any.
        TPos srcPosEnd = toSourcePosition(stMatch.row1, length(stMatch.row1));
        TPos srcPosLast = toSourcePosition(stMatch.row1, length(stMatch.row1) - 1);
        if (srcPosEnd == srcPosLast)  // trailing gap
        {
            TPos viewPosLast = toViewPosition(stMatch.row1, srcPosLast - 1);
            TPos newClippedEndPos = clippedBeginPosition(stMatch.row1) + viewPosLast + 1;
            setClippedEndPosition(stMatch.row1, newClippedEndPos);
        }

        // Hot fix when adapting end position.
        // TODO(ktrappe): Find out where the problem really occurs.
        setClippedEndPosition(stMatch.row2, clippedBeginPosition(stMatch.row2) + length(stMatch.row1));
        */
    }
    else
    {
        TPos diff = endPosition(stMatch.row1) - splitPos;
        _transformCoordinates(stMatch);
        std::cerr << stMatch << '\n';
        std::cerr << stMatch.row1 << '\n' << stMatch.row2 << '\n' << splitPos << '\n' << projSplitPos << '\n';
        setEndPosition(stMatch.row1, stMatch.end1 - diff);
        std::cerr << stMatch << '\n';
        std::cerr << stMatch.row1 << '\n' << stMatch.row2 << '\n' << splitPos << '\n' << projSplitPos << '\n';
        stMatch.end1 = endPosition(stMatch.row1);
        _transformCoordinates(stMatch);
    }
    std::cerr << stMatch << '\n';
    std::cerr << stMatch.row1 << '\n' << stMatch.row2 << '\n' << splitPos << '\n' << projSplitPos << '\n';

    stMatch.end1 = endPosition(stMatch.row1);
    stMatch.end2 = endPosition(stMatch.row2);
}

// ----------------------------------------------------------------------------
// Function getStellarIndel()
// ----------------------------------------------------------------------------

// TODO(ktrappe): Restructure and rewrite (buggy) (disabled atm)
template <typename TSequence, typename TId, typename TBreakpoint>
void _getStellarIndel(StellarMatch<TSequence, TId> & match,
                      String<TBreakpoint> & globalStellarIndels,
                      TId const & queryId,
                      TSequence & query)
{
    // std::cerr << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    typedef typename Infix<TSequence>::Type TInfix;
    typedef typename TBreakpoint::TPos TPos;
    // typedef typename StellarMatch<TSequence, TId>::TRow TRow;
    typedef Align<TInfix> TAlign;
    typedef typename Row<TAlign>::Type TRow;

    bool gapOpen = true;
    TPos startSeqPos, endSeqPos, indelStart, pos;
    indelStart = 0;

    // Match refinement
    int matchNumber, alignLen, score;
    Score<int> scoreType(1, -2, -2, -5);

    /*
    std::cerr << "Trying to build infices: " << length(match.row1) << " begin1: " << match.begin1 << " clippedBegin: " << clippedBeginPosition(match.row1) <<
        " end1: " << match.end1 << " clippedEnd: " << clippedEndPosition(match.row1) << std::endl;
    */
    TInfix seq1, seq2;
    seq1 = infix(source(match.row1), match.begin1, match.end1);

    seq2 = infix(query, match.begin2, match.end2);

    if (!match.orientation)
        reverseComplement(seq1);

    TAlign align;
    resize(rows(align), 2);
    setSource(row(align, 0), seq1);
    setSource(row(align, 1), seq2);

    StringSet<TInfix> seqs;

    appendValue(seqs, seq1);
    appendValue(seqs, seq2);

    if (length(match.row1) != length(match.row2))
    {
        std::cerr << match.row1 << '\t' << match.row2 << '\n';
        std::cerr << match << '\n';
    }
    _analyzeAlignment(match.row1, match.row2, alignLen, matchNumber);
    score = alignLen - matchNumber;

    // Refine match using banded global alignment with affine gap score
    if (score > 0)
    {
        globalAlignment(align, scoreType, -score, score, Gotoh());
        // globalAlignment(align, seqs, scoreType, -score, score, Gotoh());
        // std::cout << align << std::endl;
    }

    TRow & rowDB = row(align, 0);
    TRow & rowRead = row(align, 1);

    // Insertions
    for (pos = 0; pos < (TPos) std::max(endPosition(rowDB), endPosition(rowRead)); ++pos)
    {
        // Look for gap in rowDB
        if (isGap(rowDB, pos))
        {
            // Keep track if it is the first gap, i.e. the beginning of the insertion
            if (gapOpen)
            {
                gapOpen = false;
                indelStart = pos;
                // std::cerr << "indelStart: " << indelStart << " gapopen " << gapOpen << std::endl;
            }
        }
        else
        {
            // If an insertion just closed, i.e. !gapOpen, save the insertion in a breakpoint
            if (!gapOpen)
            {
                // Get breakpoint position (note, in case of an insertion start=end)
                startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
                TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
                TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
                TBreakpoint bp(match.id,
                               match.id,
                               match.orientation,
                               match.orientation,
                               startSeqPos,
                               startSeqPos,
                               readStartPos,
                               readEndPos,
                               queryId);
                // get insertion sequence
                TPos bPos = toSourcePosition(rowRead, indelStart);
                TPos ePos = toSourcePosition(rowRead, pos);
                if (bPos > ePos)
                    std::swap(bPos, ePos);
                TInfix inSeq = infix(query, bPos, ePos);
                setSVType(bp, static_cast<TId>("insertion"));
                setInsertionSeq(bp, inSeq);
                // _insertBreakpoint(globalStellarIndels, bp);
                // std::cerr << bp << std::endl;
                // reset gapOpen
                gapOpen = true;
            }
        }
    }
    // Get insertion at the end of the match
    if (!gapOpen)
    {
        // Get breakpoint position (note, in case of an insertion start=end)
        startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
        TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
        TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
        TBreakpoint bp(match.id,
                       match.id,
                       match.orientation,
                       match.orientation,
                       startSeqPos,
                       startSeqPos,
                       readStartPos,
                       readEndPos,
                       queryId);
        // get insertion sequence
        TPos bPos = toSourcePosition(rowRead, indelStart) + match.begin2;
        TPos ePos = toSourcePosition(rowRead, pos) + match.begin2;
        if (bPos > ePos)
            std::swap(bPos, ePos);
        TInfix inSeq = infix(query, bPos, ePos);
        setSVType(bp, static_cast<TId>("insertion"));
        setInsertionSeq(bp, inSeq);
        _insertBreakpoint(globalStellarIndels, bp);
        // reset gapOpen
        gapOpen = true;
    }

    // Deletions
    for (pos = 0; pos < (TPos) std::max(endPosition(rowDB), endPosition(rowRead)); ++pos)
    {
        // Look for gap in rowRead
        if (isGap(rowRead, pos))
        {
            if (gapOpen)
            {
                gapOpen = false;
                indelStart = pos;
            }
        }
        else
        {
            if (!gapOpen)
            {
                // Get breakpoint positions
                startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
                endSeqPos = toSourcePosition(rowDB, pos) + match.begin1;
                if (startSeqPos > endSeqPos)
                    std::swap(startSeqPos, endSeqPos);
                TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
                TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
                TBreakpoint bp(match.id,
                               match.id,
                               match.orientation,
                               match.orientation,
                               startSeqPos,
                               endSeqPos,
                               readStartPos,
                               readEndPos,
                               queryId);
                setSVType(bp, static_cast<TId>("deletion"));
                _insertBreakpoint(globalStellarIndels, bp);
                gapOpen = true;
            }
        }
    }
    // Get deletion at the end of the match
    if (!gapOpen)
    {
        startSeqPos = toSourcePosition(rowDB, indelStart) + match.begin1;
        endSeqPos = toSourcePosition(rowDB, pos) + match.begin1;
        if (startSeqPos > endSeqPos)
            std::swap(startSeqPos, endSeqPos);
        TPos readStartPos = toSourcePosition(rowRead, indelStart) + match.begin2;
        TPos readEndPos = toSourcePosition(rowRead, pos) + match.begin2;
        TBreakpoint bp(match.id,
                       match.id,
                       match.orientation,
                       match.orientation,
                       startSeqPos,
                       endSeqPos,
                       readStartPos,
                       readEndPos,
                       queryId);
        setSVType(bp, static_cast<TId>("deletion"));
        _insertBreakpoint(globalStellarIndels, bp);
        gapOpen = true;
    }

    if (!match.orientation)
        reverseComplement(seq1);
}

// /////////////////////////////////////////////////////////////////////////////
// Functions taken from Stellar code for writing Stellar parameters and read files

template <typename TId>
struct IdComparator
{
    bool operator()(TId const & id1, TId const & id2)
    {
        return std::lexicographical_compare(begin(id1, Standard()), end(id1, Standard()), begin(id2, Standard()), end(id2, Standard()));
    }
};

// /////////////////////////////////////////////////////////////////////////////
// Imports sequences from a file,
//  stores them in the StringSet seqs and their identifiers in the StringSet ids
template <typename TSequence, typename TId>
inline bool
_importSequences(CharString const & fileName,
                 CharString const & name,
                 StringSet<TSequence> & seqs,
                 StringSet<TId> & ids)
{
    typedef typename Iterator<StringSet<TId>, Standard>::Type TIdSetIterator;

    seqan::SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(fileName)))
    {
        std::cerr << "Failed to open " << name << " file.\n";
        return false;
    }

    String<Iupac> seq;
    TId id;
    TId sId;
    unsigned seqCount = 0;
    for (; !atEnd(seqFileIn); ++seqCount)
    {
        readRecord(id, seq, seqFileIn);
        appendValue(seqs, seq, Generous());

        _getShortId(sId, id);
        appendValue(ids, sId, Generous());
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
    std::cout << "Loaded " << seqCount << " " << name << " sequence" << ((seqCount > 1) ? "s." : ".") << std::endl;
    return true;
}

// /////////////////////////////////////////////////////////////////////////////
// Calculates parameters from parameters in options object and writes them to std::cout
void _writeCalculatedParams(StellarOptions & options)
{
// IOREV _notio_
    int errMinLen = (int) floor(options.epsilon * options.minLength);
    int n = (int) ceil((errMinLen + 1) / options.epsilon);
    int errN = (int) floor(options.epsilon * n);
    unsigned smin = (unsigned) _min(ceil((double)(options.minLength - errMinLen) / (errMinLen + 1)),
                                    ceil((double)(n - errN) / (errN + 1)));

    std::cout << "Calculated parameters:" << std::endl;
    if (options.qGram == (unsigned)-1)
    {
        options.qGram = (unsigned)_min(smin, 32u);
        std::cout << "  k-mer length: " << options.qGram << std::endl;
    }

    int threshold = (int) _max(1, (int) _min((n + 1) - options.qGram * (errN + 1),
                                             (options.minLength + 1) - options.qGram * (errMinLen + 1)));
    int overlap = (int) floor((2 * threshold + options.qGram - 3) / (1 / options.epsilon - options.qGram));
    int distanceCut = (threshold - 1) + options.qGram * overlap + options.qGram;
    int logDelta = _max(4, (int) ceil(log((double)overlap + 1) / log(2.0)));
    int delta = 1 << logDelta;

    std::cout << "  s^min       : " << smin << std::endl;
    std::cout << "  threshold   : " << threshold << std::endl;
    std::cout << "  distance cut: " << distanceCut << std::endl;
    std::cout << "  delta       : " << delta << std::endl;
    std::cout << "  overlap     : " << overlap << std::endl;
    std::cout << std::endl;
}

// /////////////////////////////////////////////////////////////////////////////
// Writes user specified parameters from options object to std::cout
void
_writeSpecifiedParams(StellarOptions & options)
{
// IOREV _notio_
    // Output user specified parameters
    std::cout << "User specified parameters:" << std::endl;
    std::cout << "  minimal match length             : " << options.minLength << std::endl;
    std::cout << "  maximal error rate (epsilon)     : " << options.epsilon << std::endl;
    std::cout << "  maximal x-drop                   : " << options.xDrop << std::endl;
    if (options.qGram != (unsigned)-1)
        std::cout << "  k-mer (q-gram) length            : " << options.qGram << std::endl;
    std::cout << "  search forward strand            : " << ((options.forward) ? "yes" : "no") << std::endl;
    std::cout << "  search reverse complement        : " << ((options.reverse) ? "yes" : "no") << std::endl;
    std::cout << std::endl;

    std::cout << "  verification strategy            : " << options.fastOption << std::endl;
    if (options.disableThresh != (unsigned)-1)
    {
        std::cout << "  disable queries with more than   : " << options.disableThresh << " matches" << std::endl;
    }
    std::cout << "  maximal number of matches        : " << options.numMatches << std::endl;
    std::cout << "  duplicate removal every          : " << options.compactThresh << std::endl;
    if (options.maxRepeatPeriod != 1 || options.minRepeatLength != 1000)
    {
        std::cout << "  max low complexity repeat period : " << options.maxRepeatPeriod << std::endl;
        std::cout << "  min low complexity repeat length : " << options.minRepeatLength << std::endl;
    }
    if (options.qgramAbundanceCut != 1)
    {
        std::cout << "  q-gram abundance cut ratio       : " << options.qgramAbundanceCut << std::endl;
    }
    std::cout << std::endl;
}

// /////////////////////////////////////////////////////////////////////////////
// Writes file name from options object to std::cout
void
_writeFileNames(StellarOptions & options)
{
// IOREV _notio_
    std::cout << "Database file   : " << options.databaseFile << std::endl;
    std::cout << "Query file      : " << options.queryFile << std::endl;
    std::cout << std::endl;
}

#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_MSPLAZER_STELLAR_ROUTINES_H_
