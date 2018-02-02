// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// RABEMA is a read mapping benchmark tool.  See README for an overview.
//
// This is the program for building the RABEMA gold standard.  The input to
// this program is the reference sequence as a FASTA file and a "perfect
// mapping" SAM/BAM file (see README on how to obtain such a gold standard).
// ==========================================================================

// TODO(holtgrew): Another possible way to save memory is not to store the first, always identical part of the read id.

#include <fstream>
#include <iostream>
#include <map>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/seq_io.h>

#include "curve_smoothing.h"
#include "find_hamming_simple_ext.h"
#include "find_myers_ukkonen_reads.h"
#include "io_gsi.h"
#include "sorting.h"

// ============================================================================
// Enums, Tags, Classes, Typedefs.
// ============================================================================

typedef std::map<int, TWeightedMatches> TErrorCurves;

// ----------------------------------------------------------------------------
// Helper Class IntervalizeCmp
// ----------------------------------------------------------------------------

// Comparison functor used in intervalizeErrorCurves().

struct IntervalizeCmp
{
    StringSet<CharString> const & readNameStore;

    IntervalizeCmp(StringSet<CharString> const & readNameStore) :
        readNameStore(readNameStore)
    {}

    bool operator()(unsigned lhs, unsigned rhs) const
    {
        return lessThanSamtoolsQueryName(readNameStore[lhs], readNameStore[rhs]);
    }

};

// ----------------------------------------------------------------------------
// Class ContigInterval
// ----------------------------------------------------------------------------

// Represents an interval in a contig.
struct ContigInterval
{
    // Contig the interval is in.
    size_t contigId;

    // true iff the interval is on the forward strand.
    bool isForward;

    // First value of interval
    size_t first;

    // Last value of interval.
    size_t last;

    // Default constructor so it can be used in containers.
    ContigInterval() : contigId(0), isForward(false), first(0), last(0)
    {}

    // Constructor for the record.
    ContigInterval(size_t _contigId, bool _isForward, size_t _first, size_t _last) :
        contigId(_contigId), isForward(_isForward), first(_first), last(_last)
    {}
};

// ---------------------------------------------------------------------------
// Enum DistanceMetric
// ---------------------------------------------------------------------------

enum DistanceMetric
{
    HAMMING_DISTANCE,
    EDIT_DISTANCE
};

// ---------------------------------------------------------------------------
// Class BuildGoldStandardOptions
// ---------------------------------------------------------------------------

struct BuildGoldStandardOptions
{
    // Verbosity level, quiet: 0, normal: 1, verbose: 2, very verbose: 3.
    int verbosity;

    // Whether N should match all characters without penalty.
    bool matchN;

    // Whether or not to run in oracle mode.  The input SAM/BAM file must only contain the sample location for each read
    // in this case.  This is used when building a GSI file for reads simulated by Mason where the original locations
    // are written out in a SAM/BAM file.
    bool oracleMode;

    // Maximal error rate in percent, relative to the read length.
    int maxError;

    // Indicates whether a maximal error was set.
    bool maxErrorSet;

    // Distance metric to use.
    DistanceMetric distanceMetric;

    // Path to the output file that we will write the GSI to.
    seqan::CharString outGsiPath;

    // Path to the reference contigs file.
    seqan::CharString referencePath;

    // Exactly one of the following two has to be given.
    //
    // Path to the perfect input SAM/BAM file.
    seqan::CharString inBamPath;

    BuildGoldStandardOptions() :
        verbosity(1),
        matchN(false),
        oracleMode(false),
        maxError(0),
        maxErrorSet(false),
        distanceMetric(EDIT_DISTANCE)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Function metricName()
// ---------------------------------------------------------------------------

char const * metricName(DistanceMetric met)
{
    if (met == EDIT_DISTANCE)
        return "edit";
    else  // if (met == HAMMING_DISTANCE)
        return "hamming";
}

// ---------------------------------------------------------------------------
// Function yesNo()
// ---------------------------------------------------------------------------

char const * yesNo(bool b)
{
    if (b)
        return "yes";
    else
        return "no";
}

// ----------------------------------------------------------------------------
// Function trimSeqHeaderToId()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Such a method is already in the SeqAn library, but in extras. Remove here when it is in core.
void trimSeqHeaderToId(seqan::CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
            break;
    resize(header, i);
}

// ----------------------------------------------------------------------------
// Function intervalizeAndDumpErrorCurves()
// ----------------------------------------------------------------------------

// Build intervals from the error curves and print them to stream.

template <typename TStream>
int intervalizeAndDumpErrorCurves(TStream & stream,
                                  TErrorCurves const & errorCurves,
                                  String<int> const & readAlignmentDistances,
                                  StringSet<CharString> const & readNameStore,
                                  StringSet<CharString> const & contigNameStore,
                                  BuildGoldStandardOptions const & options,
                                  bool toStdout)
{
    if (!toStdout)
        std::cerr << "\n____POINT TO INTERVAL CONVERSION______________________________________________\n\n"
                  << "Progress: ";

    // Write out the header.
    GsiHeader header;
    writeHeader(stream, header, Gsi());
    writeRecord(stream, GSI_COLUMN_NAMES, Gsi());

    // Get a list of read ids, sorted by their read name.
    String<unsigned> sortedReadIds;
    reserve(sortedReadIds, length(readNameStore));
    for (unsigned i = 0; i < length(readNameStore); ++i)
        appendValue(sortedReadIds, i);
    IntervalizeCmp cmp(readNameStore);
    std::sort(begin(sortedReadIds, Standard()), end(sortedReadIds, Standard()), cmp);

    unsigned tenPercent = errorCurves.size() / 10 + 1;
    typedef TErrorCurves::const_iterator TErrorCurvesIter;
    for (unsigned i = 0; i < length(sortedReadIds); ++i)
    {
        TErrorCurvesIter it = errorCurves.find(sortedReadIds[i]);
        if (it == errorCurves.end())
        {
            if (!toStdout)
                std::cerr << "WARNING: Something went wrong with read ids! This should not happen.\n";
            continue;
        }

        if (!toStdout)
        {
            if (tenPercent > 0u && i % tenPercent == 0u)
                std::cerr << i / tenPercent * 10 << '%';
            else if (tenPercent > 5u && i % (tenPercent / 5) == 0u)
                std::cerr << '.';
        }

        size_t readId = it->first;
        TWeightedMatches const & matches = it->second;

        // Sort the matches.  Matches with high scores (negative score and low absolute value) come first.
        TWeightedMatches sortedMatches(matches);

        // intervals[e] holds the intervals for error e of the current read.
        String<String<ContigInterval> > intervals;
        int maxError = options.oracleMode ? 0 : (int)options.maxError;
        if (options.oracleMode)
            for (unsigned j = 0; j < length(matches); ++j)
                maxError = std::max(maxError, matches[j].distance);
        resize(intervals, maxError + 1);

        // Join the intervals stored in sortedMatches.
        //
        // The position of the previous match, so we can consider only the ones with the smallest error.
        //
        // The following two vars should be != first pos and contigId.
        size_t previousPos = std::numeric_limits<size_t>::max();
        size_t previousContigId = std::numeric_limits<size_t>::max();
        typedef Iterator<TWeightedMatches>::Type TWeightedMatchesIter;
        for (TWeightedMatchesIter it = begin(sortedMatches);
             it != end(sortedMatches); ++it)
        {
            // Skip it if (it - 1) pointed to same pos (and must point to
            // one with smaller absolute distance.
            if (it->pos == previousPos && it->contigId == previousContigId)
                continue;
            // Consider all currently open intervals with a greater error than the error in *it and extend them or
            // create a new one.
            int error = options.oracleMode ? 0 : abs(it->distance);
            SEQAN_ASSERT_LEQ(error, maxError);
            for (int e = error; e <= maxError; ++e)
            {
                // Handle base case of no open interval:  Create new one.
                if (length(intervals[e]) == 0)
                {
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
                    continue;
                }
                ContigInterval & interval = back(intervals[e]);
                // Either extend the interval or create a new one.
                if (interval.contigId == it->contigId && interval.isForward == it->isForward)
                    SEQAN_ASSERT_LEQ(interval.last, it->pos);
                if (interval.contigId == it->contigId && interval.isForward == it->isForward &&
                    interval.last + 1 == it->pos)
                    back(intervals[e]).last += 1;
                else
                    appendValue(intervals[e], ContigInterval(it->contigId, it->isForward, it->pos, it->pos));
            }
            // Book-keeping.
            previousPos = it->pos;
            previousContigId = it->contigId;
        }

        // Print the resulting intervals.
        typedef Iterator<String<String<ContigInterval> > >::Type TIntervalContainerIter;
        int distance = 0;
        for (TIntervalContainerIter it = begin(intervals);
             it != end(intervals); ++it, ++distance)
        {
            typedef Iterator<String<ContigInterval> >::Type TIntervalIter;
            for (TIntervalIter it2 = begin(*it); it2 != end(*it); ++it2)
            {
                int flags = 0;
                // We appended custom prefixes to the read ids.  "/S" means single-end, "/0" means left mate, "/1" means
                // right mate.
                int mateNo = -1;
                if (back(readNameStore[readId]) == '0')
                    mateNo = 0;
                if (back(readNameStore[readId]) == '1')
                    mateNo = 1;
                SEQAN_ASSERT_EQ(readNameStore[readId][length(readNameStore[readId]) - 2], '/');
                if (mateNo == 0)
                    flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_FIRST_MATE;
                else if (mateNo == 1)
                    flags = GsiRecord::FLAG_PAIRED | GsiRecord::FLAG_SECOND_MATE;

                int gsiDistance = options.oracleMode ? readAlignmentDistances[readId] : distance;
                if (options.oracleMode && options.maxErrorSet && gsiDistance > options.maxError)
                    continue;  // Skip if cut off in Rabema oracle mode.
                GsiRecord gsiRecord(CharString(prefix(readNameStore[readId], length(readNameStore[readId]) - 2)), flags,
                                    gsiDistance, contigNameStore[it2->contigId],
                                    it2->isForward, it2->first, it2->last);
                writeRecord(stream, gsiRecord, Gsi());
            }
        }
    }
    if (!toStdout)
        std::cerr << "100% DONE\n";

    return 0;
}

// ----------------------------------------------------------------------------
// Function ceilAwayFromZero()
// ----------------------------------------------------------------------------

// Ceil away from Zero.
//
// ceilAwayFromZero(-1.5) == -2
// ceilAwayFromZero(1.5) == 2
template <typename T>
inline T ceilAwayFromZero(T x)
{
    if (x < 0)
        return floor(x);

    return ceil(x);
}

// ----------------------------------------------------------------------------
// Function buildErrorCurvePoint()
// ----------------------------------------------------------------------------

// TODO(holtgrew): This function has still a lot of cleanup potential.

// Build the error curve points around the end position of the given contig.
//
// This is pretty involved and this function easily is the most complex one.
//
// Returns rightmost border of the added points.
template <typename TContigSeq, typename TReadSeq, typename TPatternSpec, typename TReadNames>
void buildErrorCurvePoints(String<WeightedMatch> & errorCurve,
                           int & maxError,
                           TContigSeq /*const*/ & contig,
                           size_t contigId,
                           bool isForward,
                           TReadSeq /*const*/ & read,
                           size_t readId,
                           size_t endPos,
                           TReadNames const & readNames,
                           bool matchN,
                           TPatternSpec const &)
{
    typedef typename Position<TContigSeq>::Type TPosition;

    // In oracle Sam mode, the maximum error is the error at the position given in the Sam alignment.
    bool oracleMode = false;
    if (maxError == std::numeric_limits<int>::max())  // We are in oracle mode.
    {
        oracleMode = true;
        Finder<TContigSeq> finder(contig);
        Pattern<TReadSeq, TPatternSpec> pattern(read, -(int)length(read) * 40);
        bool ret = setEndPosition(finder, pattern, endPos);
        (void) ret; // If compiled without assertions.
        SEQAN_ASSERT(ret);
        maxError = -getScore(pattern);
    }

    // Debug-adjustments.
    #define ENABLE 0
    #define ALL 1
    #define READID 1

    //if (readNames[readId] == CharString("rabema_synthetic_ecoli_illumina_l100_100k_e3.000000376")) {
    //      std::cerr << "**************** read id = " << readId << " readname = " << readNames[readId] << " read = " << read << " maxError == " << maxError << std::endl;
    //}

    // The read maps with less than the given number of errors to [left, right]
    // to the contig sequence, is forward strand iff is Forwar is true.
    TPosition /*left = endPos,*/ right = endPos;

    // Skip this alignment position if not right of previously
    // right border of the interval.
    // TODO(holtgrew): Remove this piece, it's unused right now. We could/should also fix this to work for streaming over SAM files.
    //if (readId == previousReadId && contigId == previousContigId && left <= previousRightBorder) {
    //    //std::cerr << __FILE__  << ":" << __LINE__ << " Skipping because " << left << " <= " << previousRightBorder << std::endl;
    //    return previousRightBorder;
    //}

    bool ret;  // Flag used for assertions below.
    (void) ret;  // If run without assertions.
    int relativeMinScore = (int)ceilAwayFromZero(100.0 * -maxError / length(read));

    // Setup the finder and pattern.
    Finder<TContigSeq> finder(contig);
    Pattern<TReadSeq, TPatternSpec> pattern(read, -(int)length(read) * 40);
    // If configured so, match N against all other values, otherwise match
    // against none.
    _patternMatchNOfPattern(pattern, matchN);
    _patternMatchNOfFinder(pattern, matchN);
    TPosition hitBeginPosition;

    if (ENABLE && (ALL || readId == READID))
    {
        std::cerr << "**************** read id = " << readId << " readname = " << readNames[readId]
                  << " read = " << read << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " maxError = " << maxError << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " relative min score = " << relativeMinScore << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " extending to the right." << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " contig id == " << contigId << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " is forward strand? " << isForward << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " contig length = " << length(contig) << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " endPos = " << endPos << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << " infix(contig, endPos - length(read), endPos) == "
                  << infix(contig, endPos - length(read), endPos) << std::endl;
    }

    // We will first gather all results in tempMatches.  Below, we
    // will smooth the curve in these points and throw out too bad
    // entries.  Then, we will append tempMatches to errorCurve.
    String<WeightedMatch> tempMatches;

    // First, extend the interval to the right.
    {
        // Skip to original hit.
        ret = setEndPosition(finder, pattern, endPos);
        SEQAN_ASSERT(ret);
        SEQAN_ASSERT_EQ(endPos, endPosition(finder));
        SEQAN_ASSERT_GEQ(getScore(pattern), -maxError);
        while (findBegin(finder, pattern, getScore(pattern)))
            continue;  // Find leftmost begin position.
        SEQAN_ASSERT_GT(endPos, beginPosition(finder));
        SEQAN_ASSERT_EQ(getScore(pattern), getBeginScore(pattern));

        // Add original hit to the error curve points.
        int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
        appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore,
                                               beginPosition(finder)));
        hitBeginPosition = beginPosition(finder);
        if (ENABLE && (ALL || readId == READID))
        {
            std::cerr << __FILE__ << ":" << __LINE__ << " -- getScore(pattern) == " << getScore(pattern) << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches)
                      << " for read id " << readId << " (FIRST HIT)" << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read " << read
                      << " endPos = " << endPos << std::endl;
        }

        // Now extend to the first hit with too low score.
        bool foundWithTooLowScore = false;
        while (find(finder, pattern))
        {
            while (findBegin(finder, pattern, getScore(pattern)))
                continue;  // Find leftmost begin position.
            if (getScore(pattern) < -maxError && beginPosition(finder) != back(tempMatches).beginPos)
            {
                foundWithTooLowScore = true;
                if (ENABLE && (ALL || readId == READID))
                {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- Found too low score." << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- low scoring match was "
                              << WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore,
                                     beginPosition(finder)) << std::endl;
                }
                break;
            }
            int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
            appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore,
                                                   beginPosition(finder)));
            if (ENABLE && (ALL || readId == READID))
            {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches)
                          << " for read id " << readId << std::endl;
                std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder)
                          << " read " << read << std::endl;
            }
            right += 1;
        }
        // If we broke because of the score limit then collect the last not
        // yet added hit and the ones right of it until the beginPosition
        // changes.
        if (foundWithTooLowScore)
        {
            if (beginPosition(finder) == hitBeginPosition)
            {
                relativeScore = static_cast<int>(ceilAwayFromZero(100.0 * getScore(pattern) / length(read)));
                appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1,
                                                       relativeScore, beginPosition(finder)));
                TPosition currentBeginPosition = beginPosition(finder);
                // Add the rest until we hit one with a different begin position.
                //
                // Loop at most length(read) times (this limit is here
                // for the quality based DP algorithm which can suffer
                // from "infinite inserts").
                for (unsigned i = 0; find(finder, pattern) && i < length(read); ++i)
                {
                    while (findBegin(finder, pattern, getScore(pattern)))
                        continue;  // Find leftmost begin position.
                    SEQAN_ASSERT(ret);
                    if (beginPosition(finder) != currentBeginPosition)
                        break;
                    relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1,
                                                           relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID))
                    {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches)
                                  << " for read id " << readId << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder)
                                  << " read " << read << std::endl;
                    }
                    right += 1;
                }
            }
        }

        if (ENABLE && (ALL || readId == READID))
        {
            std::cerr << __FILE__ << ":" << __LINE__ << " extending to the left." << std::endl;
        }
        // Then, extend the interval to the left.
        {
            // Length by which we extend the interval.  Note that this must be
            // at least the length of the read + max error count so no special
            // treatment of islands on the left side must be done.
            TPosition kIntervalLen = length(read) + maxError;
            // Tentatively extend the interval to the left.  We will
            // extend until we hit a position with too low score.
            TPosition tentativeLeft = endPos;
            // Flag that indicates whether we found an entry with too low score.
            bool foundTooLowScore = false;
            // Flag that indicates whether we performed a loop iteration after having found a too low score.
            bool loopedAfterFoundTooLowScore = false;
            // Flag for breaking out of loop if we hit the right border of the last interval.
            bool hitLastRight = false;
            while (tentativeLeft > 0 && !loopedAfterFoundTooLowScore && !hitLastRight)
            {
                // Stop if went went to the rightmost position of the
                // last interval with the previous loop iteration.
                //if (readId == previousReadId && contigId == previousContigId && tentativeLeft == previousRightBorder)
                //    break;
                loopedAfterFoundTooLowScore = foundTooLowScore;
                TPosition oldTentativeLeft = tentativeLeft;
                // Move the tentative left position left by kIntervalLen but not
                // further left than 0.
                if (ENABLE && (ALL || readId == READID))
                    std::cerr << "tentativeLeft before = " << tentativeLeft << std::endl;
                tentativeLeft -= _min(kIntervalLen, tentativeLeft);
                if (ENABLE && (ALL || readId == READID))
                    std::cerr << "tentativeLeft after = " << tentativeLeft << std::endl;
                // Do not go further than the previous right position.
                //if (readId == previousReadId && contigId == previousContigId && tentativeLeft <= previousRightBorder) {
                //    tentativeLeft = previousRightBorder;
                //    hitLastRight = true;
                //}
                if (ENABLE && (ALL || readId == READID))
                {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- tentative left = " << tentativeLeft << std::endl;
                }
                // Search from tentative left position to the previous tentative left position.
                ret = setEndPosition(finder, pattern, tentativeLeft);
                SEQAN_ASSERT(ret);
                if (ENABLE && (ALL || readId == READID))
                {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- endPosition(finder) = "
                              << endPosition(finder) << std::endl;
                }
                if (endPosition(finder) > oldTentativeLeft)
                    break;  // Could not set position of the finder left of old tentative left.
                while (findBegin(finder, pattern, getScore(pattern)))
                    continue;  // Find leftmost begin position.
                int relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1, relativeScore,
                                                       beginPosition(finder)));
                if (ENABLE && (ALL || readId == READID))
                {
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches)
                              << " for read id " << readId << std::endl;
                    std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read "
                              << read << std::endl;
                }
                foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                while (find(finder, pattern) && endPosition(finder) != oldTentativeLeft)
                {
                    while (findBegin(finder, pattern, getScore(pattern)))
                        continue;  // Find leftmost begin position.
                    SEQAN_ASSERT(ret);
                    relativeScore = (int)ceilAwayFromZero(100.0 * getScore(pattern) / length(read));
                    appendValue(tempMatches, WeightedMatch(contigId, isForward, endPosition(finder) - 1,
                                                           relativeScore, beginPosition(finder)));
                    if (ENABLE && (ALL || readId == READID))
                    {
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- raw score is " << getScore(pattern)
                                  << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- appended " << back(tempMatches)
                                  << " for read id " << readId << std::endl;
                        std::cerr << __FILE__ << ":" << __LINE__ << " -- infix " << infix(finder) << " read "
                                  << read << std::endl;
                    }
                    foundTooLowScore = foundTooLowScore || (relativeScore < relativeMinScore);
                }
            }
            if (ENABLE && (ALL || readId == READID))
            {
                std::cerr << __FILE__ << ":" << __LINE__ << " -- after loop" << std::endl;
            }
            /*
            // Now we can be sure that the temporary matches contain an entry
            // left of the first one with a too low score that has a different
            // begin position than hitBeginPosition.
            std::sort(begin(tempMatches, Standard()), end(tempMatches, Standard()));
            ModifiedString<String<WeightedMatch>, ModReverse> rTempMatches(tempMatches);
            TPosition prefixLength = 0;
            TPosition i;
            // Search up to the first too low score.
            for (i = 0; i < length(rTempMatches); ++i) {
                if (ENABLE && (ALL || readId == READID)) {
                    std::cerr << "Considering " << rTempMatches[i] << std::endl;
                }
                if (rTempMatches[i].distance < relativeMinScore) {
                    if (rTempMatches[i].beginPos == hitBeginPosition) {
                        prefixLength += 1;
                        left -= 1;
                        i += 1;
                    }
                    break;
                }
                prefixLength += 1;
                left -= 1;
            }
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << "prefixLength == " << prefixLength << std::endl;
                std::cerr << "i == " << i << std::endl;
            }
            // Then, append until we find a position with a different begin position.
            for (; i < length(rTempMatches); ++i) {
                if (rTempMatches[i].beginPos == hitBeginPosition) {
                    prefixLength += 1;
                    left -= 1;
                } else {
                    break;
                }
            }
            // Finally, append the prefix of the given length.
            if (ENABLE && (ALL || readId == READID)) {
                std::cerr << "appending prefix(rTempMatches, " << prefixLength << ")" << std::endl;
            }
        */
        }
    }
    // Postprocessing: Sorting, smoothing, filtering.
    std::sort(begin(tempMatches, Standard()), end(tempMatches, Standard()));
    smoothErrorCurve(tempMatches);
    appendValue(tempMatches, WeightedMatch(0, 0, 0, relativeMinScore - 1, 0));  // Sentinel.
    if (oracleMode)
    {
        // In oracle Sam mode, we only want the lake with last pos endPos-1.
        String<WeightedMatch> buffer;
        bool flag = false;
        for (size_t i = 0; i < length(tempMatches); ++i)
        {
            if (tempMatches[i].distance < relativeMinScore)
            {
                if (flag)
                {
                    append(errorCurve, buffer);
                    break;
                }
                else
                {
                    clear(buffer);
                }
            }
            else
            {
                appendValue(buffer, tempMatches[i]);
                if (tempMatches[i].pos == endPos - 1)
                    flag = true;
            }
        }
    }
    else
    {
        for (size_t i = 0; i < length(tempMatches); ++i)
        {
            if (tempMatches[i].distance >= relativeMinScore)
                appendValue(errorCurve, tempMatches[i]);
        }
    }

//     if (readId == READID) {
//         std::cerr << ",-- errorCurve is (read id = " << readId << " readname = " << readNames[readId] << ")" << std::endl;
//         for (unsigned i = 0; i < length(errorCurve); ++i) {
//             std::cerr << "| " << errorCurve[i] << std::endl;
//         }
//         std::cerr << "`--" << std::endl;
//     }

//     std::cerr << __FILE__ << ":" << __LINE__ << " return " << right << std::endl;
}

// ----------------------------------------------------------------------------
// Function matchesToErrorFunction()
// ----------------------------------------------------------------------------

// Compute error curve from the reads and reference sequences.

template <typename TPatternSpec>
int matchesToErrorFunction(TErrorCurves & errorCurves,
                           String<int> & readAlignmentDistances,  // only used in case of oracle mode
                           BamFileIn & inBam,
                           StringSet<CharString> & readNameStore,
                           StringSet<CharString> const & contigNameStore,
                           FaiIndex const & faiIndex,
                           BuildGoldStandardOptions const & options,
                           TPatternSpec const & /*patternTag*/)
{
    double startTime = 0;

    if (atEnd(inBam))
        return 0;  // Do nothing if there are no more alignments in the file.

    // We store the read names and append "/0" and "/1", depending on their flag value ("/0" for first, "/1" for last
    // flag).  Sadly, there is no easy way out here.  Single-end reads are stored as "/S".
    NameStoreCache<StringSet<CharString> > readNameStoreCache(readNameStore);
    String<unsigned> readLengthStore;

    std::cerr << "\n____EXTENDING INTERVALS_______________________________________________________\n\n"
              << "Each dot represents work for 100k nucleotides in the genome.\n";

    // Stream over SAM/BAM file.  The main assumption here is that the reads are sorted by coordinate, which is check
    // when loading the header.
    //
    // For each read alignments:
    //   If aligned on forward strand:
    //     Build error curve for this read alignment on forward strand.
    //   Else:
    //     Build error curve for this read alignment on backward strand.

    startTime = sysTime();      // Time at beginning for total time display at end.
    int prevRefId = -1;      // Previous contig id.
    int prevPos = -1;
    int posOnContig = 0;
    BamAlignmentRecord record;  // Current read record.
    Dna5String contig;
    Dna5String rcContig;
    Dna5String readSeq;
    CharString readName;
    while (!atEnd(inBam))
    {
        // -------------------------------------------------------------------
        // Read SAM/BAM Record and retrieve name and sequence.
        // -------------------------------------------------------------------

        // Read SAM/BAM record.
        try
        {
            readRecord(record, inBam);
        }
        catch (seqan::ParseError const & ioErr)
        {
            std::cerr << "ERROR reading SAM/BAM record(" << ioErr.what() << ")!\n";
            return 1;
        }
        // Break if we have an unaligned SAM/BAM record.
        if (record.rID == -1)
            break;  // Done!
        // Check that the file is actually sorted by coordinate.
        if (prevRefId != -1 && (prevRefId > record.rID || (prevRefId == record.rID && prevPos > record.beginPos)))
        {
            std::cerr << "ERROR: File was not sorted by coordinate!\n";
            return 1;
        }
        // Get read name and sequence from record.
        readSeq = record.seq;  // Convert read sequence to Dna5.
        // Compute reverse complement since we align against reverse strand, SAM has aligned sequence against forward
        // strand.
        if (hasFlagRC(record))
            reverseComplement(readSeq);
        trimSeqHeaderToId(record.qName);  // Remove everything after the first whitespace.
        if (!hasFlagMultiple(record))
            append(record.qName, "/S");
        if (hasFlagMultiple(record) && hasFlagFirst(record))
            append(record.qName, "/0");
        if (hasFlagMultiple(record) && hasFlagLast(record))
            append(record.qName, "/1");
        // Translate read to read id.
        unsigned readId = 0;
        if (!getIdByName(readId, readNameStoreCache, record.qName))
        {
            readId = length(readNameStore);
            appendName(readNameStoreCache, record.qName);
            appendValue(readLengthStore, length(record.seq));
            if (options.oracleMode)
                appendValue(readAlignmentDistances, -1);
        }

        // -------------------------------------------------------------------
        // Handle Progress
        // -------------------------------------------------------------------

        // Handle progress: Change of contig.
        SEQAN_ASSERT_LEQ(prevRefId, record.rID);
        if (prevRefId != record.rID)
        {
            for (int i = prevRefId + 1; i <= record.rID; ++i)
            {
                if (i != prevRefId)
                    std::cerr << "\n";
                std::cerr << contigNameStore[record.rID] << " (" << record.rID + 1 << "/"
                          << length(contigNameStore) << ") ";
            }
            posOnContig = 0;

            // Load reference sequence on the fly using FAI index to save memory.
            unsigned faiRefId = 0;
            if (!getIdByName(faiRefId, faiIndex, contigNameStore[record.rID]))
            {
                std::cerr << "Reference sequence " << contigNameStore[record.rID] << " not known in FAI file.\n";
                return 1;
            }
            try
            {
                readSequence(contig, faiIndex, faiRefId);
            }
            catch (seqan::ParseError const & ioErr)
            {
                std::cerr << "Could not load sequence " << contigNameStore[record.rID]
                          << " from FASTA file with FAI (" << ioErr.what() << ").\n";
                return 1;
            }
            rcContig = contig;
            reverseComplement(rcContig);
        }
        // Handle progress: Position on contig.
        for (; posOnContig < record.beginPos; posOnContig += 100 * 1000)
        {
            if (posOnContig % (1000 * 1000) == 0 && posOnContig > 0)
                std::cerr << posOnContig / 1000 / 1000 << "M";
            else
                std::cerr << ".";
        }

        // -------------------------------------------------------------------
        // Extend error curve points for read.
        // -------------------------------------------------------------------

        // In oracle mode, set max error to -1, buildErrorCurvePoints() will use the error at the alignment position
        // from the SAM/BAM file.  In normal mode, convert from error rate from options to error count.
        int maxError = std::numeric_limits<int>::max();
        if (!options.oracleMode)
            maxError = static_cast<int>(floor(0.01 * options.maxError * length(record.seq)));

        // Compute end position of alignment.
        int endPos = record.beginPos + getAlignmentLengthInRef(record) - countPaddings(record.cigar);

        if (!hasFlagRC(record))
            buildErrorCurvePoints(errorCurves[readId], maxError, contig, record.rID, !hasFlagRC(record),
                                  readSeq, readId, endPos, readNameStore, options.matchN, TPatternSpec());
        else
            buildErrorCurvePoints(errorCurves[readId], maxError, rcContig, record.rID, !hasFlagRC(record),
                                  readSeq, readId, length(rcContig) - record.beginPos, readNameStore, options.matchN,
                                  TPatternSpec());
#if ENABLE
        std::cerr << "AFTER BUILD, readId == " << readId << ", read name == " << readNameStore[readId] << ", maxError == " << maxError << "\n";
#endif  // #if ENABLE
        if (options.oracleMode)
            readAlignmentDistances[readId] = maxError;

        // Update variables storing the previous read/contig id and position.
        prevRefId = record.rID;
        prevPos = record.beginPos;
    }
    std::cerr << "\n\nTook " << sysTime() - startTime << " s\n";

    // For all reads:
    //   Sort all error curve points.
    //   Fill gaps.
    //   Smooth them.
    //   Filter out low scoring ones.

    std::cerr << "\n____SMOOTHING ERROR CURVES____________________________________________________\n\n";
    startTime = sysTime();
    unsigned tenPercent = length(readLengthStore) / 10 + 1;
    std::cerr << "Progress: ";
    for (unsigned readId = 0; readId < length(readLengthStore); ++readId)
    {
        if (tenPercent > 0u && readId % tenPercent == 0u)
            std::cerr << readId / tenPercent * 10 << '%';
        else if (tenPercent > 5u && readId % (tenPercent / 5) == 0u)
            std::cerr << '.';

        std::sort(begin(errorCurves[readId], Standard()), end(errorCurves[readId], Standard()));
        fillGaps(errorCurves[readId]);
        smoothErrorCurve(errorCurves[readId]);

        // Compute relative min score for the read.
        String<WeightedMatch> filtered;
        int maxError = (int)floor(options.maxError / 100.0 * readLengthStore[readId]);
        if (options.oracleMode)
        {
            SEQAN_ASSERT_NEQ(readAlignmentDistances[readId], -1);
            maxError = readAlignmentDistances[readId];
            if (options.maxErrorSet && maxError > options.maxError)
                maxError = options.maxError;
        }
        int relativeMinScore = (int)ceilAwayFromZero(100.0 * -maxError / readLengthStore[readId]);
#if ENABLE
        std::cerr << "AFTER SMOOTHING, readId == " << readId << ", read name == " << readNameStore[readId] << ", maxError == " << maxError << "\n";
#endif  // #if ENABLE

        // Filter out low scoring ones.
        typedef typename Iterator<String<WeightedMatch> >::Type TIterator;
        for (TIterator it = begin(errorCurves[readId]); it != end(errorCurves[readId]); ++it)
        {
            if (value(it).distance >= relativeMinScore)
                appendValue(filtered, value(it));
        }
        move(errorCurves[readId], filtered);
    }
    std::cerr << "100% DONE\n"
              << "\nTook: " << sysTime() - startTime << " s\n";

    return 0;
}

// ---------------------------------------------------------------------------
// Function parseCommandLine()
// ---------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(BuildGoldStandardOptions & options, int argc, char const ** argv)
{
    // -----------------------------------------------------------------------
    // Parse Command Line Using ArgumentParser
    // -----------------------------------------------------------------------

    seqan::ArgumentParser parser("rabema_build_gold_standard");
    setShortDescription(parser, "RABEMA Gold Standard Builder");
    setVersion(parser, SEQAN_APP_VERSION " [" SEQAN_REVISION "]");
    setDate(parser, SEQAN_DATE);
    setCategory(parser, "Benchmarking");

    addUsageLine(parser,
                 "[\\fIOPTIONS\\fP] \\fB--out-gsi\\fP \\fIOUT.gsi\\fP \\fB--reference\\fP \\fIREF.fa\\fP "
                 "\\fB--in-bam\\fP \\fIPERFECT.{sam,bam}\\fP");
    addDescription(parser,
                   "This program allows one to build a RABEMA gold standard.  The input is a reference FASTA file "
                   "and a perfect SAM/BAM map (e.g. created using RazerS 3 in full-sensitivity mode).");
    addDescription(parser,
                   "The input SAM/BAM file must be \\fIsorted by coordinate\\fP.  The program will create a "
                   "FASTA index file \\fIREF.fa.fai\\fP for fast random access to the reference.");

    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable even more verbose output."));

    addSection(parser, "Input / Output");
    addOption(parser, seqan::ArgParseOption("o", "out-gsi", "Path to write the resulting GSI file to.",
                                            seqan::ArgParseArgument::OUTPUT_FILE, "GSI"));
    setValidValues(parser, "out-gsi", "gsi gsi.gz");
    setRequired(parser, "out-gsi", true);
    addOption(parser, seqan::ArgParseOption("r", "reference", "Path to load reference FASTA from.",
                                            seqan::ArgParseArgument::INPUT_FILE, "FASTA"));
    setRequired(parser, "reference", true);
    setValidValues(parser, "reference", seqan::SeqFileIn::getFileExtensions());
    addOption(parser, seqan::ArgParseOption("b", "in-bam", "Path to load the \"perfect\" SAM/BAM file from.",
                                            seqan::ArgParseArgument::INPUT_FILE, "BAM"));
    setValidValues(parser, "in-bam", BamFileIn::getFileExtensions());

    addSection(parser, "Gold Standard Parameters");
    addOption(parser, seqan::ArgParseOption("", "oracle-mode",
                                            "Enable oracle mode.  This is used for simulated data when the input "
                                            "SAM/BAM file gives exactly one position that is considered as the true "
                                            "sample position."));
    addOption(parser, seqan::ArgParseOption("", "match-N", "When set, N matches all characters without penalty."));
    addOption(parser, seqan::ArgParseOption("", "distance-metric",
                                            "Set distance metric.  Valid values: hamming, edit.  Default: edit.",
                                            seqan::ArgParseOption::STRING, "METRIC"));
    setValidValues(parser, "distance-metric", "hamming edit");
    setDefaultValue(parser, "distance-metric", "edit");

    addOption(parser, seqan::ArgParseOption("e", "max-error",
                                            "Maximal error rate to build gold standard for in percent.  This "
                                            "parameter is an integer and relative to the read length.  "
                                            "In case of oracle mode, the error rate for the read at the "
                                            "sampling position is used and \\fIRATE\\fP is used as a cutoff "
                                            "threshold.",
                                            seqan::ArgParseArgument::INTEGER, "RATE"));
    setDefaultValue(parser, "max-error", 0);

    addTextSection(parser, "Return Values");
    addText(parser, "A return value of 0 indicates success, any other value indicates an error.");

    addTextSection(parser, "Examples");

    addListItem(parser,
                "\\fBrabema_build_gold_standard\\fP \\fB-e\\fP \\fI4\\fP \\fB-o\\fP \\fIOUT.gsi\\fP \\fB-s\\fP "
                "\\fIIN.sam\\fP \\fB-r\\fP \\fIREF.fa\\fP",
                "Build gold standard from a SAM file \\fIIN.sam\\fP with all mapping locations and a FASTA "
                "reference \\fIREF.fa\\fP to GSI file \\fIOUT.gsi\\fP with a maximal error rate of \\fI4\\fP.");
    addListItem(parser,
                "\\fBrabema_build_gold_standard\\fP \\fB--distance-metric\\fP \\fIedit\\fP \\fB-e\\fP \\fI4\\fP "
                "\\fB-o\\fP \\fIOUT.gsi\\fP \\fB-b\\fP \\fIIN.bam\\fP \\fB-r\\fP \\fIREF.fa\\fP",
                "Same as above, but using Hamming instead of edit distance and BAM as the input.");
    addListItem(parser,
                "\\fBrabema_build_gold_standard\\fP \\fB--oracle-mode\\fP \\fB-o\\fP \\fIOUT.gsi\\fP \\fB-s\\fP "
                "\\fIIN.sam\\fP \\fB-r\\fP \\fIREF.fa\\fP",
                "Build gold standard from a SAM file \\fIIN.sam\\fP with the original sample position, e.g.  "
                "as exported by read simulator Mason.");

    addTextSection(parser, "Memory Requirements");
    addText(parser,
            "From version 1.1, great care has been taken to keep the memory requirements as low as possible. "
            "There memory required is two times the size of the largest chromosome plus some constant memory "
            "for each match.");
    addText(parser,
            "For example, the memory usage for 100bp human genome reads at 5% error rate was 1.7GB. Of this, "
            "roughly 400GB came from the chromosome and 1.3GB from the matches.");

    addTextSection(parser, "References");
    addText(parser,
            "M. Holtgrewe, A.-K. Emde, D. Weese and K. Reinert.  A Novel And Well-Defined Benchmarking Method "
            "For Second Generation Read Mapping, BMC Bioinformatics 2011, 12:210.");
    addListItem(parser, "\\fIhttp://www.seqan.de/rabema\\fP", "RABEMA Homepage");
    addListItem(parser, "\\fIhttp://www.seqan.de/mason\\fP", "Mason Homepage");

    // Actually do the parsing and exit on error, help display etc.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // -----------------------------------------------------------------------
    // Fill BuildGoldStandardOptions Object
    // -----------------------------------------------------------------------

    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;

    if (isSet(parser, "match-N"))
        options.matchN = true;
    if (isSet(parser, "oracle-mode"))
        options.oracleMode = true;
    if (isSet(parser, "max-error"))
    {
        getOptionValue(options.maxError, parser, "max-error");
        options.maxErrorSet = true;
    }
    CharString distanceMetric;
    getOptionValue(distanceMetric, parser, "distance-metric");
    if (distanceMetric == "hamming")
        options.distanceMetric = HAMMING_DISTANCE;

    if (isSet(parser, "out-gsi"))
        getOptionValue(options.outGsiPath, parser, "out-gsi");
    if (isSet(parser, "reference"))
        getOptionValue(options.referencePath, parser, "reference");
    if (isSet(parser, "in-bam"))
        getOptionValue(options.inBamPath, parser, "in-bam");

    return res;
}

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse command line and store results in options.
    BuildGoldStandardOptions options;
    seqan::ArgumentParser::ParseResult parseRes = parseCommandLine(options, argc, argv);
    if (parseRes != seqan::ArgumentParser::PARSE_OK)
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;

    double startTime = 0;  // For measuring time below.

    std::cerr << "==============================================================================\n"
              << "                RABEMA - Read Alignment BEnchMArk\n"
              << "==============================================================================\n"
              << "                      Building Gold Standard\n"
              << "==============================================================================\n"
              << "\n"
              << "____OPTIONS___________________________________________________________________\n\n";

    std::cerr << "Max error rate [%]    " << options.maxError << "\n"
              << "Oracle mode           " << yesNo(options.oracleMode) << '\n'
              << "Distance measure      " << metricName(options.distanceMetric) << "\n"
              << "Match Ns              " << yesNo(options.matchN) << '\n'
              << "GSI Output File       " << options.outGsiPath << '\n'
              << "SAM/BAM Input File    " << options.inBamPath << '\n'
              << "Reference File        " << options.referencePath << '\n'
              << "Verbosity             " << options.verbosity << "\n\n";

    std::cerr << "____LOADING FILES_____________________________________________________________\n\n";

    // =================================================================
    // Prepare File I/O.
    // =================================================================

    startTime = sysTime();
    std::cerr << "Reference Index       " << options.referencePath << ".fai ...";
    FaiIndex faiIndex;
    if (!open(faiIndex, toCString(options.referencePath)))
    {
        std::cerr << " FAILED (not fatal, we can just build it)\n";
        std::cerr << "Building Index        " << options.referencePath << ".fai ...";
        if (!build(faiIndex, toCString(options.referencePath)))
        {
            std::cerr << "Could not build FAI index.\n";
            return 1;
        }
        std::cerr << " OK\n";
        seqan::CharString faiPath = options.referencePath;
        append(faiPath, ".fai");
        std::cerr << "Reference Index       " << faiPath << " ...";
        try
        {
            save(faiIndex, toCString(faiPath));
        }
        catch (seqan::IOError const & ioErr)
        {
            std::cerr << "Could not write FAI index we just built (" << ioErr.what() << ").\n";
            return 1;
        }
        std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";
    }
    std::cerr << " OK (" << length(faiIndex.indexEntryStore) << " seqs)\n";

    // Open SAM file and read in header.
    BamHeader bamHeader;
    BamFileIn inBam;
    if (!open(inBam, toCString(options.inBamPath)))
    {
        std::cerr << "Could not open SAM file.\n";
        return 1;
    }
    try
    {
        readHeader(bamHeader, inBam);
    }
    catch (seqan::ParseError const & ioErr)
    {
        std::cerr << "Could not read BAM header (" << ioErr.what() << ").\n";
        return 1;
    }
    std::cerr << " OK\n";
    std::cerr << "\nTook " << sysTime() - startTime << "s\n";

    // =================================================================
    // Build point-wise error curve.
    // =================================================================
    startTime = sysTime();
    TErrorCurves errorCurves;
    int res = 0;
    StringSet<CharString> readNameStore;
    // In oracle mode, we store the distance of the alignment from the SAM file for each read.  Otherwise, this variable
    // remains unused.
    String<int> readAlignmentDistances;
    if (options.distanceMetric == EDIT_DISTANCE)
        res = matchesToErrorFunction(errorCurves, readAlignmentDistances, inBam, readNameStore,
                                     contigNames(context(inBam)), faiIndex, options, MyersUkkonenReads());
    else // options.distanceMetric == DISTANCE_HAMMING
        res = matchesToErrorFunction(errorCurves, readAlignmentDistances, inBam, readNameStore,
                                     contigNames(context(inBam)), faiIndex, options, HammingSimple());
    if (res != 0)
        return 1;

    // =================================================================
    // Convert points in error curves to intervals and write them to
    // stdout or a file.
    // =================================================================

    // The two alternatives are equivalent after opening the file.
    startTime = sysTime();
    std::cerr << "\n____WRITING OUTPUT____________________________________________________________\n\n";
    VirtualStream<char, Output> gsiStream;
    if ((options.outGsiPath == "-" && !open(gsiStream, std::cout, Nothing())) ||
        !open(gsiStream, toCString(options.outGsiPath)))
    {
        std::cerr << "Could not open output file " << options.outGsiPath << ".\n";
        return 1;
    }

    if (options.outGsiPath == "-")
        std::cerr << "Writing to stdout ...";
    else
        std::cerr << "Writing to " << options.outGsiPath << " ...";
    intervalizeAndDumpErrorCurves(gsiStream, errorCurves, readAlignmentDistances, readNameStore,
                                  contigNames(context(inBam)), options, true);
    std::cerr << " DONE\n"
              << "\n Took " << sysTime() - startTime << " s\n";

    return 0;
}
