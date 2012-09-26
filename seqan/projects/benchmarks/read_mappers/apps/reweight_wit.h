/*==========================================================================
   SeqAn - The Library for Sequence Analysis
   http://www.seqan.de 
  ==========================================================================
   Copyright (C) 2010
  
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   Lesser General Public License for more details.
  
  ==========================================================================
   Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
  ==========================================================================
   This file contains the flesh of the program reweight_wit.
  ==========================================================================*/

#ifndef BENCHMARKS_READ_MAPPERS_REWEIGHT_WIT_H_
#define BENCHMARKS_READ_MAPPERS_REWEIGHT_WIT_H_

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <cstdlib>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

#include "wit_store.h"
#include "witio.h"
#include "file_helpers.h"
#include "curve_smoothing.h"
#include "find_myers_ukkonen_ext.h"
#include "find_hamming_simple_ext.h"


struct Options
{
    // Distance function to use, also see validDistanceFunction.
    CharString distanceFunction;

    // Path to the target WIT file.
    CharString outWitFilename;

    // Path to the genome file.
    CharString genomeFilename;

    // Path to to FASTQ file with the reads.
    CharString readsFilename;

    // Path to WIT file to reweight.
    CharString inputWitFilename;

    // Maximal weighted error.
    int maxWeightedError;

    // Return true iff distanceFunction is a valid distance function.
    // Can be one of {"hamming", "edit"}.
    bool validDistanceFunction() const
    {
        if (distanceFunction == "hamming") return true;
        if (distanceFunction == "edit") return true;
        return false;
    }
};


// Compute quality-based alignment score.  The read has to be given
// since we do not have qualities in the alignment object.
template <typename TAlign>
int computeQualityAlignmentScore(TAlign const & align, Score<int, ScoreMatrix<Dna5> > const & scoreMatrix, String<Dna5Q> const & read) {
    // TODO(holtgrew): Maybe convert to iterators for performance?
    typedef typename Row<TAlign const>::Type TRow;
    typedef typename Value<TRow>::Type TAlignChar;
    typedef typename Size<TRow>::Type TSize;

    TRow & rowContig = row(align, 0);
    TRow & rowRead = row(align, 1);

    int result = 0;
    for (TSize i = 0; i < length(rowContig); ++i) {
        TAlignChar contigChar = rowContig[i];
        TAlignChar readChar = rowRead[i];
        if (isGap(rowContig, i)) {
            result -= getQualityValue(read[toSourcePosition(rowRead, i)]);
        } else if (isGap(rowRead, i)) {
            if (toSourcePosition(rowRead, i) == 0) {
                result -= getQualityValue(read[0]);
            } else if (toSourcePosition(rowRead, i) == length(read)) {
                result -= getQualityValue(read[length(read) - 1]);
            } else {
                int x = 0;
                x += getQualityValue(read[toSourcePosition(rowRead, i) - 1]);
                x += getQualityValue(read[toSourcePosition(rowRead, i)]);
                result -= ceil(1.0 * x / 2);
            }
        } else {
            result += score(scoreMatrix, readChar, contigChar) * getQualityValue(read[toSourcePosition(rowRead, i)]);
        }
    }

    return result;
}


// Build intervals from the error curves and append resulting intervals to store.
void intervalizeErrorCurves(WitStore & store,
                            String<WeightedMatch> const & matches,
                            IntervalOfReadOnContig const & originalInterval,
                            int const minScore) {
    typedef String<WeightedMatch> TWeightedMatches;
    size_t readId = originalInterval.readId;

    // Sort the matches.  Matches with high scores (negative score and
    // low absolute value) come first.
    TWeightedMatches sortedMatches(matches);
    std::sort(begin(sortedMatches, Standard()), end(sortedMatches, Standard()));
    
    // intervals[e] holds the intervals for error e of the current read.
    String<String<IntervalOfReadOnContig> > intervals;
    resize(intervals, -minScore + 1, Exact());

    // Join the intervals stored in sortedMatches.
    //
    // The position of the previous match, so we can consider only the
    // ones with the smallest error.
    //
    // The following two vars should be != first pos and contigId.
    size_t previousPos = maxValue<size_t>();
    size_t previousContigId = maxValue<size_t>();
    typedef Iterator<TWeightedMatches>::Type TWeightedMatchesIter;
    for (TWeightedMatchesIter it = begin(sortedMatches);
         it != end(sortedMatches); ++it) {
        // Skip it if (it - 1) pointed to same pos (and must point to
        // one with smaller absolute distance.
        if (it->pos == previousPos && it->contigId == previousContigId)
            continue;
        // Consider all currently open intervals with a greater error
        // than the error in *it and extend them or create a new one.
        for (int e = originalInterval.distance; e <= -minScore; ++e) {
            // Handle base case of no open interval:  Create new one.
            if (length(intervals[e]) == 0) {
                appendValue(intervals[e], IntervalOfReadOnContig(readId, e, it->contigId, it->isForward, it->pos, it->pos));
                continue;
            }
            IntervalOfReadOnContig & interval = back(intervals[e]);
            // Either extend the interval or create a new one.
            SEQAN_ASSERT(interval.last <= pos);
            if (interval.lastPos + 1 == it->pos)
                back(intervals[e]).lastPos += 1;
            else
                appendValue(intervals[e], IntervalOfReadOnContig(readId, e, it->contigId, it->isForward, it->pos, it->pos));
        }
        // Book-keeping.
        previousPos = it->pos;
        previousContigId = it->contigId;
    }

    // Append the resulting intervals to the result store.
    for (size_t i = 0; i < length(intervals); ++i)
        for (size_t j = 0; j < length(intervals[i]); ++j)
            appendValue(store, intervals[i][j]);
}


template <typename TPatternSpec>
void reweightInterval(WitStore & store,
                      IntervalOfReadOnContig const & interval,
                      String<Dna5> /*const*/ & contig,
                      String<Dna5Q> /*const*/ & read,
                      Options const & options,
                      TPatternSpec const &) {
    Finder<String<Dna5> > finder(contig);
    Pattern<String<Dna5Q>, TPatternSpec> pattern(read, -length(read));
    typedef typename Size<String<Dna5> >::Type TSize;

    // We are only considering the weighted case, N is a wildcard.
    _patternMatchNOfPattern(pattern, true);
    _patternMatchNOfFinder(pattern, true);

    // Build scoring matrix that allows N to match with all.
    int gapExtensionScore = -1;
    int gapOpenScore = -1;
    if (IsSameType<TPatternSpec, HammingSimple>::VALUE) {
        // No gaps for hamming distance.
        gapOpenScore = -length(read);
        gapExtensionScore = -length(read);
    }
    Score<int, ScoreMatrix<Dna5> > matrixScore(gapExtensionScore, gapOpenScore);
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x) {
        for (int y = 0; y < ValueSize<Dna5>::VALUE - 1; ++y)
            setScore(matrixScore, Dna5(x), Dna5(y), -1);
        setScore(matrixScore, Dna5(x), Dna5('N'), 0);
        setScore(matrixScore, Dna5('N'), Dna5(x), 0);
        setScore(matrixScore, Dna5(x), Dna5(x), 0);
    }
    for (int x = 0; x < ValueSize<Dna5>::VALUE; ++x)
        setScore(matrixScore, Dna5('N'), Dna5(x), 0);

    bool ret = setEndPosition(finder, pattern, interval.firstPos);
    SEQAN_ASSERT(ret);
    ret = findBegin(finder, pattern, getScore(pattern));
    SEQAN_ASSERT(ret);

    String<WeightedMatch> weightedMatches;

    while (find(finder, pattern) && endPosition(finder) <= interval.lastPos + 1) {
        bool ret = findBegin(finder, pattern, getScore(pattern));
        (void)ret;  // Supress warnings in Release mode.
        SEQAN_ASSERT(ret);
        SEQAN_ASSERT_GEQ(static_cast<int>(1.0 * getScore(pattern) / length(read)), -static_cast<int>(interval.distance));

        // Prepare alignment datastructures.
        Align<String<Dna5>, ArrayGaps> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(finder));
        assignSource(row(align, 1), read);

        // Perform banded Needleman-Wunsch alignment.
        StringSet<String<Dna5> > stringSet;
        appendValue(stringSet, infix(finder));
        appendValue(stringSet, read);
        int alignmentScore = globalAlignment(align, stringSet, matrixScore, 2*getScore(pattern), -2*getScore(pattern), BandedNeedlemanWunsch());
        (void)alignmentScore;  // Supress warnings in Release mode.
        SEQAN_ASSERT_EQ(alignmentScore, getScore(pattern));

        // Compute quality-based score of alignment.  We pass the
        // score matrix to allow for N-is-wildcard mode.
        int qualityValue = computeQualityAlignmentScore(align, matrixScore, read);

        // Compute relative quality score and append weighted match to list of
        // weighted matches.
        int relativeQualityScore = -static_cast<int>(ceil(-100.0 * qualityValue / length(read)));
        // TODO(holtgrew): Reorder begin pos in WeightedMatch constructor.
        WeightedMatch weightedMatch(interval.contigId, interval.isForward, endPosition(finder) - 1, relativeQualityScore, beginPosition(finder));
        appendValue(weightedMatches, weightedMatch);
    }

    // Smooth the error curve to get rid of falsely separating bad
    // scores.  Note that we do not need to fill gaps here.
    smoothErrorCurve(weightedMatches);

    // Build intervals from data points and append to result WitStore.
    intervalizeErrorCurves(store, weightedMatches, interval, -options.maxWeightedError);
}


template <typename TPatternSpec>
void reweightWitStoreForContig(WitStore & reweightedStore,
                               WitStore /*const*/ & store,
                               StringSet<String<Dna5> > /*const*/ & contigs,
                               StringSet<String<Dna5Q> > /*const*/ & reads,
                               size_t contigId,
                               Options const & options,
                               TPatternSpec const &)
{
    // TODO(holtgrew): Change to iterate over intervals of contig ids, as in wit builder.  The way it is done here with filtering is slow!
    
    typedef typename WitStore::TIntervalStore TIntervalStore;
    typedef typename Iterator<TIntervalStore, Standard>::Type TIntervalIterator;

    String<Dna5> rcContig = contigs[contigId];
    reverseComplement(rcContig);
    
    for (TIntervalIterator it = begin(store.intervals, Standard()) + 1; it != end(store.intervals, Standard()); ++it) {
        if (value(it - 1).readId == value(it).readId &&
            value(it - 1).contigId == value(it).contigId &&
            value(it - 1).firstPos == value(it).firstPos) {
            continue;  // Skip all of non-maximal weight.
        }
        IntervalOfReadOnContig interval = value(it - 1);
        if (interval.contigId == contigId) {
            if (interval.isForward)
                reweightInterval(reweightedStore, interval, contigs[interval.contigId], reads[interval.readId], options, TPatternSpec());
            else
                reweightInterval(reweightedStore, interval, rcContig, reads[interval.readId], options, TPatternSpec());
        }
    }
    // Do not forget last interval.
    if (back(store.intervals).contigId == contigId) {
        if (back(store.intervals).isForward)
            reweightInterval(reweightedStore, back(store.intervals), contigs[back(store.intervals).contigId],
                             reads[back(store.intervals).readId], options, TPatternSpec());
        else
            reweightInterval(reweightedStore, back(store.intervals), rcContig,
                             reads[back(store.intervals).readId], options, TPatternSpec());
    }
}


template <typename TPatternSpec>
void reweightWitStore(WitStore & store,
                      StringSet<String<Dna5> > /*const*/ & contigs,
                      StringSet<String<Dna5Q> > /*const*/ & reads,
                      Options const & options,
                      TPatternSpec const &)
{
    SEQAN_ASSERT_GT(length(store.intervals), 0u);

    // Sort intervals in store by (contigId, readId, distance).
    sortWitRecords(store, SortDistance());
    sortWitRecords(store, SortFirstPos());
    sortWitRecords(store, SortReadId());
    sortWitRecords(store, SortContigId());

    // We will insert the reweighted intervals into reweightedStore and later
    // move this back into store again.
    WitStore reweightedStore;
    reweightedStore.readNames = store.readNames;
    reweightedStore.contigNames = store.contigNames;
    for (size_t contigId = 0; contigId < length(contigs); ++contigId) {
        reweightWitStoreForContig(reweightedStore, store, contigs, reads, contigId, options, TPatternSpec());
    }
    move(store, reweightedStore);
}

#endif  // BENCHMARKS_READ_MAPPERS_REWEIGHT_WIT_H_
