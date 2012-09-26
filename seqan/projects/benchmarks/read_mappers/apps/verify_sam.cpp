/*
  Simple "verifier" for Sam files that checks whether all matches of
  reads in a contig are actually fund.

  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
*/

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>

#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

#include "curve_smoothing.h"
#include "witio.h"

using namespace seqan;

char const * kUsageString = "%s CONTIGS.fasta MAP.sam READ-ID\n";


// Represents an interval in a contig.
struct ContigInterval {
    // Contig the interval is in.
    size_t contigId;
    // First value of interval
    size_t first;
    // Last value of interval.
    size_t last;
    
    ContigInterval()
    {}

    ContigInterval(size_t _contigId, size_t _first, size_t _last) :
        contigId(_contigId), first(_first), last(_last)
    {}
};


template <typename Stream>
Stream &operator<<(Stream &out, const ContigInterval &c) {
    out << "ContigInterval(contigId=" << c.contigId << ", first="
        << c.first << ", last=" << c.last << ")";
    return out;
}


// Map the given read to the contigs in the fragment store using an
// online search algorithm.
template <typename TString, typename TFragmentStore>
void mapReadOnline(TWeightedMatches & matches,
                   TString const & read, TFragmentStore & fragments,
                   double maxErrorRate) {
    typedef typename TFragmentStore::TContigSeq        TContigSeq;
    TString readCopy(read);
    int minScore = (int)-floor(maxErrorRate / 100.0 * length(read));
    std::cout << "  min score == " << minScore << std::endl;
    for (size_t contigId = 0; contigId < length(fragments.contigStore); ++contigId) {
        std::cout << "  contig " << fragments.contigNameStore[contigId] << " (#" << contigId << ")..." << std::endl;
        TContigSeq & contig = fragments.contigStore[contigId].seq;
        // Find on forward strand.
        {
            Finder<TContigSeq> finder(contig);
            Pattern<TString, Myers<FindInfix> > pattern(readCopy, minScore);
            while (find(finder, pattern)) {
                bool ret = findBegin(finder, pattern);
                SEQAN_ASSERT(ret);
                int theScore = getScore(pattern);
                appendValue(matches, WeightedMatch(contigId, true, endPosition(finder), theScore, beginPosition(finder)));
            }
        }
        // Find on backward strand.
        {
            reverseComplement(readCopy);
            Finder<TContigSeq> finder(contig);
            Pattern<TString, Myers<FindInfix> > pattern(readCopy, minScore);
            while (find(finder, pattern)) {
                bool ret = findBegin(finder, pattern);
                SEQAN_ASSERT(ret);
                int theScore = getScore(pattern);
                appendValue(matches, WeightedMatch(contigId, false, endPosition(finder), theScore, beginPosition(finder)));
            }
        }
    }
    std::sort(begin(matches, Standard()), end(matches, Standard()));
    fillGaps(matches);
}


// Convert point-wise list of matches to an interval list.
template <typename TFragmentStore>
void intervalizeMatches(String<std::map<size_t, std::pair<size_t, bool> > > & intervalMaps,
                        TWeightedMatches & matches, TFragmentStore const & fragments) {
    typedef size_t TSize;
    TSize previousPos = maxValue<TSize>();
    TSize previousContigId = maxValue<TSize>();
    typedef Iterator<TWeightedMatches>::Type TWeightedMatchesIter;
    String<ContigInterval> intervals;
    for (TWeightedMatchesIter it = begin(matches); it != end(matches); ++it) {
        if (it->pos == previousPos && it->contigId == previousContigId)
            continue;
        if (length(intervals) == 0) {
            appendValue(intervals, ContigInterval(it->contigId, it->pos, it->pos));
            continue;
        }
        if (back(intervals).last + 1 == it->pos) {
            back(intervals).last += 1;
        } else {
            appendValue(intervals, ContigInterval(it->contigId, it->pos, it->pos));
        }
        previousPos = it->pos;
        previousContigId = it->contigId;
    }

    // Build result.
    resize(intervalMaps, length(fragments.contigNameStore)); 
    for (size_t i = 0; i < length(intervalMaps); ++i)
        intervalMaps[i].clear();
    for (size_t i = 0; i < length(intervals); ++i) {
        std::cout << "  Adding interval [" << intervals[i].first << ", " << intervals[i].last << "]"
                  << " on contig " << fragments.contigNameStore[intervals[i].contigId] << " (#"
                  << intervals[i].contigId << ")" << std::endl;
        intervalMaps[intervals[i].contigId][intervals[i].last] = std::make_pair(intervals[i].first, false);
    }
}


template <typename TFragmentStore>
void checkMatches(String<std::map<size_t, std::pair<size_t, bool> > > & intervalMaps,
                  TFragmentStore & fragments, size_t readId) {
    typedef std::map<size_t, std::pair<size_t, bool> > TMap;

    typedef typename TFragmentStore::TContigStore      TContigStore;
    typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
    typedef typename Value<TContigStore>::Type         TContig;
    typedef typename TFragmentStore::TContigSeq        TContigSeq;
    typedef typename TFragmentStore::TReadSeq          TReadSeq;
    typedef typename Value<TAlignedReadStore>::Type    TAlignedRead;
    typedef typename TAlignedRead::TPos                TAlignedReadPos;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> > TContigGaps;

    bool foundError = false;
    size_t alignedReadCount = 0;
    for (size_t i = 0; i < length(fragments.alignedReadStore); ++i) {
        if (fragments.alignedReadStore[i].readId != readId) continue;
        alignedReadCount += 1;
        size_t contigId = fragments.alignedReadStore[i].contigId;

//         TContigSeq & contig = fragments.contigStore[contigId].seq;
        TContigGaps contigGaps(fragments.contigStore[contigId].seq, fragments.contigStore[contigId].gaps);

        size_t lastGapPos = _max(fragments.alignedReadStore[i].beginPos, fragments.alignedReadStore[i].endPos);
        size_t firstGapPos = _min(fragments.alignedReadStore[i].beginPos, fragments.alignedReadStore[i].endPos);

        size_t last = positionGapToSeq(contigGaps, lastGapPos);
        size_t first = positionGapToSeq(contigGaps, firstGapPos);

        TMap::iterator it;
        it = intervalMaps[contigId].lower_bound(last);
        if (it == intervalMaps[contigId].end()) {
            std::cout << "ERROR, found no interval for aligned read!" << std::endl
                      << "  contig id = " << contigId << std::endl
                      << "  name is     " << fragments.contigNameStore[contigId] << std::endl
                      << "  first pos = " << first << std::endl
                      << "  last pos =  " << last << std::endl
                      << "  begin pos = " << fragments.alignedReadStore[i].beginPos << std::endl
                      << "  end pos =   " << fragments.alignedReadStore[i].endPos << std::endl;
            foundError = true;
            continue;
        }
        if (it->second.first > last) {
            std::cout << "ERROR, found WRONG interval for aligned read!" << std::endl
                      << "  Aligned read is" << std::endl
                      << "    contig id = " << contigId << std::endl
                      << "    name is     " << fragments.contigNameStore[contigId] << std::endl
                      << "    begin pos = " << fragments.alignedReadStore[i].beginPos << std::endl
                      << "    end pos =   " << fragments.alignedReadStore[i].endPos << std::endl
                      << "  Interval is [" << it->second.first << ", " << it->first << "]" << std::endl;
            foundError = true;
            continue;
        }
        it->second.second = true;
    }
    std::cout << "Total number of alignments for read " << readId << " in Sam file: " << alignedReadCount << std::endl;

    typedef std::map<size_t, std::pair<size_t, bool> > TMap;
    size_t totalIntervals = 0;
    size_t intervalsHit = 0;
    for (size_t contigId = 0; contigId < length(intervalMaps); ++contigId) {
        TMap const & intervalMap = intervalMaps[contigId];
        totalIntervals += intervalMap.size();
        for (TMap::const_iterator it = intervalMap.begin(); it != intervalMap.end(); ++it) {
            if (!it->second.second) {
                foundError = true;
                std::cout << "ERROR, the following interval in contig "
                          << fragments.contigNameStore[contigId]
                          << " (#" << contigId << ") was not hit: ["
                          << it->second.first << ", " << it->first
                          << "]" << std::endl;
            } else {
                intervalsHit += 1;
            }
        }
    }

    std::cout << "Intervals hit: " << intervalsHit << " / " << totalIntervals << std::endl;
    
    if (!foundError)
        std::cout << "Sam file contains all hits and no errorneous ones." << std::endl;
    else
        std::cout << "Sam file NOT CORRECT (see above)." << std::endl;
}


int main(int argc, const char *argv[]) {
    // Verify command line arguments.
    if (argc != 4) {
        fprintf(stderr, "ERROR: Invalid number of arguments!\n");
        fprintf(stderr, kUsageString, argv[0]);
        return 1;
    }
    // Get shortcuts to command line arguments.
    char const * kContigsPath = argv[1];
    char const * kMapsPath = argv[2];

    // Load fragment store from file.
    typedef FragmentStore<> TFragmentStore;
    TFragmentStore fragments;
    // Load contigs and Sam file.  Reads are already in the Sam file.
    std::cout << "Loading fragment store..." << std::endl;
    if (!loadContigs(fragments, kContigsPath)) {
        std::cerr << "Could not read contigs from file " << kContigsPath << std::endl;
        return 1;
    }
    // Read Sam file.
    {
        std::fstream fstrm(kMapsPath, std::ios_base::in | std::ios_base::binary);
        if (!fstrm.is_open()) {
            std::cerr << "Could not open Sam file." << std::endl;
            return 1;
        }
        read(fstrm, fragments, Sam());
    }

    // Search reads in contigs.
    const size_t readId = atoi(argv[3]);
    const double maxErrorRate = 10.0;  // TODO(holtgrew): Make configurable through args.
    std::cout << "Searching for reads in contig..." << std::endl;
    std::cout << "  read id =   " << readId << std::endl;
    std::cout << "  read name = " << fragments.readNameStore[readId] << std::endl;
    std::cout << "  read seq =  " << fragments.readSeqStore[readId] << std::endl;
    Value<TFragmentStore::TReadSeqStore>::Type rcRead(fragments.readSeqStore[readId]);
    reverseComplement(rcRead);
    std::cout << "  rev comp =  " << rcRead << std::endl;
    TWeightedMatches matches;
    mapReadOnline(matches, fragments.readSeqStore[readId], fragments, maxErrorRate);

    // Smoothing error curve.
    std::cout << "Smoothing error curve..." << std::endl;
    smoothErrorCurve(matches);

    std::cout << "--" << std::endl;
    for (size_t i = 0; i < length(matches); ++i) {
        std::cout << matches[i] << std::endl;
    }
    std::cout << "--" << std::endl;

    // Build intervals from error curve.
    std::cout << "Building intervals from aligned read positions..." << std::endl;
    // Maps right interval border to pair with (left border, (weight, has been hit flag)).
    String<std::map<size_t, std::pair<size_t, bool> > > intervalMaps;
    intervalizeMatches(intervalMaps, matches, fragments);

    // Check whether RazerS hits all lakes.
    std::cout << "Checking whether all lakes are hit by Sam file..." << std::endl;
    checkMatches(intervalMaps, fragments, readId);

    return 0;
}
