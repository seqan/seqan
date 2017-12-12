#ifndef RAZERS_MATCH_FILTER_H_
#define RAZERS_MATCH_FILTER_H_

#include <seqan/graph_types/graph_idmanager.h>

#include "razers.h"

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
class MatchFilter
{
public:
    // Number of reads.
    unsigned readCount;
    // If a read has more than matchThreshold matches then it gets a histogram.
    unsigned matchThreshold;
    // Fraction of reads expect to have histograms for.
    double frac;
    // Count matches for each read
    String<unsigned> hitCount;
    // Map from read number to histogram id.
    std::unordered_map<unsigned, unsigned> readIdToHistogramId;
    // Id manager for histogram.
    IdManager<unsigned> idManager;
    // Histograms.
    String<String<unsigned> > histograms;
    // Ids of the reads that are purged.
    String<unsigned> purgedReadIds;
    // The callback context object.
    Holder<TCallback> callback;
    // Read ID offset.
    unsigned readOffset;
    // The read sequences.
    TReadSeqSet const & readSeqs;
    // The options object.
    RazerSCoreOptions<TOptionsSpec> const & options;

    MatchFilter(unsigned readCount_, unsigned matchThreshold_, double frac_, TCallback & callback_, unsigned readOffset_, TReadSeqSet const & readSeqs_, RazerSCoreOptions<TOptionsSpec> const & options_) :
        readCount(readCount_), matchThreshold(matchThreshold_), frac(frac_), callback(callback_), readOffset(readOffset_), readSeqs(readSeqs_), options(options_)
    {
        resize(hitCount, readCount, 0, Exact());
        reserve(histograms, unsigned(frac * readCount));
    }

};

const int NO_NEW_LIMIT = -1;
const int CAN_BE_PURGED = -2;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline unsigned
_createHistogram(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId)
{
    unsigned result = obtainId(filter.idManager);
    if (result >= length(filter.histograms))
        resize(filter.histograms, length(filter.histograms) + 1);
    resize(filter.histograms[result], 1 + (int)(filter.options.errorRate * length(filter.readSeqs[readId])), 0, Exact());
    SEQAN_ASSERT_LT(result, length(filter.histograms));
    return result;
}

template <typename TMatchFilter>
inline void
_freeHistogram(TMatchFilter & filter, unsigned histogramId)
{
    clear(filter.histograms[histogramId]);
    releaseId(filter.idManager, histogramId);
}

template <typename TMatchFilter>
inline void
_incrementCount(TMatchFilter & filter, unsigned histogramId, int score)
{
    SEQAN_ASSERT_LEQ(score, 0);
    ++filter.histograms[histogramId][-score];
}

template <typename TMatchFilter>
inline int
_newLimit(TMatchFilter const & filter, unsigned histogramId)
{
    // TODO(holtgrew): This could be speeded up if using prefix sum data structures for histograms.
    typedef typename Iterator<String<unsigned> const, Standard>::Type TIter;

    int max = filter.options.maxHits;
    if (max == 0)
        return NO_NEW_LIMIT;

    // we need at least one more than max before we can purge
    if (filter.options.purgeAmbiguous)
        ++max;

    String<unsigned> const & hist = filter.histograms[histogramId];
    TIter itBeg = begin(hist, Standard());
    TIter itEnd = end(hist, Standard());

    for (TIter it = itBeg; it != itEnd; ++it)
    {
        max -= *it;
        if (max <= 0)
            return it - itBeg;
    }
    return NO_NEW_LIMIT;
}

template <typename TMatchFilter>
inline int
_newLimitDistRange(TMatchFilter const & filter, unsigned histogramId)
{
    // TODO(holtgrew): This could be speeded up if using prefix sum data structures for histograms.
    typedef typename Iterator<String<unsigned> const, Standard>::Type TIter;

    String<unsigned> const & hist = filter.histograms[histogramId];

    int max = filter.options.maxHits;
    TIter itBeg = begin(hist, Standard());
    TIter itEnd = end(hist, Standard());
    TIter it = itBeg;

    if (max == 0)
    {
        // adapt maximal errors
        while (it != itEnd && *it == 0)
            ++it;

        return _min((it - itBeg) + (int)filter.options.scoreDistanceRange, itEnd - itBeg);
    }

    // we need at least one more than max before we can purge
    if (filter.options.purgeAmbiguous)
        ++max;

    while (it != itEnd && *it == 0)
        ++it;

    // look in a range [minErrors..minErrors+scoreDistanceRange)
    if (itEnd > it + filter.options.scoreDistanceRange)
        itEnd = it + filter.options.scoreDistanceRange;

    for (; it != itEnd; ++it)
    {
        max -= *it;
        if (max <= 0)
        {
            // do we have enough matches in [0..scoreDistanceRange) ?
            if (filter.options.purgeAmbiguous && it - itBeg < (int)filter.options.scoreDistanceRange)
                return CAN_BE_PURGED;
            else
                return it - itBeg;
        }
    }
    return it - itBeg;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
registerRead(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId, int score)
{
    // std::cerr << "registering read " << readId << std::endl;
    if (filter.hitCount[readId - filter.readOffset] == std::numeric_limits<unsigned>::max())
        return;

    if (++filter.hitCount[readId - filter.readOffset] < filter.matchThreshold)
        return;

    // TODO(holtgrew): Maybe global read to histogram map; faster?

    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    unsigned histogramId = 0;
    if (filter.hitCount[readId - filter.readOffset] == filter.matchThreshold)
    {
        // std::cerr << "new histogram for read " << readId << std::endl;
        histogramId = _createHistogram(filter, readId);
        filter.readIdToHistogramId[readId] = histogramId;
    }
    else
    {
        // std::cerr << "updating histogram for read " << readId << std::endl;
        typedef typename std::unordered_map<unsigned, unsigned>::iterator TIterator;
        TIterator it = filter.readIdToHistogramId.find(readId);
        SEQAN_ASSERT(it != filter.readIdToHistogramId.end());
        histogramId = it->second;
    }

    // Insert value into histogram.
    _incrementCount(filter, histogramId, score);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline bool
processRead(MatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned readId)
{
    typedef typename std::unordered_map<unsigned, unsigned>::iterator TIterator;

    if (filter.hitCount[readId - filter.readOffset] < filter.matchThreshold)
        return false;

    // std::cerr << "processing read " << readId << std::endl;
    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    TIterator it = filter.readIdToHistogramId.find(readId);
    if (it == filter.readIdToHistogramId.end())
        return false;  // Must have been disabled before.

    unsigned histogramId = it->second;

    // Perform actions.
    int newLimit;

    if (filter.options.scoreDistanceRange == 0)
    {
        newLimit = _newLimit(filter, histogramId);
        if (newLimit == NO_NEW_LIMIT)
            return false;

        if (filter.options.purgeAmbiguous)
        {
//            std::cerr << "PURGED " << readId << std::endl;
            appendValue(filter.purgedReadIds, readId);
            newLimit = 0;
        }
    }
    else
    {
        newLimit = _newLimitDistRange(filter, histogramId);
        if (newLimit == CAN_BE_PURGED)
        {
            // std::cerr << "PURGED " << readId << std::endl;
            appendValue(filter.purgedReadIds, readId);
            newLimit = 0;
        }
    }

//    std::cerr << "LIMITING " << readId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[readId - filter.readOffset] << "\t" << newLimit << std::endl;
    limitRead(value(filter.callback), readId, newLimit - 1);
    filter.options.errorCutOff[readId] = newLimit;
    if (newLimit == 0)
    {
        _freeHistogram(filter, histogramId);
        filter.readIdToHistogramId.erase(readId);
        filter.hitCount[readId - filter.readOffset] = std::numeric_limits<unsigned>::max();
        return true;
    }
    return false;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_MATCH_FILTER_H_
