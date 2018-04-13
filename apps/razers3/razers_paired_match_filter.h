#ifndef RAZERS_PAIRED_MATCH_FILTER_H_
#define RAZERS_PAIRED_MATCH_FILTER_H_

#include <seqan/graph_types/graph_idmanager.h>

#include "razers.h"

namespace seqan {

template <typename TThreadLocalStorage>
class FilterPatternLSetMaxErrorsWrapper;
template <typename TThreadLocalStorage>
class FilterPatternRSetMaxErrorsWrapper;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
class PairedMatchFilter
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
    std::unordered_map<unsigned, unsigned> pairIdToHistogramId;
    // Id manager for histogram.
    IdManager<unsigned> idManager;
    // Histograms.
    String<String<unsigned> > histograms;
    // Ids of the reads that are purged.
    String<unsigned> purgedPairIds;
    // The callback context object.
    Holder<TCallback> callback;
    // Read ID offset.
    unsigned readOffset;
    // The read sequences.
    TReadSeqSet const & readSeqs;
    // The options object.
    RazerSCoreOptions<TOptionsSpec> const & options;

    PairedMatchFilter(unsigned readCount_, unsigned matchThreshold_, double frac_, TCallback & callback_, unsigned readOffset_, TReadSeqSet const & readSeqs_, RazerSCoreOptions<TOptionsSpec> const & options_) :
        readCount(readCount_), matchThreshold(matchThreshold_), frac(frac_), callback(callback_), readOffset(readOffset_), readSeqs(readSeqs_), options(options_)
    {
        resize(hitCount, readCount, 0, Exact());
        reserve(histograms, unsigned(frac * readCount));
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline unsigned
_createHistogram(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId)
{
    unsigned result = obtainId(filter.idManager);
    if (result >= length(filter.histograms))
        resize(filter.histograms, length(filter.histograms) + 1);
    resize(filter.histograms[result], 1 + (int)(filter.options.errorRate * (length(filter.readSeqs[2 * pairId]) + length(filter.readSeqs[2 * pairId + 1]))), 0, Exact());
    SEQAN_ASSERT_LT(result, length(filter.histograms));
    return result;
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline void
registerRead(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId, int score)
{
    // std::cerr << "registering read " << pairId << std::endl;
    if (filter.hitCount[pairId - filter.readOffset] == std::numeric_limits<unsigned>::max())
        return;

    if (++filter.hitCount[pairId - filter.readOffset] < filter.matchThreshold)
        return;

    // TODO(holtgrew): Maybe global read to histogram map; faster?

    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    unsigned histogramId = 0;
    if (filter.hitCount[pairId - filter.readOffset] == filter.matchThreshold)
    {
        // std::cerr << "new histogram for read " << pairId << std::endl;
        histogramId = _createHistogram(filter, pairId);
        filter.pairIdToHistogramId[pairId] = histogramId;
    }
    else
    {
        // std::cerr << "updating histogram for read " << pairId << std::endl;
        typedef typename std::unordered_map<unsigned, unsigned>::iterator TIterator;
        TIterator it = filter.pairIdToHistogramId.find(pairId);
        SEQAN_ASSERT(it != filter.pairIdToHistogramId.end());
        histogramId = it->second;
    }

    // Insert value into histogram.
    _incrementCount(filter, histogramId, score);
}

template <typename TOptionsSpec, typename TReadSeqSet, typename TCallback>
inline bool
processRead(PairedMatchFilter<TOptionsSpec, TReadSeqSet, TCallback> & filter, unsigned pairId)
{
    typedef typename std::unordered_map<unsigned, unsigned>::iterator TIterator;

    if (filter.hitCount[pairId - filter.readOffset] < filter.matchThreshold)
        return false;

    // std::cerr << "processing read " << pairId << std::endl;
    // Get histogram id, insert new histogram if necessary, exit if no histogram yet.
    TIterator it = filter.pairIdToHistogramId.find(pairId);
    if (it == filter.pairIdToHistogramId.end())
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
//            std::cerr << "PURGED " << pairId << std::endl;
            appendValue(filter.purgedPairIds, pairId);
            newLimit = 0;
        }
    }
    else
    {
        newLimit = _newLimitDistRange(filter, histogramId);
        if (newLimit == CAN_BE_PURGED)
        {
            // std::cerr << "PURGED " << pairId << std::endl;
            appendValue(filter.purgedPairIds, pairId);
            newLimit = 0;
        }
    }

//    std::cerr << "LIMITING " << pairId << "\t" << filter.histograms[histogramId][0] << "\t" << filter.hitCount[pairId - filter.readOffset] << "\t" << newLimit << std::endl;
    FilterPatternLSetMaxErrorsWrapper<TCallback> wrapperL(value(filter.callback));
    FilterPatternRSetMaxErrorsWrapper<TCallback> wrapperR(value(filter.callback));
    setMaxErrors(wrapperL, pairId, newLimit - 1);
    setMaxErrors(wrapperR, pairId, newLimit - 1);
    if (filter.options.errorCutOff[2 * pairId] > newLimit)
        filter.options.errorCutOff[2 * pairId] = newLimit;
    if (filter.options.errorCutOff[2 * pairId + 1] > newLimit)
        filter.options.errorCutOff[2 * pairId + 1] = newLimit;

    if (newLimit == 0)
    {
        _freeHistogram(filter, histogramId);
        filter.pairIdToHistogramId.erase(pairId);
        filter.hitCount[pairId - filter.readOffset] = std::numeric_limits<unsigned>::max();
        return true;
    }
    return false;
}

}  // namespace seqan

#endif  // #ifndef RAZERS_PAIRED_MATCH_FILTER_H_
