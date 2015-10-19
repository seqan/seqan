// ==========================================================================
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================

#ifndef GENERALSTATS_H
#define GENERALSTATS_H


struct AdapterTrimmingStats
{
    std::vector<std::vector<unsigned>> removedLength;
    std::vector<unsigned> numRemoved;
    unsigned overlapSum;
    unsigned minOverlap, maxOverlap;

    AdapterTrimmingStats() : overlapSum(0),
        minOverlap(std::numeric_limits<unsigned>::max()), maxOverlap(0) {};

    AdapterTrimmingStats& operator+= (AdapterTrimmingStats const& rhs)
    {
        overlapSum += rhs.overlapSum;
        minOverlap = minOverlap < rhs.minOverlap ? minOverlap : rhs.minOverlap;
        maxOverlap = maxOverlap < rhs.maxOverlap ? rhs.maxOverlap : maxOverlap;
        {
            const auto len = rhs.removedLength.size();
            if (removedLength.size() < len)
                removedLength.resize(std::max(removedLength.size(), len));
            for (unsigned int i = 0;i < len;++i)
            {
                const auto len2 = rhs.removedLength[i].size();
                if (removedLength[i].size() < len2)
                    removedLength[i].resize(len2);
                for (unsigned k = 0;k < len2;++k)
                    removedLength[i][k] += rhs.removedLength[i][k];
            }
        }

        {
            const auto len = rhs.numRemoved.size();
            if (numRemoved.size() < len)
                numRemoved.resize(len);
            for (unsigned int i = 0;i < len;++i)
                numRemoved[i] += rhs.numRemoved[i];
        }
        return *this;
    }
    void clear()
    {
        overlapSum = 0;
        minOverlap = std::numeric_limits<unsigned>::max();
        maxOverlap = 0;
    }
};

struct GeneralStats
{
    unsigned removedN;       //Number of deleted sequences due to N's
    unsigned removedDemultiplex;
    unsigned removedQuality;
    unsigned long uncalledBases;//Number of uncalled bases (evtl. Masked) in surviving sequences
    unsigned removedShort;  //Number of deleted sequences due to shortness.
    unsigned int readCount;
    double processTime;
    double ioTime;
    std::vector<unsigned int> matchedBarcodeReads;
    AdapterTrimmingStats adapterTrimmingStats;

    GeneralStats(): removedN(0), removedDemultiplex(0), removedQuality(0), uncalledBases(0), removedShort(0), readCount(0), processTime(0), ioTime(0) {};
    GeneralStats(unsigned int N, unsigned int numAdapters) : GeneralStats() 
    { 
        matchedBarcodeReads.resize(N); 
        adapterTrimmingStats.numRemoved.resize(numAdapters);
    };
    GeneralStats(const GeneralStats& rhs) = default;
    GeneralStats(GeneralStats&& rhs) = default;
    //{
    //    removedN = rhs.removedN;
    //    removedDemultiplex = rhs.removedDemultiplex;
    //    removedQuality = rhs.removedQuality;
    //    uncalledBases = rhs.uncalledBases;
    //    removedShort = rhs.removedShort;
    //    readCount = rhs.readCount;
    //    processTime = rhs.processTime;
    //    ioTime = rhs.ioTime;
    //    matchedBarcodeReads = rhs.matchedBarcodeReads;
    //    adapterTrimmingStats = rhs.adapterTrimmingStats;
    //};
    GeneralStats& operator=(const GeneralStats& rhs) = default;
    GeneralStats& operator=(GeneralStats&& rhs) = default;
    //{
    //    removedN = rhs.removedN;
    //    removedDemultiplex = rhs.removedDemultiplex;
    //    removedQuality = rhs.removedQuality;
    //    uncalledBases = rhs.uncalledBases;
    //    removedShort = rhs.removedShort;
    //    readCount = rhs.readCount;
    //    processTime = rhs.processTime;
    //    ioTime = rhs.ioTime;
    //    matchedBarcodeReads = rhs.matchedBarcodeReads;
    //    adapterTrimmingStats = rhs.adapterTrimmingStats;
    //    return *this;
    //}

    GeneralStats& operator+=(const GeneralStats& rhs)
    {
        removedN += rhs.removedN;
        removedDemultiplex += rhs.removedDemultiplex;
        removedQuality += rhs.removedQuality;
        uncalledBases += rhs.uncalledBases;
        removedShort += rhs.removedShort;
        readCount += rhs.readCount;
        processTime += rhs.processTime;
        ioTime += rhs.ioTime;
        if (matchedBarcodeReads.size() != rhs.matchedBarcodeReads.size())
            matchedBarcodeReads.resize(rhs.matchedBarcodeReads.size());
        matchedBarcodeReads = matchedBarcodeReads + rhs.matchedBarcodeReads;
        adapterTrimmingStats += rhs.adapterTrimmingStats;
        return *this;
    }
};

#endif
