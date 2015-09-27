#ifndef GENERALSTATS_H_
#define GENERALSTATS_H_


struct AdapterTrimmingStats
{
    AdapterTrimmingStats() : a1count(0), a2count(0), overlapSum(0),
        minOverlap(std::numeric_limits<unsigned>::max()), maxOverlap(0) {};

    AdapterTrimmingStats& operator+= (AdapterTrimmingStats const& rhs)
    {
        a1count += rhs.a1count;
        a2count += rhs.a2count;
        overlapSum += rhs.overlapSum;
        minOverlap = minOverlap < rhs.minOverlap ? minOverlap : rhs.minOverlap;
        maxOverlap = maxOverlap < rhs.maxOverlap ? rhs.maxOverlap : maxOverlap;
        return *this;
    }
    void clear()
    {
        a1count = 0;
        a2count = 0;
        overlapSum = 0;
        minOverlap = std::numeric_limits<unsigned>::max();
        maxOverlap = 0;
    }

    unsigned a1count, a2count;
    unsigned overlapSum;
    unsigned minOverlap, maxOverlap;
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
    std::vector<unsigned> matchedBarcodeReads;
    AdapterTrimmingStats adapterTrimmingStats;

    GeneralStats() : removedN(0), removedDemultiplex(0), uncalledBases(0), removedShort(0), readCount(0), processTime(0), ioTime(0) {};
};

#endif
