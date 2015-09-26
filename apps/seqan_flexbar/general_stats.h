#ifndef GENERALSTATS_H_
#define GENERALSTATS_H_

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

    GeneralStats() : removedN(0), removedDemultiplex(0), uncalledBases(0), removedShort(0), readCount(0), processTime(0), ioTime(0) {};
};

#endif
