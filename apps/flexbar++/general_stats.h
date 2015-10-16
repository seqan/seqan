// ==========================================================================
//                              generalProcessing.h
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Benjamin Menkuec <benjamin@menkuec.de>
// ==========================================================================

#ifndef GENERALSTATS_H
#define GENERALSTATS_H


struct AdapterTrimmingStats
{
    std::vector<std::vector<unsigned>> removedLength;
    std::vector<unsigned> numRemoved;
    unsigned a1count, a2count;
    unsigned overlapSum;
    unsigned minOverlap, maxOverlap;

    AdapterTrimmingStats() : a1count(0), a2count(0), overlapSum(0),
        minOverlap(std::numeric_limits<unsigned>::max()), maxOverlap(0) {};

    AdapterTrimmingStats& operator+= (AdapterTrimmingStats const& rhs)
    {
        a1count += rhs.a1count;
        a2count += rhs.a2count;
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
        a1count = 0;
        a2count = 0;
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
