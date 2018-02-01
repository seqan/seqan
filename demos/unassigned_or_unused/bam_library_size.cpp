// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Tobias Rausch <rausch@embl.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Shows how to use the bam_io module to estimate the insert size from a SAM
// or BAM file.
// ==========================================================================

#include <iostream>
#include <vector>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#if SEQAN_HAS_ZLIB

using namespace seqan;

// Stores information about a library.

struct LibraryInfo
{
    // This is the notation introduced by Illumina.
    //
    // F+: >R1 R2>
    // F-: >R2 R1>
    // R+: >R1 R2<
    // R-: <R2 R1>
    enum Orientation
    {
        F_PLUS  = 0,
        F_MINUS = 1,
        R_PLUS  = 2,
        R_MINUS = 3
    };

    unsigned median;
    double   stdDev;
    unsigned maxNormalISize;
    Orientation defaultOrient;

    LibraryInfo() :
        median(0), stdDev(0), maxNormalISize(0), defaultOrient(F_PLUS)
    {}
};

inline int
getStrandIndependentOrientation(BamAlignmentRecord const & rec)
{
    if (!hasFlagRC(rec))
    {
        if (!(hasFlagNextRC(rec)))
        {
            return (rec.beginPos < rec.pNext) ? LibraryInfo::F_PLUS : LibraryInfo::F_MINUS;
        }
        else
        {
            return (rec.beginPos < rec.pNext) ? LibraryInfo::R_PLUS : LibraryInfo::R_MINUS;
        }
    }
    else
    {
        if (!hasFlagNextRC(rec))
        {
            return (rec.beginPos > rec.pNext) ? LibraryInfo::R_PLUS : LibraryInfo::R_MINUS;
        }
        else
        {
            return (rec.beginPos > rec.pNext) ? LibraryInfo::F_PLUS : LibraryInfo::F_MINUS;
        }
    }
}

bool endsWith(CharString const & str, CharString const & x)
{
    typedef Size<CharString>::Type TSize;

    TSize len = std::min(length(str), length(x));
    TSize pos = length(str) - len;
    return suffix(str, pos) == x;
}

bool performEstimation(LibraryInfo & libInfo, BamFileIn & bamFileIn)
{
    // Vector of all insert sizes.
    typedef std::vector<unsigned int> TVecISize;
    TVecISize vecISize;

    // Read Header.
    BamHeader header;
    readHeader(header, bamFileIn);

    // Orientations.
    unsigned orientationCounters[4] = {0, 0, 0, 0};

    // Running mean.
    uint64_t runningMean = 0;
    unsigned pairCount = 0;

    // Count.
    BamAlignmentRecord record;
    while (!atEnd(bamFileIn))
    {
        readRecord(record, bamFileIn);

        // Skip all records except for first mate of properly mapped pairs.
        if (!hasFlagMultiple(record) || hasFlagUnmapped(record) || hasFlagNextUnmapped(record) ||
            hasFlagSecondary(record) || hasFlagFirst(record))
            continue;

        vecISize.push_back(abs(record.tLen));
        runningMean += abs(record.tLen);
        pairCount += 1;

        orientationCounters[getStrandIndependentOrientation(record)] += 1;
    }

    // Get default orientation.
    unsigned orientMax = 0;
    for (unsigned i = 0; i < 4; ++i)
        if (orientationCounters[i] > orientationCounters[orientMax])
            orientMax = i;
    LibraryInfo::Orientation orientations[4] = { LibraryInfo::F_PLUS, LibraryInfo::F_MINUS, LibraryInfo::R_PLUS, LibraryInfo::R_MINUS };
    libInfo.defaultOrient = orientations[orientMax];

    // Trim off the chimera peak in mate-pair libraries.
    if ((libInfo.defaultOrient == LibraryInfo::R_MINUS) && (pairCount > 0) && (runningMean / pairCount) >= 1000u)
    {
        typedef TVecISize::const_iterator TVecISizeIter;
        TVecISize vecISizeTmp;
        for (TVecISizeIter it = vecISize.begin(); it < vecISize.end(); ++it)
            if (*it > 1000)
                vecISizeTmp.push_back(*it);
        std::swap(vecISize, vecISizeTmp);
    }


    // Check that this is a proper paired-end library
    if (vecISize.empty())
        return true;

    // Get library stats.
    //
    // Start with median.
    typedef TVecISize::iterator TVecISizeIter;
    TVecISizeIter begin = vecISize.begin();
    TVecISizeIter end = vecISize.end();
    std::nth_element(begin, begin + (end - begin) / 2, end);
    begin = vecISize.begin();
    end = vecISize.end();
    libInfo.median = *(begin + (end - begin) / 2);
    // Standard deviation is next.
    //
    // SD calculation cutoffs are 7 SDs to the left and right assuming 10% library deviation.
    libInfo.stdDev = 0;
    double cutoffMax = libInfo.median + 7 * 0.1 * libInfo.median;
    double cutoffMin = libInfo.median - 7 * 0.1 * libInfo.median;
    if ((cutoffMin < 0) || (cutoffMax < cutoffMin))
        cutoffMin = 0;
    unsigned int count = 0;
    for (; begin < end; ++begin)
    {
        if ((*begin >= cutoffMin) && (*begin <= cutoffMax))
        {
            libInfo.stdDev += (*begin - libInfo.median) * (*begin - libInfo.median);
            ++count;
        }
    }
    if (count == 0u)  // prevent div-by-zero below
        count = 1;
    libInfo.stdDev = sqrt(libInfo.stdDev / count);
    libInfo.maxNormalISize = static_cast<unsigned>(libInfo.median + 3 * libInfo.stdDev);

    return true;
}

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "Invalid arguments!\n"
                  << "USAGE: bam_library_size {IN.sam,IN.bam}\n";
        return 1;
    }

    LibraryInfo libInfo;

    BamFileIn bamFileIn;
    if (!open(bamFileIn, argv[1]))
    {
        std::cerr << "Could not open input SAM/BAM file " << argv[1] << "\n";
        return 1;
    }

    if (!performEstimation(libInfo, bamFileIn))
        return 1;

    // Print result.
    std::cout << "Library Information\n\n"
              << "path:                       " << argv[1] << "\n"
              << "median:                     " << libInfo.median << "\n"
              << "standard deviation:         " << libInfo.stdDev << "\n"
              << "maximum normal insert size: " << libInfo.maxNormalISize << "\n";

    std::cout << "orientation:                ";
    switch (libInfo.defaultOrient)
    {
    case LibraryInfo::F_PLUS:
        std::cout << "F+ R1 ---> ---> R2\n";
        break;

    case LibraryInfo::F_MINUS:
        std::cout << "F- R1 ---> ---> R2\n";
        break;

    case LibraryInfo::R_PLUS:
        std::cout << "R+ R1 ---> <--- R2\n";
        break;

    case LibraryInfo::R_MINUS:
        std::cout << "R- R1 <--- ---> R2\n";
        break;
    }

    return 0;
}

#else

int main(int, char const **)
{
    std::cerr << "bam_library_size can only be compiled correctly with zlib." << std::endl;
    return 0;
}

#endif  // #if SEQAN_HAS_ZLIB
