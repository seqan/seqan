// ==========================================================================
//                                 BASIL
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_LIBRARY_INFO_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_LIBRARY_INFO_H_

#include <functional>

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BamLibraryInfo
// ----------------------------------------------------------------------------

// Represents the orientation and insert size of paired reads.

struct BamLibraryInfo
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

    double median;
    double stdDev;
    unsigned maxNormalISize;
    Orientation defaultOrient;

    // Average read length.
    int avgReadLen;

    BamLibraryInfo() : median(0), stdDev(0), maxNormalISize(0), defaultOrient(F_PLUS), avgReadLen(0)
    {}

    static char const * orientationStr(Orientation o);
};

// ----------------------------------------------------------------------------
// Class BamLibraryEstimator
// ----------------------------------------------------------------------------

// Allows to detect LibraryInfo information from BAM file.

class BamLibraryEstimator
{
public:
    // Number of records to evaluate.
    int numRecords;

    BamLibraryEstimator(int numRecords = 100*1000) : numRecords(numRecords)
    {}

    // Run estimation with and without progress indication update functor.
    int run(BamLibraryInfo & result, char const * path);
    int run(BamLibraryInfo & result, char const * path, std::function<void(int)> fn);
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getOrientationStr()
// ----------------------------------------------------------------------------

char const * getOrientationStr(BamLibraryInfo::Orientation o);

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_LIBRARY_INFO_H_
