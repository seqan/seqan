// ==========================================================================
//                                   ANISE
// ==========================================================================
// Copyright (c) 2011-2012, Manuel Holtgrewe, FU Berlin
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

#include "anise_options.h"

#include <seqan/file.h>

namespace   // anonymous namespace
{

// --------------------------------------------------------------------------
// Function getYesNo()
// --------------------------------------------------------------------------

char const * getYesNo(bool v)
{
    return v ? "yes" : "no";
}

// --------------------------------------------------------------------------
// Function getVerbosityStr()
// --------------------------------------------------------------------------

char const * getVerbosityStr(int v)
{
    switch (v)
    {
        case 0:
            return "quiet";
        case 1:
            return "normal";
        case 2:
            return "verbose";
        default:
            return "very-verbose";
    }
}

// --------------------------------------------------------------------------
// Function getTammiMethodStr()
// --------------------------------------------------------------------------

char const * getTammiMethodStr(rep_sep::ReadSeparatorOptions::Method m)
{
    if (m == rep_sep::ReadSeparatorOptions::TAMMI_SIMPLE)
        return "SIMPLE";
    else
        return "PHRED";
}

// --------------------------------------------------------------------------
// Function printOptions()
// --------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & stream, AniseOptions const & options)
{
    stream << "VERBOSITY                         \t" << getVerbosityStr(options.verbosity) << "\n"
           << "NUMBER OF THREADS                 \t" << options.numThreads << "\n"
           << "\n"
           << "INPUT REFERENCE FILE              \t" << options.inputReference << "\n"
           << "INPUT VCF FILE                    \t" << options.inputVcf << "\n"
           << "INPUT MAPPING FILE                \t" << options.inputMapping << "\n"
           << "OUTPUT FASTA FILE                 \t" << options.outputFasta << "\n"
           << "OUTPUT MAPPING FILE               \t" << options.outputMapping << "\n"
           << "OUTPUT DEBUG DIR                  \t" << options.outputDebugDir << "\n"
           << "TEMPORARY FILES DIR               \t" << options.tmpDir << "\n"
           << "\n"
           << "RECURSION MAX STEPS               \t" << options.recursionMaxSteps << "\n"
           << "REALIGN ASSEMBLY                  \t" << getYesNo(options.recursionMaxSteps) << "\n"
           << "\n"
           << "AUTO LIBRARY INFO                 \t" << getYesNo(options.autoLibraryInfo) << "\n"
           << "AUTO LIBRARY NUM RECORDS          \t" << options.autoLibraryNumRecords << "\n"
           << "CMD LINE LIBRARY INFO\n"
           << "  FRAGMENT LEN MEDIAN             \t" << options.libraryInfo.median << "\n"
           << "  FRAGMENT LEN STDDEV             \t" << options.libraryInfo.stdDev << "\n"
           // << "  FRAGMENT LEN MAX NORMAL         \t" << options.libraryInfo.maxNormalISize << "\n"
           << "  DEFAULT ORIENTATION             \t" << getOrientationStr(options.libraryInfo.defaultOrient) << "\n"
           << "FRAGMENT SIZE FACTOR              \t" << options.fragmentSizeFactor << "\n"
           << "MAX FRAGMENT SIZE                 \t" << options.maxFragmentSize() << "\n"
           << "\n"\
           << "ASSEMBLY SITE WINDOW RADIUS       \t" << options.assemblySiteWindowRadius << "\n"
           << "ASSEMBLY SITE FRINGE RADIUS       \t" << options.assemblySiteFringeRadius << "\n"
           << "\n"
           << "STOP INITIAL READ COUNT           \t" << options.stopInitialReadCount << "\n"
           << "STOP READ COUNT                   \t" << options.stopReadCount << "\n"
           << "STOP COVERAGE                     \t" << options.stopCoverage << "\n"
           << "\n"
           << "REALIGNMENT BANDWIDTH             \t" << options.realignmentBandwidth << "\n"
           << "REALIGNMENT BORDER                \t" << options.realignmentBorder << "\n"
           << "\n"
           << "REPEAT SEPARATION OPTIONS\n"
           << "  SEPARATE REPEATS                \t" << getYesNo(options.separateRepeats) << "\n"
           << "  TAMMI METHOD TO USE             \t" << getTammiMethodStr(options.readSepOptions.tammiMethod) << "\n"
           << "  PER BASE ERROR                  \t" << options.readSepOptions.pErr << "\n"
           << "  MAX RANDOM CORRELATION          \t" << options.readSepOptions.maxRandomCorrelation << "\n"
           << "  TAU MIN                         \t" << options.readSepOptions.tauMin << "\n"
           << "  R MIN                           \t" << options.readSepOptions.rMin << "\n"
           << "  START COMPRESSION AT            \t" << options.readSepOptions.startCompressionAt << "\n"
           << "  SPLIT D MIN                     \t" << getYesNo(options.readSepOptions.splitDMin) << "\n"
           << "\n"
           << "READ MAPPING ERROR RATE           \t" << options.readMappingErrorRate << "\n"
           << "READ MAPPING BATCH SIZE           \t" << options.readMappingBatchSize << "\n"
           << "\n"
           << "OVERLAPPER MIN OVERLAP RATIO      \t" << options.overlapperMinOverlapRatio << "\n"
           << "OVERLAPPER MAX ERROR RATE         \t" << options.overlapperMaxErrorRate << "\n"
           << "\n"
           << "MSA SCORE MATCH                   \t" << options.msaScoreMatch << "\n"
           << "MSA SCORE MISMATCH                \t" << options.msaScoreMismatch << "\n"
           << "MSA SCORE GAP OPEN                \t" << options.msaScoreGapOpen << "\n"
           << "MSA SCORE GAP EXTEND              \t" << options.msaScoreGapExtend << "\n"
           << "\n"
           << "CONSENSUS MIN BASE SUPPORT        \t" << options.consensusMinBaseSupport << "\n"
           << "CONSENSUS MIN CONTIG LENGTH RATE  \t" << options.consensusMinContigLengthRate << "\n";
    return stream;
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class AniseOptions
// ----------------------------------------------------------------------------

void AniseOptions::print(std::ostream & out) const
{
    out << *this;
}
