// ==========================================================================
//                              OEA EXTRACTOR
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SHARED_ORPHAN_EXTRACTOR_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SHARED_ORPHAN_EXTRACTOR_H_

#include <string>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "progress_indicator.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class OrphanExtractor
// --------------------------------------------------------------------------

// Extract orphans from a indexed BAM file.

class OrphanExtractor
{
public:

    OrphanExtractor(char const * bamInFile, char const * fastqOutFile, int verbosity = 1) :
            bamInFile(bamInFile), fastqOutFile(fastqOutFile), verbosity(verbosity)
    {}

    // Returns the number of orphans.
    size_t run()
    {
        return createFastq();
    }

private:

    // Create the orphans.fq file.  Returns seqan::maxValue<size_t>() if exists already.
    size_t createFastq();

    // The BamStream to use for reading the reads and the BamIndex to use for jumping.
    seqan::BamStream bamStream;
    seqan::BamIndex<seqan::Bai> bamIndex;
    // The SequenceStream to use for writing orphans.
    seqan::SequenceStream orphansStream;

    // The options to use.
    std::string bamInFile;
    std::string fastqOutFile;
    int verbosity;
};

size_t OrphanExtractor::createFastq()
{
    size_t numOrphans = 0;

    // Open orphans file for writing.
    open(orphansStream, fastqOutFile.c_str(), seqan::SequenceStream::WRITE);
    if (!isGood(orphansStream))
        throw std::runtime_error("Could not open orphans FASTQ file for writing!");

    // Open BAM and BAI file.
    open(bamStream, bamInFile.c_str());
    if (!isGood(bamStream))
        throw std::runtime_error("Could not open input mapping for orphans extraction!");
    seqan::CharString baiFilename = bamInFile;
    append(baiFilename, ".bai");
    seqan::BamIndex<seqan::Bai> bamIndex;
    if (read(bamIndex, toCString(baiFilename)) != 0)
        throw std::runtime_error("Could not open BAI file for orphans extraction!");

    // Jump to orphans.
    if (!jumpToOrphans(bamStream, bamIndex))
        throw std::runtime_error("Could not jump to orphans in BAM file.");

    // Computation for progress display and progress bar for this.
    unsigned const MIB = 1024 * 1024;
    __int64 pos = positionInFile(bamStream) / MIB;
    __int64 size = fileSize(bamStream) / MIB;
    unsigned const MAX_BATCH_SIZE = 1000;
    unsigned batchSize = 0;
    ProgressBar pb(std::cerr, pos, size, verbosity == 1);
    pb.setLabel("  orphan extraction [MiB of BAM]");
    pb.updateDisplay();

    // Actually read the orphans.
    seqan::BamAlignmentRecord recordL, recordR;
    while (!atEnd(bamStream))
    {
        // Read pair of orphans.
        if (readRecord(recordL, bamStream) != 0)
            throw std::runtime_error("Could not read first mate record from BAM file.");
        if (!hasFlagUnmapped(recordL) || !hasFlagNextUnmapped(recordL))
            throw std::runtime_error("First mate record is not an orphans!");
        if (readRecord(recordR, bamStream) != 0)
            throw std::runtime_error("Could not read second mate record from BAM file.");
        if (!hasFlagUnmapped(recordR) || !hasFlagNextUnmapped(recordR))
            throw std::runtime_error("Second mate record is not an orphans!");

        // Update progress bar.
        if (++batchSize > MAX_BATCH_SIZE)
        {
            batchSize = 0;
            pos = positionInFile(bamStream);
            pb.advanceTo(pos / MIB);
        }

        // Perform sanity checks.
        if (recordL.qName != recordR.qName)
        {
            std::stringstream ss;
            ss << "\nERROR: Invalid pairing for reads " << recordL.qName << ", " << recordR.qName << ". "
               << "Are the records sorted by QNAME?";
            throw std::runtime_error(ss.str());
        }

        // Enforce ordering.
        if (hasFlagFirst(recordR))
        {
            using std::swap;
            swap(recordL, recordR);
        }

        // Write out left and right mate of orphan pair.
        if (writeRecord(orphansStream, recordL.qName, recordL.seq, recordL.qual) != 0)
            throw std::runtime_error("Could not write first mate to orphans FASTQ file.");
        if (writeRecord(orphansStream, recordR.qName, recordR.seq, recordR.qual) != 0)
            throw std::runtime_error("Could not write second mate to orphans FASTQ file.");

        numOrphans += 2;
    }

    // Finish progress bar display.
    pb.finish();

    // Manually close orphans stream to protect against crashes later that might have kep the stream unflushed.
    close(orphansStream);
    
    return numOrphans;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SHARED_ORPHAN_EXTRACTOR_H_
