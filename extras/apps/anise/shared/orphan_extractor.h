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

    // The BamFileIn to use for reading the reads and the BamIndex to use for jumping.
    seqan::BamFileIn bamFileIn;
    seqan::BamIndex<seqan::Bai> bamIndex;
    // The SequenceStream to use for writing orphans.
    seqan::SeqFileOut orphansFileOut;

    // The options to use.
    std::string bamInFile;
    std::string fastqOutFile;
    int verbosity;
};

size_t OrphanExtractor::createFastq()
{
    size_t numOrphans = 0;

    // Open orphans file for writing.
    if (!open(orphansFileOut, fastqOutFile.c_str()))
        throw std::runtime_error("Could not open orphans FASTQ file for writing!");

    // Open BAM and BAI file.
    if (!open(bamFileIn, bamInFile.c_str()))
        throw std::runtime_error("Could not open input mapping for orphans extraction!");
    seqan::CharString baiFilename = bamInFile;
    append(baiFilename, ".bai");
    seqan::BamIndex<seqan::Bai> bamIndex;
    if (!open(bamIndex, toCString(baiFilename)))
        throw std::runtime_error("Could not open BAI file for orphans extraction!");

    // Jump to orphans.
    bool hasAlignments = false;
    if (!jumpToOrphans(bamFileIn, hasAlignments, bamIndex))
        throw std::runtime_error("Could not jump to orphans in BAM file.");
    (void)hasAlignments;

    // Computation for progress display and progress bar for this.
    unsigned const MIB = 1024 * 1024;
    __int64 pos = position(bamFileIn) / MIB;
    __int64 size = fileSize(bamFileIn) / MIB;
    unsigned const MAX_BATCH_SIZE = 1000;
    unsigned batchSize = 0;
    ProgressBar pb(std::cerr, pos, size, verbosity == 1);
    pb.setLabel("  orphan extraction [MiB of BAM]");
    pb.updateDisplay();

    // Actually read the orphans.
    seqan::BamAlignmentRecord recordL, recordR;
    while (!atEnd(bamFileIn))
    {
        // Read pair of orphans.
        try
        {
            readRecord(recordL, bamFileIn);
            readRecord(recordR, bamFileIn);
        }
        catch (seqan::ParseError const & e)
        {
            throw std::runtime_error("Could not read second mate record from BAM file.");
        }
        if (!hasFlagUnmapped(recordL) || !hasFlagNextUnmapped(recordL))
            throw std::runtime_error("First mate record is not an orphans!");
        if (!hasFlagUnmapped(recordR) || !hasFlagNextUnmapped(recordR))
            throw std::runtime_error("Second mate record is not an orphans!");

        // Update progress bar.
        if (++batchSize > MAX_BATCH_SIZE)
        {
            batchSize = 0;
            pos = position(bamFileIn);
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
        try
        {
            writeRecord(orphansFileOut, recordL.qName, recordL.seq, recordL.qual);
            writeRecord(orphansFileOut, recordR.qName, recordR.seq, recordR.qual);
        }
        catch (seqan::IOError const & e)
        {
            throw std::runtime_error("Could not write second to orphans FASTQ file.");
        }

        numOrphans += 2;
    }

    // Finish progress bar display.
    pb.finish();

    // Manually close orphans stream to protect against crashes later that might have kep the stream unflushed.
    close(orphansFileOut);
    
    return numOrphans;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SHARED_ORPHAN_EXTRACTOR_H_
