// ==========================================================================
//                                 BASIL
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

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_BAM_FILTER_PIPELINE_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_BAM_FILTER_PIPELINE_H_

#include <memory>
#include <stdexcept>
#include <functional>
#include <vector>

#include <seqan/bam_io.h>

#include "library_info.h"
#include "thread_safe_queue.h"

// ============================================================================
// Forwards
// ============================================================================

// Forward to the implementation of BamFilter.
class BamFilterImpl;
// Forward to the implementation of BamReader.
class BamReaderImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Typedefs
// ----------------------------------------------------------------------------

// A typedef for the BamIOContext specialization we use here.

typedef seqan::BamIOContext<seqan::StringSet<seqan::CharString> > TBamIOContext;

// ----------------------------------------------------------------------------
// Class BamReaderOptions
// ----------------------------------------------------------------------------

// Configuration for the BAM reading.

struct BamReaderOptions
{
    // The path to the file name to read from.
    seqan::CharString bamFileName;
    // Read that many chunks into one buffer.
    unsigned chunkSize;

    BamReaderOptions() : chunkSize(20*1000)
    {}
};

// ----------------------------------------------------------------------------
// Class BamFilterOptions
// ----------------------------------------------------------------------------

// Configuration for the BAM filter pipeline.

struct BamFilterOptions
{
    // Maximal coverage to allow.
    int maxCoverage;
    // Minimal alignment quality to accept.
    int minAlignmentQuality;
    // Smallest clipping size to consider.
    int clippingMinLength;
    // Number of threads to use for realigment.
    int realignmentNumThreads;
    // Chunk size for realignment.
    int realignmentChunkSize;
    // Maximal alignment length.
    int maxAlignmentLength;

    // Library size information.
    BamLibraryInfo libraryInfo;
    // Factor to get from std dev to allowed fragment size.
    double fragmentSizeFactor;

    // Reference filename.
    seqan::CharString refFileName;

    // Disable realignment.
    bool disableRealignment;

    // Computes the maximal fragment size from the configured options.
    int maxFragmentSize() const
    {
        return libraryInfo.median + fragmentSizeFactor * libraryInfo.stdDev;
    }

    BamFilterOptions() :
            maxCoverage(0), minAlignmentQuality(0), clippingMinLength(0), realignmentNumThreads(1),
            realignmentChunkSize(0), maxAlignmentLength(0), fragmentSizeFactor(0), disableRealignment(false)
    {}
};

// ----------------------------------------------------------------------------
// Class BamFilter
// ----------------------------------------------------------------------------

// Hold the state of the BAM filter pipeline.
//
// The function makePipline() returns a Pipeline<> with the filter.

class BamFilter
{
public:
    BamFilter(TBamIOContext & context, BamFilterOptions const & options);
    ~BamFilter();  // for pimpl

    // Creates the pipeline segment for the BAM filter.
    Pipeline<std::vector<seqan::BamAlignmentRecord *>> makePipeline(int bufferSize);

private:
    std::unique_ptr<BamFilterImpl> impl;
};

// ----------------------------------------------------------------------------
// Class BamReader
// ----------------------------------------------------------------------------

class BamReader
{
public:
    BamReader(BamReaderOptions const & options);
    ~BamReader();  // for pimpl

    // Whether or not the reader is at the end.
    bool atEnd() const;
    // Read records into the queue back.
    void read(std::vector<seqan::BamAlignmentRecord *> & out);

    // Return the BamIOContext of the _bamStream.
    TBamIOContext & bamIOContext();
    // Return the internal BamFileIn.
    seqan::BamFileIn const & bamFileIn();

    // Set a callback that is called regularly after having read some records.
    void setProgressCallback(std::function<void()> fun);

private:
    std::unique_ptr<BamReaderImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_BAM_FILTER_PIPELINE_H_
