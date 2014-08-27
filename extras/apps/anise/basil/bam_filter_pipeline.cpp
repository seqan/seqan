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

#include "bam_filter_pipeline.h"

#include "filter_high_coverage.h"
#include "filter_discordant.h"
#include "filter_clean_pairs.h"
#include "filter_remove_clipping.h"
#include "filter_realignment.h"
#include "filter_reorder.h"
#include "filter_shared_types.h"

// ----------------------------------------------------------------------------
// Class BamReaderImpl
// ----------------------------------------------------------------------------

class BamReaderImpl
{
public:
    BamReaderImpl(BamReaderOptions const & options) :
            done(false), options(options), _bamStream(toCString(options.bamFileName)), recordsRead(0)
    {
        if (!isGood(_bamStream))
            throw BamFilterException("Could not open BAM file!");
    }

    bool atEnd() const;

    void read(std::vector<seqan::BamAlignmentRecord *> & out);

    TBamIOContext & bamIOContext();
    seqan::BamStream const & bamStream();
    void setProgressCallback(std::function<void()> fun);

private:
    bool done;  // true if read first orphan
    BamReaderOptions options;
    seqan::BamStream _bamStream;
    std::function<void()> progressCallback;

    unsigned recordsRead;
};

bool BamReaderImpl::atEnd() const
{
    return (done || seqan::atEnd(_bamStream));
}

TBamIOContext & BamReaderImpl::bamIOContext()
{
    return _bamStream.bamIOContext;
}

seqan::BamStream const & BamReaderImpl::bamStream()
{
    return _bamStream;
}

void BamReaderImpl::setProgressCallback(std::function<void()> fun)
{
    progressCallback = fun;
}

void BamReaderImpl::read(std::vector<seqan::BamAlignmentRecord *> & out)
{
    out.clear();

    for (unsigned i = 0; i < options.chunkSize && !atEnd(); ++i, ++recordsRead)
    {
        std::unique_ptr<seqan::BamAlignmentRecord> record(new seqan::BamAlignmentRecord);
        if (readRecord(*record, _bamStream) != 0)
            throw BamFilterException("Problem reading BAM record.");

        // static int j = 0;
        // printf("%d creations\n", ++j);

        if (record->rID == seqan::BamAlignmentRecord::INVALID_REFID)
        {
            done = true;
            break;  // Break at unaligned record.
        }

        if (hasFlagSecondary(*record) || hasFlagQCNoPass(*record) || hasFlagDuplicate(*record))
            continue;  // Ignore secondary and flagged records.

        if (recordsRead % (10 * 1000) == 0u)
            progressCallback();

        // fprintf(stderr, "ALLOCATED\t%p\n", record.get());
        out.push_back(record.release());
    }
}

// ----------------------------------------------------------------------------
// Class BamReader
// ----------------------------------------------------------------------------

BamReader::BamReader(BamReaderOptions const & options) : impl(new BamReaderImpl(options))
{}

BamReader::~BamReader()
{}

bool BamReader::atEnd() const
{
    return impl->atEnd();
}

void BamReader::read(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->read(out);
}

TBamIOContext & BamReader::bamIOContext()
{
    return impl->bamIOContext();
}

seqan::BamStream const & BamReader::bamStream()
{
    return impl->bamStream();
}

void BamReader::setProgressCallback(std::function<void()> fun)
{
    impl->setProgressCallback(fun);
}

// ----------------------------------------------------------------------------
// Class BamFilterImpl
// ----------------------------------------------------------------------------

class BamFilterImpl
{
public:
    typedef std::vector<seqan::BamAlignmentRecord *> TBuffer;

    BamFilterImpl(TBamIOContext & bamIOContext, BamFilterOptions const & options) :
            options(options),
            highCoverageFilter(options.maxCoverage),
            discordantPairsFilter(options.maxFragmentSize(), options.minAlignmentQuality),
            cleanPairsFilter(options.maxFragmentSize() + options.maxAlignmentLength),
            removeClippingFilter(options.maxFragmentSize() + options.maxAlignmentLength,
                                 options.clippingMinLength)
    {
        if (!options.disableRealignment)
        {
            realignmentFilter.reset(new RealignmentFilter(options.maxFragmentSize(), options.libraryInfo.defaultOrient,
                                                          options.refFileName, bamIOContext, options.realignmentNumThreads,
                                                          options.realignmentChunkSize));
            reorderFilter.reset(new ReorderFilter(options.maxFragmentSize() + options.maxAlignmentLength));
        }
    }

    Pipeline<TBuffer> makePipeline(int bufferSize);

private:

    template <typename TFilter>
    std::function<TBuffer(TBuffer &, bool)> makeOneSegment(TFilter & filter, std::string name);

    BamFilterOptions options;

    // The filter objects.
    HighCoverageFilter highCoverageFilter;
    DiscordantPairsFilter discordantPairsFilter;
    CleanPairsFilter cleanPairsFilter;
    RemoveClippingFilter removeClippingFilter;
    std::unique_ptr<RealignmentFilter> realignmentFilter;
    std::unique_ptr<ReorderFilter> reorderFilter;
};

Pipeline<std::vector<seqan::BamAlignmentRecord *>> BamFilterImpl::makePipeline(int bufferSize)
{
    Pipeline<std::vector<seqan::BamAlignmentRecord *>> result;
    result.append(makeOneSegment(highCoverageFilter, "high-coverage"), bufferSize)
            .append(makeOneSegment(discordantPairsFilter, "discordant-pairs"), bufferSize)
            .append(makeOneSegment(cleanPairsFilter, "clean-pairs"), bufferSize)
            .append(makeOneSegment(removeClippingFilter, "remove-clipping"), bufferSize);
    if (!options.disableRealignment)
    {
        result.append(makeOneSegment(*realignmentFilter, "realignment"), bufferSize);
        result.append(makeOneSegment(*reorderFilter, "reorder"), bufferSize);
    }
    return result;
}

template <typename TFilter>
std::function<typename BamFilterImpl::TBuffer(typename BamFilterImpl::TBuffer &, bool)>
BamFilterImpl::makeOneSegment(TFilter & filterObj, std::string name)
{
    // Driver code for the HighCoverageFilter.
    auto passThroughFilter =
            [&,name](TBuffer & inBuffer, bool isFinish)
            {
                TBuffer result;
                if (!isFinish)
                    filterObj.filter(result, inBuffer);
                else
                    filterObj.finish(result);
                return result;
            };
    return passThroughFilter;
}

// ----------------------------------------------------------------------------
// Class BamFilter
// ----------------------------------------------------------------------------

BamFilter::BamFilter(TBamIOContext & context, BamFilterOptions const & options) : impl(new BamFilterImpl(context, options))
{}

BamFilter::~BamFilter()
{}

Pipeline<std::vector<seqan::BamAlignmentRecord *>> BamFilter::makePipeline(int bufferSize)
{
    return impl->makePipeline(bufferSize);
}
