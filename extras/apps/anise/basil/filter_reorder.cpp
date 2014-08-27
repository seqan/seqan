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

#include "filter_reorder.h"

#include <algorithm>
#include <deque>

namespace  // anonymous
{

// ----------------------------------------------------------------------------
// Function ltReorderRecord)
// ----------------------------------------------------------------------------

bool ltReorderRecord(seqan::BamAlignmentRecord const * lhs,
                     seqan::BamAlignmentRecord const * rhs)
{
    return std::make_pair(lhs->rID, lhs->beginPos) < std::make_pair(rhs->rID, rhs->beginPos);
}

}  // anonymous namespace

// ----------------------------------------------------------------------------
// Class ReorderFilterImpl
// ----------------------------------------------------------------------------

class ReorderFilterImpl
{
public:
    ReorderFilterImpl(int maxBeginPosDistance) : maxBeginPosDistance(maxBeginPosDistance)
    {}

    void filter(std::vector<seqan::BamAlignmentRecord *> & out,
                std::vector<seqan::BamAlignmentRecord *> const & in);

    void finish(std::vector<seqan::BamAlignmentRecord *> & out);

private:
    typedef std::deque<seqan::BamAlignmentRecord *> TBuffer;

    TBuffer buffer;

    // Maximal distance between the begin positions of a valid pair.  This can be the sum of the maximal fragment size
    // and the maximal alignment length.
    int maxBeginPosDistance;

    bool spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                               seqan::BamAlignmentRecord const & rhs,
                               int maxBeginPosDistance) const;
};

void ReorderFilterImpl::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                               std::vector<seqan::BamAlignmentRecord *> const & in)
{
    // Sort input order, we keep buffer sorted.
    std::vector<seqan::BamAlignmentRecord *> inCopy(in);
    std::sort(inCopy.begin(), inCopy.end(), ltReorderRecord);
    // Merge input buffer into buffer.
    std::vector<seqan::BamAlignmentRecord *> bufferCopy(buffer.begin(), buffer.end());
    buffer.resize(inCopy.size() + bufferCopy.size());
    std::merge(bufferCopy.begin(), bufferCopy.end(), inCopy.begin(), inCopy.end(),
               buffer.begin(), ltReorderRecord);

    // Write out from buffer as far as possible.
    while (!buffer.empty() && spanAboveFragmentSize(*buffer.front(), *buffer.back(), maxBeginPosDistance))
    {
        // auto ptr = buffer.front();
        // std::cerr << "REORDER WRITING\t" << ptr->qName << "/" << (hasFlagLast(*ptr) + 1) << "\t" << ptr->beginPos << "\n";
        out.push_back(buffer.front());
        buffer.pop_front();
    }
}

void ReorderFilterImpl::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    // for (auto ptr : buffer)
    //     std::cerr << "REORDER WRITING\t" << ptr->qName << "/" << (hasFlagLast(*ptr) + 1) << "\t" << ptr->beginPos << "\n";
    std::copy(buffer.begin(), buffer.end(), std::back_inserter(out));
    buffer.clear();
}

bool ReorderFilterImpl::spanAboveFragmentSize(seqan::BamAlignmentRecord const & lhs,
                                              seqan::BamAlignmentRecord const & rhs,
                                              int maxBeginPosDistance) const
{
    if (lhs.rID != rhs.rID)
        return true;
    return abs(lhs.beginPos - rhs.beginPos) > maxBeginPosDistance;
}

// ----------------------------------------------------------------------------
// Class ReorderFilter
// ----------------------------------------------------------------------------

ReorderFilter::ReorderFilter(int maxCoverage) : impl(new ReorderFilterImpl(maxCoverage))
{}

ReorderFilter::~ReorderFilter()
{}

void ReorderFilter::filter(std::vector<seqan::BamAlignmentRecord *> & out,
                           std::vector<seqan::BamAlignmentRecord *> const & in)
{
    impl->filter(out, in);
}

void ReorderFilter::finish(std::vector<seqan::BamAlignmentRecord *> & out)
{
    impl->finish(out);
}
