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

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_CLIPPING_CLUSTERING_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_CLIPPING_CLUSTERING_H_

#include <memory>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

// ============================================================================
// Forwards
// ============================================================================

// Forward to the implementation of ClippingClusterAlgo.
class ClippingClusterAlgoImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ClippingClusterOptions
// ----------------------------------------------------------------------------

// Configuration for the clipping clustering.

struct ClippingClusterOptions
{
    // Maximal assumed alignment length.
    int maxAlignmentLength;
    // The window size for counting clippings.
    int minCoverageWindowLength;
    // Smallest number of clipping positions to look for.
    int minCoverage;
    // Clippings below this length are ignored in the clustring.
    int minClippingLength;

    ClippingClusterOptions() :
            maxAlignmentLength(0), minCoverageWindowLength(0), minCoverage(0), minClippingLength(0)
    {}
};

// ----------------------------------------------------------------------------
// Class ClippingClusterRecord
// ----------------------------------------------------------------------------

// These records make up the clipping clusters.

struct ClippingClusterRecord
{
    // Region of the cluster.
    seqan::GenomicRegion region;
    // The center of the record.
    int centerPos;

    // Number of supporting records to the left/right/both.
    int leftWeight, rightWeight, totalWeight;

    explicit
    ClippingClusterRecord(int rID = -1, int beginPos = -1, int endPos = -1, int centerPos = -1,
                          int leftWeight = 0, int rightWeight = 0, int totalWeight = 0) :
            centerPos(centerPos), leftWeight(leftWeight), rightWeight(rightWeight), totalWeight(totalWeight)
    {
        region.seqId = rID;
        region.beginPos = beginPos;
        region.endPos = endPos;
    }

    bool operator<(ClippingClusterRecord const & other) const
    {
        if (region.seqId == -1 && other.region.seqId == -1)
            return false;  // equal
        if (other.region.seqId == -1)
            return true;  // smaller than dummy
        return std::make_pair(std::make_pair(region.seqId, region.beginPos), region.endPos) <
            std::make_pair(std::make_pair(other.region.seqId, other.region.beginPos), other.region.endPos);
    }

    // Returns geometric mean of +1 pseudocounts of the left and right weight.
    double score() const
    {
        return sqrt((1.0 + leftWeight) * (1.0 + rightWeight));
    }
};

// ----------------------------------------------------------------------------
// Class ClippingClusterAlgo
// ----------------------------------------------------------------------------

// Class for the iterative computation of clipping position clusters.

class ClippingClusterAlgo
{
public:
    ClippingClusterAlgo(ClippingClusterOptions const & options);
    ~ClippingClusterAlgo();  // for pimpl

    // Register next BamAlignmentRecords.
    void push(std::vector<seqan::BamAlignmentRecord *> const & records);

    // Finish computation.
    void finish();

    // Return the generated clipping clusters.
    std::vector<ClippingClusterRecord> result();

private:
    std::unique_ptr<ClippingClusterAlgoImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_CLIPPING_CLUSTERING_H_
