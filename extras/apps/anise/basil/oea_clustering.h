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
// Iterative OEA clustering step.
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_BASIL_OEA_CLUSTERING_H_
#define SANDBOX_HERBARIUM_APPS_BASIL_OEA_CLUSTERING_H_

#include <memory>
#include <vector>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include "library_info.h"

// ============================================================================
// Forwards
// ============================================================================

// Forward to the implementation of OeaClusterAlgo.
class OeaClusterAlgoImpl;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Enum OeaClusterSelection
// ----------------------------------------------------------------------------

// Enum for cluster selection approach.

enum class OeaClusterSelection
{
    CHAINING,
    SET_COVER
};

// ----------------------------------------------------------------------------
// Class OeaClusterOptions
// ----------------------------------------------------------------------------

// Configuration for the OEA clustering step.

struct OeaClusterOptions
{
    // Library info as given from command line.  Automatically determined if not given.
    BamLibraryInfo libraryInfo;
    // Factor to get from std dev to allowed fragment size.
    double fragmentSizeFactor;
    // Whether or not to enable debug output.
    bool debug;

    // Configuration for cluster selection.
    OeaClusterSelection clusterSelection;
    // Minimal number of OEA anchors for a tentative site to be considered.
    int oeaMinSupport;
    // Minimal number of OEA anchors on each side for a tentative site to be considered.
    int oeaMinSupportEachSide;

    OeaClusterOptions() : fragmentSizeFactor(3), debug(false)
    {}

    int maxFragmentSize() const
    {
        return libraryInfo.median + fragmentSizeFactor * libraryInfo.stdDev;
    }
};

// ----------------------------------------------------------------------------
// Class OeaClusterRecord
// ----------------------------------------------------------------------------

// A record describing an OEA record.  A vector of these is the result of OeaClusterAlgo.

struct OeaClusterRecord
{
    explicit
    OeaClusterRecord(int rID = -1, int beginPos = -1, int endPos = -1, int centerPos = -1,
                     int leftWeight = 0, int rightWeight = 0, int totalWeight = 0) :
            centerPos(centerPos), leftWeight(leftWeight), rightWeight(rightWeight), totalWeight(totalWeight)
    {
        region.seqId = rID;
        region.beginPos = beginPos;
        region.endPos = endPos;
    }

    bool operator<(OeaClusterRecord const & other) const
    {
        return std::make_pair(std::make_pair(region.seqId, region.beginPos), region.endPos) <
            std::make_pair(std::make_pair(other.region.seqId, other.region.beginPos), other.region.endPos);
    }

    // Returns geometric mean of +1 pseudocounts of the left and right weight.
    double score() const
    {
        return sqrt((1.0 + leftWeight) * (1.0 + rightWeight));
    }

    // Region of the cluster.
    seqan::GenomicRegion region;
    // The center of the record.
    int centerPos;

    // Number of supporting records to the left/right/both.
    int leftWeight, rightWeight, totalWeight;
};

// ----------------------------------------------------------------------------
// Class OeaClusterAlgo
// ----------------------------------------------------------------------------

// Algorithm for clustering of OEA anchor records.
//
// This can be used to iteratively feed BamAlignmentRecords into the clustering algorithm.  Make sure to call finish()
// when the record stream is done.  The result can be retrieved using result().

class OeaClusterAlgo
{
public:
    OeaClusterAlgo(OeaClusterOptions const & options);
    ~OeaClusterAlgo();  // for pimpl

    // Register next BamAlignmentRecords.
    void push(std::vector<seqan::BamAlignmentRecord *> const & records);

    // Finish computation.
    void finish();

    // Return the generated OEA clusters.
    std::vector<OeaClusterRecord> result();

private:
    std::unique_ptr<OeaClusterAlgoImpl> impl;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_BASIL_OEA_CLUSTERING_H_
