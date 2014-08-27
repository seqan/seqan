// ==========================================================================
//                                   ANISE
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_SET_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_SET_H_

#include <iosfwd>
#include <vector>

#include "rep_sep/read.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FeatureReadSet
// ----------------------------------------------------------------------------

// A set of reads, input for clique enumeration.

struct FeatureReadSet
{
    size_t size() const { return reads.size(); }
    bool empty() const { return reads.empty(); }

    unsigned numActualReads { 0 };  // total number of real reads
    unsigned numFeatures { 0 };  // number of features
    std::vector<Read> reads;  // sorted by (contigID, beginPos).

    std::vector<Read>::iterator begin() { return reads.begin(); }
    std::vector<Read>::iterator end() { return reads.end(); }
    std::vector<Read>::const_iterator begin() const { return reads.begin(); }
    std::vector<Read>::const_iterator end() const { return reads.end(); }

    void print(std::ostream & out) const;
};

inline void swap(FeatureReadSet & lhs, FeatureReadSet & rhs)
{
    using std::swap;
    swap(lhs.numActualReads, rhs.numActualReads);
    swap(lhs.numFeatures, rhs.numFeatures);
    swap(lhs.reads, rhs.reads);
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function removeReads()
// ----------------------------------------------------------------------------

// Remove reads where removeContig[contigID] is true and renumber the contigID members afterwards.

void removeReads(FeatureReadSet & out,
                 FeatureReadSet const & in,
                 std::vector<bool> removeContig);

// ----------------------------------------------------------------------------
// Function fixFeatureCounts()
// ----------------------------------------------------------------------------

// Fixes feature counts.

void fixFeatureCounts(FeatureReadSet & readSet, FeatureReadSet const & atomicReadSet);

// ----------------------------------------------------------------------------
// Function distributeUnplacedReads()
// ----------------------------------------------------------------------------

// Distribute yet unplaced reads from atomic read set to read set.

void distributeUnplacedReads(FeatureReadSet & readSet,
                             FeatureReadSet const & atomicReadSet);

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_SET_H_
