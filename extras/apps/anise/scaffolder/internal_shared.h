// ==========================================================================
//                                  ANISE
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_INTERNAL_SHARED_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_INTERNAL_SHARED_H_

#include <iostream>

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class ReductionGraphLabel
// --------------------------------------------------------------------------

// Bundle the information required for transitive reduction.

struct ReductionGraphLabel
{
    // Kind of the edge type.
    enum EdgeType
    {
        CONTIG,
        MATE,
        INFERRED  // from greedy path merging algo
    };
    EdgeType edgeType;

    // Whether or not this is a contig edge.
    bool isContig() const
    {
        return (edgeType == CONTIG);
    }
    // Length mean and standard deviation.
    double lengthMean;
    double lengthStdDev;
    // The index in the input links if mate edge, contig id in case of contig, maxValue<unsigned>() in case of INFERRED.
    unsigned idx;
    // Count of the represented mates.
    int count;
    // Weight of the link, only for mate edges.
    double weight;

    explicit
    ReductionGraphLabel(EdgeType edgeType = MATE, double lengthMean = 0, double lengthStdDev = 0,
                        unsigned idx = -1, int count = -1, double _weight = -1) :
            edgeType(edgeType), lengthMean(lengthMean), lengthStdDev(lengthStdDev), idx(idx),
            count(count), weight(_weight)
    {}

    void print(std::ostream & out) const
    {
        out << "ReductionGraphLabel(edgeType=" << edgeType << ", lengthMean=" << lengthMean
            << ", lengthStdDev=" << lengthStdDev << ", idx="
            << idx << ", count=" << count << ", weight=" << weight <<")";
    }
};

inline std::ostream & operator<<(std::ostream & out, ReductionGraphLabel const & label)
{
    label.print(out);
    return out;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_INTERNAL_SHARED_H_
