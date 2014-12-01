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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_MATE_LINK_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_MATE_LINK_H_

#include <iosfwd>
#include <utility>

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class MateEdgeLabel
// ----------------------------------------------------------------------------

struct MateEdgeLabel
{
    // Length mean and standard deviation.
    double lengthMean, lengthStdDev;
    // Number of links represented by this edge.
    int count;
    // Weight, only used for sorting, not for statistical computations.
    double weight;

    explicit
    MateEdgeLabel(double lengthMean = 0, double lengthStdDev = 0, int count = 1, double weight = 1.0) :
            lengthMean(lengthMean), lengthStdDev(lengthStdDev), count(count), weight(weight)
    {}

    std::ostream & print(std::ostream & out) const;
};

std::ostream & operator<<(std::ostream & out, MateEdgeLabel const & label);

// ----------------------------------------------------------------------------
// Class MateLink
// ----------------------------------------------------------------------------

struct MateLink
{
    unsigned source, target;
    MateEdgeLabel label;

    explicit MateLink(unsigned source = -1, unsigned target = -1) :
            source(source), target(target), label()
    {}

    MateLink(unsigned source, unsigned target, MateEdgeLabel label) :
            source(source), target(target), label(label)
    {}

    std::ostream & print(std::ostream & out) const;

    bool operator<(MateLink const & other) const
    {
        return std::make_pair(source, target) < std::make_pair(other.source, other.target);
    }
};

std::ostream & operator<<(std::ostream & out, MateLink const & link);

// ----------------------------------------------------------------------------
// Class ContigEdgeLabel
// ----------------------------------------------------------------------------

struct ContigEdgeLabel
{
    unsigned length;

    explicit ContigEdgeLabel(unsigned length) : length(length)
    {}
};

std::ostream & operator<<(std::ostream & out, ContigEdgeLabel const & label);

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_MATE_LINK_H_
