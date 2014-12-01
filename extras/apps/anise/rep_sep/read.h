// ==========================================================================
//                                   ANISE
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_H_

#include "rep_sep/weighted_feature_vector.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Read
// ----------------------------------------------------------------------------

// Represents a read or super read.

struct Read
{
    static const unsigned INVALID = (unsigned)-1;

    Read() = default;

    explicit Read(unsigned id, unsigned contigID = INVALID, int beginPos = -1, int endPos = -1) :
            id(id), contigID(contigID), beginPos(beginPos), endPos(endPos)
    {}

    // Return number of conflicts (smaller counts).
    int conflicts(Read const & other) const;

    // Updating merging.  Up to ignoreConflicts are ignored upon merging.
    //
    // Note that the resulting weights in clique merging do not represent actual counts since we merge with redundancy there.
    Read & mergeWithThis(Read const & other, int ignoreConflicts = 0, bool ignoreOverlap = false, bool ignoreContigID = false);

    void addSubRead(unsigned readID);

    // Returns whether *this overlaps with other and the overlap rate is at least minRate * l where l is the smaller
    // read length of *this and other.
    bool overlaps(Read const & other, double minRate = 0) const
    {
        if (contigID != other.contigID)
            return false;
        if (!(beginPos < other.endPos && other.beginPos < endPos))
            return false;

        // Compute minimal overlap size.
        int minLen = minRate * std::min(other.length(), length()) + 0.0001;
        return (overlapLength(other) >= minLen);
    }

    // Returns length.
    int length() const { return endPos - beginPos; }

    // Returns length of overlap with other.
    int overlapLength(Read const & other) const
    {
        if (!(beginPos < other.endPos && other.beginPos < endPos))
            return 0;
        return std::min(endPos, other.endPos) - std::max(beginPos, other.beginPos);
    }

    // ID of this read.
    unsigned id { INVALID };
    // Contig ID.
    unsigned contigID { INVALID };
    // Begin and end position of the read with respect to the reference.
    int beginPos { -1 };
    int endPos   { -1 };
    // IDs of reads this read consists of.
    std::vector<unsigned> subReads;
    // The read features.
    WeightedFeatureVector features;
};

inline void swap(Read & lhs, Read & rhs)
{
    using std::swap;
    swap(lhs.id, rhs.id);
    swap(lhs.contigID, rhs.contigID);
    swap(lhs.beginPos, rhs.beginPos);
    swap(lhs.endPos, rhs.endPos);
    swap(lhs.subReads, rhs.subReads);
    swap(lhs.features, rhs.features);
}

std::ostream & operator<<(std::ostream & out, Read const & read);

// Functor for comparing two reads by their ID.

bool ltID(Read const & lhs, Read const & rhs);

// Functor for comparing two reads by their contig ID.

bool ltContigID(Read const & lhs, Read const & rhs);

// Functor for comparing two reads by their coordinate.

bool ltCoordinate(Read const & lhs, Read const & rhs);

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_READ_H_
