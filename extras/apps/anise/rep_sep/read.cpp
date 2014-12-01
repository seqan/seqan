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

#include "read.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <tuple>

namespace rep_sep {

// --------------------------------------------------------------------------
// Class Read
// --------------------------------------------------------------------------

int Read::conflicts(Read const & other) const
{
    return features.conflicts(other.features);
}

void Read::addSubRead(unsigned readID)
{
    subReads.insert(std::lower_bound(subReads.begin(), subReads.end(), readID), readID);
}

Read & Read::mergeWithThis(Read const & other, int ignoreConflicts, bool ignoreOverlap, bool ignoreContigID)
{
    if (!ignoreContigID && (contigID == INVALID || other.contigID == INVALID || contigID != other.contigID))
        throw std::runtime_error("Cannot merge: Incompatible read contigs.");
    if (!ignoreOverlap && (other.beginPos >= endPos || beginPos >= other.endPos))
        throw std::runtime_error("Cannot merge: Reads do not overlap!");

    beginPos = std::min(beginPos, other.beginPos);
    endPos = std::max(endPos, other.endPos);
    if (beginPos == -1)
        beginPos = other.beginPos;
    if (endPos == -1)
        endPos = other.endPos;

    // Merge subReads.
    {
        decltype(subReads) tmp;
        std::set_union(subReads.begin(), subReads.end(), other.subReads.begin(), other.subReads.end(),
                       std::back_inserter(tmp));
        swap(tmp, subReads);
    }

    // Merge features.
    features.mergeWithThis(other.features, ignoreConflicts);

    return *this;
}

std::ostream & operator<<(std::ostream & out, Read const & read)
{
    out << "Read(id=" << read.id << ", contigID=" << read.contigID
        << ", beginPos=" << read.beginPos << ", endPos=" << read.endPos
        << ", subReads (#=" << read.subReads.size() << ")={";
    for (unsigned i = 0; i < read.subReads.size(); ++i)
    {
        if (i > 0)
            out << ", ";
        out << read.subReads[i];
    }
    out << "}, features=" << read.features << ")";
    return out;
}

bool ltID(Read const & lhs, Read const & rhs)
{
    return (lhs.id < rhs.id);
}

bool ltContigID(Read const & lhs, Read const & rhs)
{
    return (lhs.contigID < rhs.contigID);
}

bool ltCoordinate(Read const & lhs, Read const & rhs)
{
    return (std::make_tuple(lhs.contigID, lhs.beginPos, -lhs.endPos) <
            std::make_tuple(rhs.contigID, rhs.beginPos, -rhs.endPos));
}

}  // namespace rep_sep
