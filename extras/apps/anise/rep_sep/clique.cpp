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

#include "clique.h"

#include <algorithm>
#include <iostream>

#include <boost/dynamic_bitset.hpp>

#include <seqan/basic.h>

namespace rep_sep {

// --------------------------------------------------------------------------
// Class Clique
// --------------------------------------------------------------------------

unsigned Clique::intersectionSize(boost::dynamic_bitset<> const & otherBitSet) const
{
    SEQAN_ASSERT_EQ(bitSet.size(), otherBitSet.size());
    return (bitSet & otherBitSet).count();
}

void Clique::intersection(boost::dynamic_bitset<> & result,
                          boost::dynamic_bitset<> const & otherBitSet) const
{
    SEQAN_ASSERT_EQ(bitSet.size(), otherBitSet.size());
    result = bitSet;
    result &= otherBitSet;
}

void Clique::addRead(Read const & read)
{
    SEQAN_CHECK(read.contigID == contigID, "contigs must match");
    beginPos = std::min(beginPos, read.beginPos);
    endPos = std::max(endPos, read.endPos);
    if (beginPos == -1)
        beginPos = read.beginPos;
    if (endPos == -1)
        endPos = read.endPos;
    _size += !bitSet[read.id];
    SEQAN_ASSERT_LT(read.id, bitSet.size());
    bitSet.set(read.id, true);
}

bool Clique::includes(Clique const & other) const
{
    if (size() < other.size())
        return false;  // too small to subsume
    return (other.bitSet.is_subset_of(bitSet));
}

std::ostream & operator<<(std::ostream & out, Clique const & clique)
{
    out << "Clique(contigID=" << clique.contigID << ", beginPos=" << clique.beginPos
        << ", endPos=" << clique.endPos << ", readIDs (#=" << clique.bitSet.size() << ")={";
    for (unsigned i = 0; i < clique.bitSet.size(); ++i)
        if (clique.bitSet[i])
            out << i << ", ";
    out << "})";
    return out;
}

}  // namespace rep_sep
