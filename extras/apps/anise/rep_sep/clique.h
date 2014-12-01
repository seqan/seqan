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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_H_

#include "rep_sep/read.h"
#include "rep_sep/read_set.h"
#include "rep_sep/string_packed_pop_count.h"

#include <boost/dynamic_bitset.hpp>

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Clique
// ----------------------------------------------------------------------------

// An entry in CliqueSet.

struct Clique
{
    static const unsigned INVALID = (unsigned)-1;

    Clique() = default;

    explicit Clique(unsigned numReads) : Clique()
    {
        bitSet.resize(numReads);
    }

    Clique(unsigned contigID, int beginPos, int endPos, unsigned numReads) :
            contigID(contigID), beginPos(beginPos), endPos(endPos)
    {
        bitSet.resize(numReads);
    }

    Clique(Read const & read, unsigned numReads) :
            contigID(read.contigID), beginPos(read.beginPos), endPos(read.endPos)
    {
        bitSet.resize(numReads);
        SEQAN_ASSERT_LT(read.id, numReads);
        bitSet.set(read.id);
        _size = 1;
    }

    Clique(unsigned contigID,
           boost::dynamic_bitset<> const & bitSet,
           FeatureReadSet const & featureReadSet) :
            contigID(contigID), bitSet(bitSet)
    {
        _size = bitSet.count();
        if (!size())
            return;

        beginPos = seqan::maxValue<int>();
        endPos = seqan::minValue<int>();
        for (unsigned i = 0; i < bitSet.size(); ++i)
            if (bitSet[i])
            {
                beginPos = std::min(beginPos, featureReadSet.reads[i].beginPos);
                endPos = std::max(endPos, featureReadSet.reads[i].endPos);
            }
    }

    bool includes(Clique const & other) const;

    // Returns intersection size of readIDs with otherReadIDs.
    unsigned intersectionSize(boost::dynamic_bitset<> const & otherBitSet) const;
    // Returns inserection of readIDs with otherReadIDs.
    void intersection(boost::dynamic_bitset<> & result,
                      boost::dynamic_bitset<> const & otherBitSet) const;

    // Insert new read into clique.
    void addRead(Read const & read);

    // Returns number of reads in clique.
    size_t size() const { return _size; }

    // Genomic region.
    unsigned contigID { INVALID };
    int beginPos { -1 };
    int endPos { -1 };
    // Unsigned, size.
    unsigned _size { 0 };
    // Bit vector with read set.
    boost::dynamic_bitset<> bitSet;
};

std::ostream & operator<<(std::ostream & out, Clique const & clique);

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_H_
