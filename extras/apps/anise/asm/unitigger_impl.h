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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_UNITIGGER_IMPL_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_ASM_UNITIGGER_IMPL_H_

#include <iosfwd>
#include <map>
#include <set>
#include <vector>

#include <seqan/sequence.h>

#include "asm/contig_graph.h"

namespace assembler {

// ============================================================================
// Forwards
// ============================================================================

class Overlap;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class UnitiggerImpl
// ----------------------------------------------------------------------------

// Abstract base class for unitigger implementations.

class UnitiggerImpl
{
public:
    virtual std::unique_ptr<ContigGraph> run(seqan::StringSet<seqan::Dna5String> const & seqs,
                                             std::vector<Overlap> const & overlaps) const = 0;
};

// ----------------------------------------------------------------------------
// Class ApproximateReadLayout
// ----------------------------------------------------------------------------

struct ApproximateReadLayout
{
    // TODO(holtgrew): Use Contig::PositionedRead?
    struct LayoutEntry
    {
        LayoutEntry() = default;

        LayoutEntry(unsigned readID, int pos, int len) : readID(readID), pos(pos), len(len)
        {}

        unsigned readID { (unsigned)-1 };
        int pos { -1 };
        int len { 0 };

        bool operator<(LayoutEntry const & other) const
        {
            return makeTuple() < other.makeTuple();
        }

        void print(std::ostream & out) const;
        
    private:

        std::tuple<int, int, unsigned> makeTuple() const
        {
            return std::make_tuple(pos, len, readID);
        }
    };

    void add(LayoutEntry const & entry);

    // Shift all positions to 0.
    void shift();

    int readPos(unsigned readID) const;

    // Convert to Contig using the approximate layout information.
    //
    // NB: The generated overlaps will not be assigned an ID since the overlap information is only approximate.
    Contig toContig(unsigned contigID) const;

    void print(std::ostream & out) const;

    // Entries, sorted by position.
    std::set<LayoutEntry> entries;
    // Positions by readID.
    std::map<unsigned, int> positions;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace assembler

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_ASM_UNITIGGER_IMPL_H_
