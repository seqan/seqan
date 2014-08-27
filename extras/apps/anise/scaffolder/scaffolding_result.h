// ==========================================================================
//                                  ANISE
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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_SCAFFOLDING_RESULT_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_SCAFFOLDING_RESULT_H_

#include <iosfwd>
#include <vector>
#include <utility>

namespace scaffolder {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class ScaffoldingResult
// ----------------------------------------------------------------------------

// The result of the path merging algorithm consists of a list of scaffolds.  Each scaffold is a list of contigs (given
// as ids) with a begin position on the line.  Each scaffold starts at position 0, the positions indicate overlaps.

struct PositionedContig
{
    // Position of the contig.
    int pos;
    // Standard deviation in position.
    int posSD;
    // ID of the contig.
    unsigned id;
    // Length of the contig.
    unsigned length;

    explicit
    PositionedContig(int pos, int posSD, unsigned id, unsigned length) :
            pos(pos), posSD(posSD), id(id), length(length)
    {}

    void print(std::ostream & out) const;

    bool operator<(PositionedContig const & other) const
    {
        return std::make_pair(pos, id) < std::make_pair(other.pos, id);
    }

    bool operator==(PositionedContig const & other) const
    {
        return std::make_pair(pos, id) == std::make_pair(other.pos, id);
    }

    // Returns true if overlaps tentatively with other.
    bool tentativeOverlap(PositionedContig const & other, int k = 3) const;
};

std::ostream & operator<<(std::ostream & out, PositionedContig const & pc);

// ----------------------------------------------------------------------------
// Class ScaffoldingResult
// ----------------------------------------------------------------------------

struct ScaffoldingResult
{
    // A Scaffold is a list [(position, contigID)].
    typedef std::vector<PositionedContig> TScaffold;

    // The result of the scaffolding is a list of scaffolds.x
    std::vector<TScaffold> scaffolds;

    void clear()
    {
        scaffolds.clear();
    }

    void print(std::ostream & out) const;

    // Shifts the scaffolds to position 0.
    void shiftScaffolds();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace scaffolder

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_SCAFFOLDER_SCAFFOLDING_RESULT_H_
