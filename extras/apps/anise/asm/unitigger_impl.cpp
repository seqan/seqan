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

#include "unitigger_impl.h"

#include <seqan/basic.h>

#include <iostream>

namespace assembler {

// ----------------------------------------------------------------------------
// Class ApproximateReadLayout::LayoutEntry
// ----------------------------------------------------------------------------

void ApproximateReadLayout::LayoutEntry::print(std::ostream & out) const
{
    out << "ApproximateReadLayout::LayoutEntry(readID=" << readID
        << ", pos=" << pos << ", len=" << len << ")";
}

// ----------------------------------------------------------------------------
// Class ApproximateReadLayout
// ----------------------------------------------------------------------------

void ApproximateReadLayout::add(LayoutEntry const & entry)
{
    entries.insert(entry);
    positions[entry.readID] = entry.pos;
}

void ApproximateReadLayout::shift()
{
    if (entries.empty())
        return;
    int shift = entries.begin()->pos;

    decltype(entries) tmp;
    std::transform(entries.begin(), entries.end(), std::inserter(tmp, tmp.end()),
                   [shift](LayoutEntry entry) { entry.pos -= shift; return entry; });
    for (auto const pair : positions)
        positions[pair.first] = pair.second - shift;
    swap(tmp, entries);
}

int ApproximateReadLayout::readPos(unsigned readID) const
{
    SEQAN_ASSERT_MSG(positions.count(readID), "readID=%u", readID);
    return positions.find(readID)->second;
}

// Convert to Contig using the approximate layout information.
//
// NB: The generated overlaps will not be assigned an ID since the overlap information is only approximate.
Contig ApproximateReadLayout::toContig(unsigned contigID) const
{
    Contig result;
    result.id = contigID;

    std::transform(entries.begin(), entries.end(), std::back_inserter(result.posReads),
                   [](LayoutEntry const & entry) {
                       return Contig::PositionedRead(entry.readID, entry.pos, entry.len);
                   });
    return result;
}

void ApproximateReadLayout::print(std::ostream & out) const
{
    out << "APPROXIMATE READ LAYOUT\n\n";
    for (auto const & entry : entries)
    {
        entry.print(std::cerr);
        std::cerr << "\n";
    }
    out << "/END OF APPROXIMATE READ LAYOUT\n";
}

}  // namespace assembler
