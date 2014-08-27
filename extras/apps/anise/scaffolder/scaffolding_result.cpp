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

#include "scaffolding_result.h"

#include <algorithm>
#include <iostream>

namespace scaffolder {

// ----------------------------------------------------------------------------
// Class ScaffoldingResult
// ----------------------------------------------------------------------------

void PositionedContig::print(std::ostream & out) const
{
    out << "PositionedContig(pos=" << pos << ", posSD=" << posSD << ", id="
        << id << ", length=" << length << ")";
}

bool PositionedContig::tentativeOverlap(PositionedContig const & other, int k) const
{
    return !(pos + k * posSD + (int)length < other.pos - k * posSD);
}

std::ostream & operator<<(std::ostream & out, PositionedContig const & pc)
{
    pc.print(out);
    return out;
}

// ----------------------------------------------------------------------------
// Class ScaffoldingResult
// ----------------------------------------------------------------------------

void ScaffoldingResult::print(std::ostream & out) const
{
    out << "SCAFFOLDING RESULT\n";
    unsigned i = 0;
    for (auto const & scaffold : scaffolds)
    {
        std::cerr << i++ << "\t";
        for (auto const & posC : scaffold)
            std::cerr << posC << " ";
            std::cerr << "\n";
    }
}

// Shifts the scaffolds to position 0.
void ScaffoldingResult::shiftScaffolds()
{
    // Shift all scaffolds so they start at 0.
    for (auto & scaffold : scaffolds)
    {
        auto minPos = std::min_element(scaffold.begin(), scaffold.end(),
                                       [](PositionedContig lhs, PositionedContig rhs) {
                                           return lhs.pos < rhs.pos;
                                       })->pos;
        std::for_each(scaffold.begin(), scaffold.end(),
                      [minPos](PositionedContig & x) { x.pos -= minPos; });
    }
}

}  // namespace scaffolder
