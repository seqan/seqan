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

#include "mate_link.h"

#include <ostream>

namespace scaffolder {

// ----------------------------------------------------------------------------
// Class MateEdgeLabel
// ----------------------------------------------------------------------------

std::ostream & MateEdgeLabel::print(std::ostream & out) const
{
    return out << "MateEdgeLabel(lengthMean=" << lengthMean
               << ", lengthStdDev=" << lengthStdDev
               << ", count=" << count
               << ", weight=" << weight << ")";
}

std::ostream & operator<<(std::ostream & out, MateEdgeLabel const & label)
{
    return label.print(out);
}

// ----------------------------------------------------------------------------
// Class MateLink
// ----------------------------------------------------------------------------

std::ostream & MateLink::print(std::ostream & out) const
{
    return out << "MateLink(" << source << ", " << target << ", " << label << ")";
}

std::ostream & operator<<(std::ostream & out, MateLink const & link)
{
    return link.print(out);
}

// ----------------------------------------------------------------------------
// Class ContigEdgeLabel
// ----------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & out, ContigEdgeLabel const & label)
{
    return out << "ContigEdgeLabel(" << label.length << ")";
}

}  // namespace scaffolder
