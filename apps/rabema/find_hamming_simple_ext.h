// ==========================================================================
//                      RABEMA Read Alignment Benchmark
// ==========================================================================
// Copyright (C) 2010-1012 Manuel Holtgrewe, FU Berlin
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Extensions for the naive hamming finder.
// ==========================================================================

#ifndef SEQAN_APPS_RABEMA_FIND_HAMMING_SIMPLE_EXT_H_
#define SEQAN_APPS_RABEMA_FIND_HAMMING_SIMPLE_EXT_H_

#include <seqan/find.h>

// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(seqan::Finder<THaystack, void> & finder,
                           seqan::Pattern<TNeedle, seqan::HammingSimple> & pattern,
                           const TPosition & pos)
{
//     std::cerr << "setEndPosition(finder, pattern, " << pos << ")" << std::endl;
    // Compute delta, such that we start searching at pos - delta.
    TPosition delta = length(needle(pattern));
    if (delta > pos)
        delta = pos;
//     std::cerr << "delta == " << delta << std::endl;

    // Set end position in the finder to pos - delta.
    finder.data_length = length(needle(pattern));
    setPosition(finder, pos - delta);
    finder.data_endPos = pos - delta;
//     std::cerr << "beginPosition(finder) == " << beginPosition(finder) << std::endl;
//     std::cerr << "endPosition(finder) == " << endPosition(finder) << std::endl;

    // Clear the pattern, and search until we are at pos.
    bool result;
    while ((result = find(finder, pattern)) && endPosition(finder) < pos)
    {
//         std::cerr << "Skipping over end pos " << endPosition(finder) << std::endl;
        continue;
    }
//     std::cerr << "XXX beginPosition(finder) == " << beginPosition(finder) << std::endl;
//     std::cerr << "XXX endPosition(finder) == " << endPosition(finder) << std::endl;
    return result;
}

#endif  // SEQAN_APPS_RABEMA_FIND_HAMMING_SIMPLE_EXT_H_
