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
// Extensions for the myers ukkonen pattern.
// ==========================================================================

#ifndef SEQAN_APPS_RABEMA_FIND_MYERS_UKKONEN_EXT_H_
#define SEQAN_APPS_RABEMA_FIND_MYERS_UKKONEN_EXT_H_

#include <seqan/find.h>

using namespace seqan;

// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(Finder<THaystack, void> & finder,
                           Pattern<TNeedle, Myers<FindInfix> > & pattern,
                           const TPosition & pos)
{
    // Compute delta, such that we start searching at pos - delta.
    TPosition delta = pattern.needleSize + _min(pattern.needleSize, static_cast<size_t>(-scoreLimit(pattern))) + 1;
    if (delta > pos)
        delta = pos;

    // Set end position in the finder to pos - delta.
    setPosition(finder, pos - delta);

    // Clear the pattern, and search until we are at pos.
    _patternInit(pattern, finder);
    bool result;
    while ((result = find(finder, pattern)) &&
           endPosition(finder) < pos)
        continue;
    return result;
}

#endif  // SEQAN_APPS_RABEMA_FIND_MYERS_UKKONEN_EXT_H_
