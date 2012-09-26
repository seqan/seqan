/*
  Extensions for the myers ukkonen pattern.
*/

#ifndef FIND_MYERS_UKKONEN_EXT_H_
#define FIND_MYERS_UKKONEN_EXT_H_

#include <seqan/find.h>

using namespace seqan;

// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(Finder<THaystack, void> & finder,
                           Pattern<TNeedle, Myers<FindInfix> > & pattern,
                           const TPosition & pos) {
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

#endif  // FIND_MYERS_UKKONEN_EXT_H_
