/*
  Extensions for the naive hamming finder.
*/

#ifndef FIND_HAMMING_SIMPLE_EXT_H_
#define FIND_HAMMING_SIMPLE_EXT_H_

// Set the end position of the pattern in the finder.
template <typename THaystack, typename TNeedle, typename TPosition>
inline bool setEndPosition(Finder<THaystack, void> & finder,
                           Pattern<TNeedle, HammingSimple> & pattern,
                           const TPosition & pos) {
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
    while ((result = find(finder, pattern)) && endPosition(finder) < pos) {
//         std::cerr << "Skipping over end pos " << endPosition(finder) << std::endl;
        continue;
    }
//     std::cerr << "XXX beginPosition(finder) == " << beginPosition(finder) << std::endl;
//     std::cerr << "XXX endPosition(finder) == " << endPosition(finder) << std::endl;
    return result;
}

#endif  // FIND_HAMMING_SIMPLE_EXT_H_
