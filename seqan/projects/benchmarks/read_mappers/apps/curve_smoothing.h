#ifndef WIT_BUILDER_CURVE_SMOOTHING_H_
#define WIT_BUILDER_CURVE_SMOOTHING_H_

#include <algorithm>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

#include "witio.h"

using namespace seqan;

/*
  Fills gaps in the given error curve.  If there are two adjacent
  entries in the string that have the same start position but not
  adjacent end positions then we add an entry for each end position
  in between.

  errorCurve must be sorted by (contigId, endPosition).
 */
void fillGaps(String<WeightedMatch> & errorCurve) {
    typedef Iterator<String<WeightedMatch>, Standard>::Type TIterator;

    // Guard against not enough items in string.
    if (length(errorCurve) < 2)
        return;

    // The insertion buffer.
    String<WeightedMatch> buffer;
    for (TIterator it = begin(errorCurve, Standard()); (it + 1) != end(errorCurve, Standard()); ++it) {
        if (value(it).contigId != value(it + 1).contigId)
            continue;  // Skip if contigId is not equal.
        if (value(it).isForward != value(it + 1).isForward)
            continue;  // Skip if direction is not equal.
        if (value(it).pos + 1 == value(it + 1).pos)
            continue;  // Skip if the end positions are contiguous.
        if (value(it).beginPos != value(it + 1).beginPos)
            continue;  // Skip if the begin positions are not equal.
        // If the code reaches this position then we have a gap and fill it.
        int gapScore = _min(value(it).distance, value(it + 1).distance);
        for (size_t i = value(it).pos + 1; i < value(it + 1).pos; ++i)
            appendValue(buffer, WeightedMatch(value(it).contigId, value(it).isForward, i, gapScore, value(it).beginPos));
    }

    // Append buffer and sort inserted items to the right position.
    append(errorCurve, buffer);
    std::sort(begin(errorCurve), end(errorCurve));
}


/*
  The algorithm works as follows.  For each contiguous interval of
  tuples with the same start position, do the following.  From the
  left to the right and from the right to the left to the first
  absolute maximum in each direction, make the curve monotonously
  increasing.  Between the rightmost and the leftmost maximum, set the
  values to the maximal value.
 */
void smoothErrorCurve(String<WeightedMatch> &errorCurve)
{
    if (length(errorCurve) == 0u) return;
    // TODO(holtgrew): Standard should be the standard iterator, yes?
    typedef Iterator<String<WeightedMatch>, Standard>::Type TIterator;
    typedef Position<String<WeightedMatch> >::Type TPosition;

    // For all intervals with the same beginPosition.
    WeightedMatchBeginPosNeqOrContigIdNeq pred;
    TIterator itBegin = begin(errorCurve, Standard());
    TIterator itEnd = std::adjacent_find(itBegin, end(errorCurve, Standard()), pred);
    do {
//         std::cerr << "begin is " << (itBegin - begin(errorCurve, Standard())) << std::endl;
//         std::cerr << "end is " << (itEnd - begin(errorCurve, Standard())) << std::endl;
        if (itEnd != end(errorCurve, Standard()))
            itEnd += 1;
        SEQAN_ASSERT_NEQ(itBegin, itEnd);

        if (itBegin != itEnd) {
            // Find maximal score value.
            int maxScore = value(itBegin).distance;
            for (TIterator it = itBegin; it != itEnd; ++it)
                maxScore = _max(maxScore, value(it).distance);

//             std::cerr << "--- BEFORE itBegin..itEnd" << std::endl;
//             for (TIterator it = itBegin; it != itEnd; ++it) {
//                 std::cerr << value(it) << std::endl;
//             }
//             std::cerr << "---" << std::endl;

//             std::cerr << "maxScore is " << maxScore << std::endl;

            // Search leftmost and rightmost maximum.
            TIterator itLeftmost = 0, itRightmost = 0;
            for (TIterator it = itBegin; it != itEnd; ++it) {
                if (value(it).distance != maxScore)
                    continue;
                if (itLeftmost == 0)
                    itLeftmost = it;
                itRightmost = it;
            }

//             std::cerr << "leftmost is " << (itLeftmost - begin(errorCurve, Standard())) << " (" << value(itLeftmost).distance << ")" << std::endl;
//             std::cerr << "rightmost is " << (itRightmost - begin(errorCurve, Standard())) << " (" << value(itRightmost).distance << ")" << std::endl;

            // Level between leftmost and rightmost maxima.
            for (TIterator it = itLeftmost; it != itRightmost + 1; ++it)
                value(it).distance = maxScore;

            // Make monotonously decreasing from begin to leftmost.
            int currentMax = value(itBegin).distance;
            for (TIterator it = itBegin; it != itLeftmost; ++it) {
                currentMax = _max(currentMax, value(it).distance);
//                 std::cerr << "value(" << (it - begin(errorCurve, Standard())) << ") = " << currentMax << std::endl;
                value(it).distance = currentMax;
            }

            // Make monotously decreasing from end to rightmost.
            TPosition offset = itBegin - begin(errorCurve, Standard());
            TPosition rightmostPos = offset + itRightmost - itBegin;
            TPosition intervalLength = itEnd - itBegin;
            typedef ModifiedString<Segment<String<WeightedMatch>, InfixSegment>, ModReverse> TModifiedString;
            TModifiedString rRightInterval(infix(errorCurve, rightmostPos, offset + intervalLength));
            typedef Iterator<TModifiedString>::Type TModifiedStringIterator;
            currentMax = front(rRightInterval).distance;
            for (TModifiedStringIterator it = begin(rRightInterval, Standard()); it != end(rRightInterval, Standard()); ++it) {
                currentMax = _max(currentMax, value(it).distance);
//                 std::cerr << "value(" << value(it) << ") = " << currentMax << std::endl;
                value(it).distance = currentMax;
            }
        }

//         std::cerr << "--- AFTER itBegin..itEnd" << std::endl;
//         for (TIterator it = itBegin; it != itEnd; ++it) {
//             std::cerr << value(it) << std::endl;
//         }
//         std::cerr << "---" << std::endl;

        // Proceed.
        itBegin = itEnd;
        itEnd = std::adjacent_find(itBegin, end(errorCurve, Standard()), pred);
    } while (itBegin != end(errorCurve));
}

#endif  // WIT_BUILDER_CURVE_SMOOTHING_H_
