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
// Definition of the class WeightedMatch and curve smoothing operations on
// strings of WeightedMatch objects.
// ==========================================================================

#ifndef SEQAN_APPS_RABEMA_CURVE_SMOOTHING_H_
#define SEQAN_APPS_RABEMA_CURVE_SMOOTHING_H_

#include <algorithm>

#include <seqan/modifier.h>
#include <seqan/sequence.h>

using namespace seqan;

// ============================================================================
// Enums, Tags, Classes, Typedefs.
// ============================================================================

// ----------------------------------------------------------------------------
// Class WeightedMatch
// ----------------------------------------------------------------------------

// Represents a read match to a contig with an error.
//
// WeightedMatch(c, p, d) represents a match of a read at position p in contig c with distance d.

class WeightedMatch
{
public:
    size_t contigId;
    bool isForward;
    size_t pos;
    int distance;
    // Optional begin position, used in wit builder to smooth matches.
    // Is not written out or used in the less than operator.
    size_t beginPos;

    WeightedMatch() {}

    WeightedMatch(size_t _contigId, bool _isForward, size_t _pos, int _distance, size_t _beginPos) :
        contigId(_contigId), isForward(_isForward), pos(_pos), distance(_distance), beginPos(_beginPos)
    {}

    // Order lexicographically by (contigId, pos, -distance).
    bool operator<(WeightedMatch const & other) const
    {
        if (contigId < other.contigId) return true;

        if (contigId == other.contigId && isForward > other.isForward) return true;

        if (contigId == other.contigId && isForward == other.isForward &&
            pos < other.pos) return true;

        if (contigId == other.contigId && isForward == other.isForward &&
            pos == other.pos && distance > other.distance) return true;

        return false;
    }

    bool operator==(WeightedMatch const & other) const
    {
        if (contigId != other.contigId) return false;

        if (isForward != other.isForward) return false;

        if (pos != other.pos) return false;

        if (distance != other.distance) return false;

        if (beginPos != other.beginPos) return false;

        return true;
    }

};

typedef String<WeightedMatch> TWeightedMatches;

// ----------------------------------------------------------------------------
// Helper Class WeightedMatchBeginPosNeqOrContigIdNeq
// ----------------------------------------------------------------------------

struct WeightedMatchBeginPosNeqOrContigIdNeq :
    std::binary_function<WeightedMatch, WeightedMatch, bool>
{
    bool operator()(WeightedMatch const & arg1, WeightedMatch const & arg2)
    {
        return arg1.beginPos != arg2.beginPos || arg1.contigId != arg2.contigId || arg1.isForward != arg2.isForward;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator<<()                                        [WeightedMatch]
// ----------------------------------------------------------------------------

// Stream output for WeightedMatch objects, for debugging.
template <typename TStream>
TStream & operator<<(TStream & out, WeightedMatch const & m)
{
    out << "(" << m.contigId << ", " << (m.isForward ? "F, " : "R, ")
    << m.pos << ", " << m.distance << ", " << m.beginPos << ")";
    return out;
}

// ----------------------------------------------------------------------------
// Function fillGaps()
// ----------------------------------------------------------------------------

// Fills gaps in the given error curve.  If there are two adjacent entries in the string that have the same start
// position but not adjacent end positions then we add an entry for each end position in between.

// errorCurve must be sorted by (contigId, endPosition).

void fillGaps(String<WeightedMatch> & errorCurve)
{
    typedef Iterator<String<WeightedMatch>, Standard>::Type TIterator;

    // Guard against not enough items in string.
    if (length(errorCurve) < 2)
        return;

    // The insertion buffer.
    String<WeightedMatch> buffer;
    for (TIterator it = begin(errorCurve, Standard()); (it + 1) != end(errorCurve, Standard()); ++it)
    {
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

// ----------------------------------------------------------------------------
// Function smoothErrorCurve()
// ----------------------------------------------------------------------------

// Smoothing of error curves, for the application of match equivalence.
//
// The algorithm works as follows.  For each contiguous interval of tuples with the same start position, do the
// following.  From the left to the right and from the right to the left to the first absolute maximum in each
// direction, make the curve monotonously increasing.  Between the rightmost and the leftmost maximum, set the values to
// the maximal value.

void smoothErrorCurve(String<WeightedMatch> & errorCurve)
{
    if (length(errorCurve) == 0u)
        return;

    // TODO(holtgrew): Standard should be the standard iterator, yes?
    typedef Iterator<String<WeightedMatch>, Standard>::Type TIterator;
    typedef Position<String<WeightedMatch> >::Type TPosition;

    // For all intervals with the same beginPosition.
    WeightedMatchBeginPosNeqOrContigIdNeq pred;
    TIterator itBegin = begin(errorCurve, Standard());
    TIterator itEnd = std::adjacent_find(itBegin, end(errorCurve, Standard()), pred);
    do
    {
//         std::cerr << "begin is " << (itBegin - begin(errorCurve, Standard())) << std::endl;
//         std::cerr << "end is " << (itEnd - begin(errorCurve, Standard())) << std::endl;
        if (itEnd != end(errorCurve, Standard()))
            itEnd += 1;
        SEQAN_ASSERT_NEQ(itBegin, itEnd);

        if (itBegin != itEnd)
        {
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
            for (TIterator it = itBegin; it != itEnd; ++it)
            {
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
            for (TIterator it = itBegin; it != itLeftmost; ++it)
            {
                currentMax = _max(currentMax, value(it).distance);
//                 std::cerr << "value(" << (it - begin(errorCurve, Standard())) << ") = " << currentMax << std::endl;
                value(it).distance = currentMax;
            }

            // Make monotonously decreasing from end to rightmost.
            TPosition offset = itBegin - begin(errorCurve, Standard());
            TPosition rightmostPos = offset + itRightmost - itBegin;
            TPosition intervalLength = itEnd - itBegin;
            typedef ModifiedString<Infix<String<WeightedMatch> >::Type, ModReverse> TModifiedString;
            TModifiedString rRightInterval(infix(errorCurve, rightmostPos, offset + intervalLength));
            typedef Iterator<TModifiedString, Standard>::Type TModifiedStringIterator;
            currentMax = front(rRightInterval).distance;
            for (TModifiedStringIterator it = begin(rRightInterval, Standard()); it != end(rRightInterval, Standard()); ++it)
            {
                currentMax = _max(currentMax, value(it).distance);
//                 std::cerr << "value(" << value(it) << ") = " << currentMax << std::endl;
                value(host(it)).distance = currentMax;
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
    }
    while (itBegin != end(errorCurve));
}

#endif  // SEQAN_APPS_RABEMA_CURVE_SMOOTHING_H_
