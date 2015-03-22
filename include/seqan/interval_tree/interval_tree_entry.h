// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Facade header for module interval_tree.
// ==========================================================================

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_TREE_ENTRY_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_TREE_ENTRY_H_

#include <seqan/basic.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class IntervalTreeEntry
// ----------------------------------------------------------------------------

template <typename TValue>
class IntervalTreeEntry
{
    typedef typename Position<TValue>::Type TPosition_;

public:

    // The value wrapped in the entry. Must follow the IntervalConcept concept.
    TValue value;
    // The maximal end position of this nodes and all of its children.
    TPosition_ maxEnd;

    IntervalTreeEntry() : value(), maxEnd() {}

    IntervalTreeEntry(TValue value, TPosition_ maxEnd) : value(value), maxEnd(maxEnd) {}

    template <typename TPos> bool _allLeftOf(TPos point) const
    {
        return (maxEnd <= (TPosition_)point);
    }

    template <typename TPos> bool _isLeftOf(TPos point) const
    {
        return (endPos(value) <= (TPosition_)point);
    }

    template <typename TPos> bool _isRightOf(TPos point) const
    {
        return ((TPosition_)point < beginPos(value));
    }

    template <typename TPos> bool _contains(TPos point) const
    {
        return (beginPos(value) <= (TPosition_)point && (TPosition_)point < endPos(value));
    }

    template <typename TPos> bool _overlapsWith(TPos posBegin, TPos posEnd) const
    {
        return ((TPosition_)posBegin < endPos(value) && beginPos(value) < (TPosition_)posEnd);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/**
 * @mfn IntervalTreeEntry#Position
 * @brief The position type of the <tt>IntervalTreeEntry</tt>, forwards to <tt>Position&lt;TValue&gt;</tt>.
 *
 * @signature Position<TIntervalTreeEntry>::Type;
 *
 * @tparam TIntervalTreeEntry The IntervalTreeEntry to query for its position type.
 */
template <typename TValue>
struct Position<IntervalTreeEntry<TValue> >
{
    typedef typename Position<TValue>::Type Type;
};

template <typename TValue>
struct Position<IntervalTreeEntry<TValue> const>
{
    typedef typename Position<TValue const>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_TREE_ENTRY_H_
