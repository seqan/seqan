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
// Static, immutable version of interval tree.
// ==========================================================================

// TODO(holtgrewe): Implement neighbor search in case of no overlapping intervals.

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_H_

#include <algorithm>
#include <vector>
#include <utility>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include "base_interval_tree.h"
#include "interval_tree_entry.h"
#include "static_interval_tree_iterator.h"

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class StaticIntervalTree
// ----------------------------------------------------------------------------

// Tags for specializing IntervalTree<>.

struct Static_;
typedef Tag<Static_> Static;

/*!
 * @class StaticIntervalTree
 * @headerfile <seqan/interval_tree.h>
 * @brief Static @link IntervalTree @endlink implementation.
 *
 * @signature template <typename TValue>
 *            class IntervalTree<TValue, Static>;
 *
 * @tparam TValue Value to store in the tree.
 * @tparam TSpec  Specializing tag.
 *
 * The implementation is based on Cormen's augmented binary search trees as described in <a
 * href="http://en.wikipedia.org/wiki/Interval_tree#Augmented_tree">this Wikipedia article</a>.
 */
template <typename TValue>
class IntervalTree<TValue, Static>
{
public:
    IntervalTree() {}

    template <typename TContainer>
    IntervalTree(TContainer const & container)
    {
        initialize(container);
    }

    std::vector<IntervalTreeEntry<TValue> > entries;

private:

    template <typename TContainer> void initialize(TContainer const & container);

    typename Position<IntervalTreeEntry<TValue> >::Type
    computeMaxEndProperties(typename std::vector<TValue>::size_type beginIdx,
                            typename std::vector<TValue>::size_type endIdx);
};

template <typename TValue>
struct CmpIntervalTreeEntryByBeginPos_
{
    bool operator()(IntervalTreeEntry<TValue> const & lhs, IntervalTreeEntry<TValue> const & rhs)
    {
        return (beginPos(lhs.value) < beginPos(rhs.value) ||
                (beginPos(lhs.value) == beginPos(rhs.value) && endPos(lhs.value) < endPos(rhs.value)));
    }
};

template <typename TValue>
template <typename TContainer>
void IntervalTree<TValue, Static>::initialize(TContainer const & container)
{
    // Fill this->entries[i].value from the values in container.
    entries.resize(length(container));
    typename std::vector<IntervalTreeEntry<TValue> >::iterator entryIt = entries.begin();
    for (typename Iterator<TContainer const, Standard>::Type it = begin(container, Standard());
         it != end(container, Standard()); ++it, ++entryIt)
    {
        entryIt->value = *it;
        entryIt->maxEnd = endPos(*it);
    }

    // Sort this->entries by the begin position of their value member.
    std::sort(begin(entries, Standard()), end(entries, Standard()),
              CmpIntervalTreeEntryByBeginPos_<TValue>());

    // Compute the maxEnd members of this->entries.
    typedef typename std::vector<TValue>::size_type TSize;
    computeMaxEndProperties((TSize)0, (TSize)length(entries));
}

template <typename TValue>
typename Position<IntervalTreeEntry<TValue> >::Type
IntervalTree<TValue, Static>::computeMaxEndProperties(
        typename std::vector<TValue>::size_type beginIdx,
        typename std::vector<TValue>::size_type endIdx)
{
    typedef typename Position<IntervalTreeEntry<TValue> >::Type TPos;
    typedef typename std::vector<TValue>::size_type             TSize;

    if (beginIdx == endIdx)
        return MaxValue<TPos>::VALUE;

    TSize centerIdx = (endIdx + beginIdx) / 2;
    IntervalTreeEntry<TValue> & entry = entries[centerIdx];

    if (beginIdx + 1 == endIdx)
        return entry.maxEnd;

    entry.maxEnd = std::max(entry.maxEnd,
                            std::max(computeMaxEndProperties(beginIdx, centerIdx),
                                     computeMaxEndProperties(centerIdx, endIdx)));
    return entry.maxEnd;
}

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/**
 * @mfn IntervalTree#Size
 * @brief The size type of the <tt>IntervalTree</tt>.
 */

template <typename TValue>
struct Size<IntervalTree<TValue, Static> >
{
    typedef typename std::vector<TValue>::size_type Type;
};

template <typename TValue>
struct Size<IntervalTree<TValue, Static> const>
{
    typedef typename std::vector<TValue>::size_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

/**
 * @mfn IntervalTree#Iterator
 * @brief The iterator type of the <tt>IntervalTree</tt>.
 */

template <typename TValue, typename TIterSpec>
struct Iterator<IntervalTree<TValue, Static>, TIterSpec>
{
    typedef Iter<IntervalTree<TValue, Static>, void> Type;
};

template <typename TValue, typename TIterSpec>
struct Iterator<IntervalTree<TValue, Static> const, TIterSpec>
{
    typedef Iter<IntervalTree<TValue, Static> const, void> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// TODO(holtgrewe): Get rid of non-const length once all such functions are removed
template <typename TValue>
typename Size<IntervalTree<TValue, Static> >::Type
length(IntervalTree<TValue, Static> const & tree)
{
    return tree.entries.size();
}

template <typename TValue>
typename Size<IntervalTree<TValue, Static> >::Type
length(IntervalTree<TValue, Static> & tree)
{
    return tree.entries.size();
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

// TODO(holtgrewe): For some reason, using TIterSpec template parameter instead of Standard/Rooted created ambiguous function calls and references to missing class specialization of IteratorDefaultImp_, same below for end().

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type
begin(IntervalTree<TValue, Static> const & tree, Standard const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type(tree.entries.begin());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static>, Standard>::Type
begin(IntervalTree<TValue, Static> & tree, Standard const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static>, Standard>::Type(tree.entries.begin());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static> const, Rooted>::Type
begin(IntervalTree<TValue, Static> const & tree, Rooted const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type(tree.entries.begin());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static>, Rooted>::Type
begin(IntervalTree<TValue, Static> & tree, Rooted const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static>, Standard>::Type(tree.entries.begin());
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type
end(IntervalTree<TValue, Static> const & tree, Standard const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type(tree.entries.end());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static>, Standard>::Type
end(IntervalTree<TValue, Static> & tree, Standard const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static>, Standard>::Type(tree.entries.end());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static> const, Rooted>::Type
end(IntervalTree<TValue, Static> const & tree, Rooted const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static> const, Standard>::Type(tree.entries.end());
}

template <typename TValue>
typename Iterator<IntervalTree<TValue, Static>, Rooted>::Type
end(IntervalTree<TValue, Static> & tree, Rooted const & /*tag*/)
{
    return typename Iterator<IntervalTree<TValue, Static>, Standard>::Type(tree.entries.end());
}

// ----------------------------------------------------------------------------
// Function findOverlappingWithPoint()
// ----------------------------------------------------------------------------

template <typename TValue, typename TResult>
void _findOverlappingWithPoint(IntervalTree<TValue, Static> const & tree,
                               typename Size<IntervalTree<TValue, Static> const >::Type beginIdx,
                               typename Size<IntervalTree<TValue, Static> const >::Type endIdx,
                               typename Size<IntervalTree<TValue, Static> const >::Type centerIdx,
                               typename Position<IntervalTree<TValue, Static> const >::Type point,
                               TResult & result)
{
    typedef typename Iterator<IntervalTree<TValue, Static> const>::Type TIter;

    if (beginIdx >= endIdx) // handle base case of empty interval
        return;

    IntervalTreeEntry<TValue> const & entry = tree.entries[centerIdx];

    if (entry._allLeftOf(point)) // point is right of the rightmost point of any interval in this entry
        return;
    
    if (beginIdx < centerIdx) // recurse left
        _findOverlappingWithPoint(tree, beginIdx, centerIdx, beginIdx + (centerIdx - beginIdx) / 2,
                                  point, result);
    
    if (entry._contains(point)) // check this entry
        appendValue(result, TIter(tree.entries.begin() + centerIdx));
    
    if (entry._isRightOf(point)) // point is left of the start of the interval, can't to the right
        return;
    
    if (centerIdx + 1 < endIdx) // recurse right
        _findOverlappingWithPoint(tree, centerIdx + 1, endIdx, (centerIdx + 1) + (endIdx - (centerIdx + 1)) / 2,
                                  point, result);
}

// Guarantee: result is sorted by (beginPos, endPos), by search algorithm.

template <typename TValue, typename TPos, typename TResult>
void findOverlappingWithPoint(IntervalTree<TValue, Static> const & tree,
                              TPos point,
                              TResult & result)
{
    typedef IntervalTree<TValue, Static> const                     TIntervalTree;
    typedef typename Position<TIntervalTree>::Type                 TPos2;
    typedef typename Size<TIntervalTree>::Type                     TSize;
    typedef typename Iterator<TIntervalTree const, Standard>::Type TIter;

    _findOverlappingWithPoint(tree, (TSize)0, (TSize)tree.entries.size(), (TSize)(tree.entries.size() / 2), (TPos2)point, result);
}

// ----------------------------------------------------------------------------
// Function findOverlappingWithInterval()
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function findOverlappingWithPoint()
// ----------------------------------------------------------------------------

template <typename TValue, typename TResult>
void _findOverlappingWithInterval(IntervalTree<TValue, Static> const & tree,
                                  typename Size<IntervalTree<TValue, Static> const >::Type beginIdx,
                                  typename Size<IntervalTree<TValue, Static> const >::Type endIdx,
                                  typename Size<IntervalTree<TValue, Static> const >::Type centerIdx,
                                  typename Position<IntervalTree<TValue, Static> const >::Type posBegin,
                                  typename Position<IntervalTree<TValue, Static> const >::Type posEnd,
                                  TResult & result)
{
    typedef typename Iterator<IntervalTree<TValue, Static> const>::Type TIter;

    if (beginIdx >= endIdx) // handle base case of empty interval
        return;

    IntervalTreeEntry<TValue> const & entry = tree.entries[centerIdx];

    if (entry._allLeftOf(posBegin)) // posBegin is right of the rightmost point of any interval in this node
        return;

    if (beginIdx < centerIdx) // recurse left
        _findOverlappingWithInterval(tree, beginIdx, centerIdx, beginIdx + (centerIdx - beginIdx) / 2, posBegin,
                                     posEnd, result);
    
    if (entry._overlapsWith(posBegin, posEnd)) // check this node
        appendValue(result, TIter(tree.entries.begin() + centerIdx));

    if (entry._isRightOf(posEnd - 1)) // last interval entry is left of the start of the interval, can't to the right
        return;

    if (centerIdx + 1 < endIdx) // recurse right
        _findOverlappingWithInterval(tree, centerIdx + 1, endIdx, (centerIdx + 1) + (endIdx - (centerIdx + 1)) / 2,
                                     posBegin, posEnd, result);
}

// Guarantee: result is sorted by (beginPos, endPos), by search algorithm.

template <typename TValue, typename TPos, typename TResult>
void findOverlappingWithInterval(IntervalTree<TValue, Static> const & tree,
                                 TPos posBegin,
                                 TPos posEnd,
                                 TResult & result)
{
    typedef IntervalTree<TValue, Static> const                     TIntervalTree;
    typedef typename Position<TIntervalTree>::Type                 TPos2;
    typedef typename Size<TIntervalTree>::Type                     TSize;
    typedef typename Iterator<TIntervalTree const, Standard>::Type TIter;

    _findOverlappingWithInterval(tree, (TSize)0, (TSize)tree.entries.size(), (TSize)(tree.entries.size() / 2),
                                 (TPos2)posBegin, (TPos2)posEnd, result);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_STATIC_INTERVAL_TREE_H_
