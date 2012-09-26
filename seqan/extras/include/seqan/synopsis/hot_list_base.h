// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef SEQAN_SYNOPSIS_HOT_LIST_BASE_H_
#define SEQAN_SYNOPSIS_HOT_LIST_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.HotList
..cat:Synopsis Data Structures
..summary:Data structure for (approximately) keeping a list of the k most frequent items in a data stream.
..signature:HotList<TValue, TSpec>
..param.TValue:Type of the items to keep in hot list.
..param.TSpec:Tag used for the specialization.
..include:seqan/synopsis.h
 */

template <typename TValue, typename TSpec, typename TCallback = Nothing>
class HotList;

// Helper comparator functor, used for sorting items in getItems().
template <typename T1, typename T2>
struct CmpTripleI2DescI3DescI1Asc_
        : std::binary_function<Triple<T1, T2, T2>, Triple<T1, T2, T2>, bool>
{
    typedef Triple<T1, T2, T2> TTriple;
    bool operator()(TTriple const & a, TTriple const & b) const
    {
        return (a.i2 > b.i2) || (a.i2 == b.i2 && a.i3 > b.i3) || (a.i2 == b.i2 && a.i3 == b.i3 && a.i1 < b.i1);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T:Class.HotList
///.Metafunction.Value.class:Class.HotList

template <typename TValue, typename TSpec, typename TCallback>
struct Value<HotList<TValue, TSpec, TCallback> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec, typename TCallback>
struct Value<HotList<TValue, TSpec, TCallback> const>
{
    typedef TValue const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T:Class.HotList
///.Metafunction.Size.class:Class.HotList

template <typename TValue, typename TSpec, typename TCallback>
struct Size<HotList<TValue, TSpec, TCallback> >
{
    typedef unsigned Type;
};

template <typename TValue, typename TSpec, typename TCallback>
struct Size<HotList<TValue, TSpec, TCallback> const>
{
    typedef unsigned Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function hotListItemAdded()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Add a HotListCallback Concept?

// A function of this name and interface (except that callbackContext is the
// TCallback from the user's HotList specialization) is called whenever an
// item is added to a hot list.  value is the value that is removed and the
// lower and upper bounds for the count are given.
template <typename TValue, typename TCount>
inline void
hotListItemAdded(Nothing const & /*callbackContext*/,
                 TValue const & /*value*/,
                 TCount const & /*lowerBound*/,
                 TCount const & /*upperBound*/)
{}

// ----------------------------------------------------------------------------
// Function hotListItemRemoved()
// ----------------------------------------------------------------------------

// A function of this name and interface (except that callbackContext is the
// TCallback from the user's HotList specialization) is called whenever an
// item is removed from a hot list.  value is the value that is removed and
// the lower and upper bounds for the count are given.
template <typename TValue, typename TCount>
inline void
hotListItemRemoved(Nothing const & /*callbackContext*/,
                   TValue const & /*value*/,
                   TCount const & /*lowerBound*/,
                   TCount const & /*upperBound*/)
{}


// ----------------------------------------------------------------------------
// Function _callbackContext()
// ----------------------------------------------------------------------------

// Return the callback context.
template <typename TValue, typename TSpec, typename TCallback>
inline TCallback &
_callbackContext(HotList<TValue, TSpec, TCallback> & hotList)
{
    return value(hotList._callbackContext);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

///.Function.clear.param.object.type:Class.HotList
///.Function.clear.class:Class.HotList

// ----------------------------------------------------------------------------
// Function registerItem()
// ----------------------------------------------------------------------------

/**
.Function.registerItem:
..class:Class.HotList
..cat:Synopsis Data Structures
..summary:Register an item as seen with a @Class.HotList@.
..signature:registerItem(hotList, value)
..param.hotList:The @Class.HotList@ to register the value with.
...type:Class.HotList
..param.value:The value to observe as registered.
..include:seqan/synopsis.h
 */

// ----------------------------------------------------------------------------
// Function removeItem()
// ----------------------------------------------------------------------------

/**
.Function.removeItem
..class:Class.HotList
..cat:Synopsis Data Structures
..summary:Manually remove an item from a @Class.HotList@.
..signature:removeItem(hotList, value)
..param.hotList:The @Class.HotList@ to remove the value from.
...type:Class.HotList
..param.value:The value to remove from @Class.HotList@.
..include:seqan/synopsis.h
 */

// ----------------------------------------------------------------------------
// Function getItems()
// ----------------------------------------------------------------------------

/**
.Function.getItems
..class:Class.HotList
..cat:Synopsis Data Structures
..summary:Get the frequent items.
..signature:getItems(items, hotList)
..param.items:List of value/frequency estimation @Class.Triple|Triples@.
...type:nolink:$String<Triple<TValue, TSize, TSize> >$
..param.hotList:The @Class.HotList@ to remove the value from.
...type:Class.HotList
..remarks:The result is ordered descendingly by estimated frequency, ties are broken by value.
..include:seqan/synopsis.h
 */

}  // namespace seqan

#endif  // SEQAN_SYNOPSIS_HOT_LIST_BASE_H_
