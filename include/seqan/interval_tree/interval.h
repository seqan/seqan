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

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_H_

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
// Class Interval
// ----------------------------------------------------------------------------

/*!
 * @class IntervalConcept
 * @brief Concept for half-open intervals.
 */
// beginPos()
// endPos()
// cargo()

/*!
 * @class Interval
 * @headerfile <seqan/interval_tree.h>
 * @implements IntervalConcept
 * @brief Half-open interval.
 *
 * @signature template <typename TCargo, typename TPos>
 *            class Interval;
 *
 * @tparam TCargo Type of the cargo  to store.
 * @tparam TPos   Type to use for positions.
 */
template <typename TCargo, typename TPos>
class Interval
{
public:
    TCargo cargo;
    TPos beginPos;
    TPos endPos;

    Interval() : cargo(), beginPos(), endPos()
    {}

    Interval(TCargo const & cargo, TPos beginPos, TPos endPos) :
            cargo(cargo), beginPos(beginPos), endPos(endPos)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TCargo, typename TPos>
struct Position<Interval<TCargo, TPos> >
{
    typedef TPos Type;
};

template <typename TCargo, typename TPos>
struct Position<Interval<TCargo, TPos> const>
{
    typedef TPos Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TCargo, typename TPos>
struct Value<Interval<TCargo, TPos> >
{
    typedef TCargo Type;
};

template <typename TCargo, typename TPos>
struct Value<Interval<TCargo, TPos> const>
{
    typedef TCargo Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TReference, typename TPos>
struct Reference<Interval<TReference, TPos> >
{
    typedef TReference & Type;
};

template <typename TReference, typename TPos>
struct Reference<Interval<TReference, TPos> const>
{
    typedef TReference const & Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TCargo, typename TPos>
TPos beginPos(Interval<TCargo, TPos> const & interval)
{
    return interval.beginPos;
}

template <typename TCargo, typename TPos>
TPos endPos(Interval<TCargo, TPos> const & interval)
{
    return interval.endPos;
}

template <typename TCargo, typename TPos>
typename Reference<Interval<TCargo, TPos> const>::Type cargo(Interval<TCargo, TPos> const & interval)
{
    return interval.cargo;
}

template <typename TCargo, typename TPos>
typename Reference<Interval<TCargo, TPos> >::Type cargo(Interval<TCargo, TPos> & interval)
{
    return interval.cargo;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_INTERVAL_H_
