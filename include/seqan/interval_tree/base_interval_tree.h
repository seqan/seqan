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

#ifndef INCLUDE_SEQAN_INTERVAL_TREE_BASE_INTERVAL_TREE_H_
#define INCLUDE_SEQAN_INTERVAL_TREE_BASE_INTERVAL_TREE_H_

#include <seqan/basic.h>

#include <seqan/misc/interval_tree.h>

namespace seqan
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class IntervalTree
// ----------------------------------------------------------------------------

/*!
 * @class IntervalTree
 * @implements ContainerConcept
 * @headerfile <seqan/interval_tree.h>
 * @brief Data structure for fast range queries on a set of intervals.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class IntervalTree;
 *
 * @tparam TValue Value to store in the tree.
 * @tparam TSpec  Specializing tag.
 */
template <typename TValue, typename TSpec>
class IntervalTree;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

/**
 * @mfn IntervalTree#Position
 * @brief The position type of the <tt>IntervalTree</tt>, forwards to <tt>Position&lt;TValue&gt;</tt>.
 */
template <typename TValue, typename TSpec>
struct Position<IntervalTree<TValue, TSpec> >
{
    typedef typename Position<TValue>::Type Type;
};

template <typename TValue, typename TSpec>
struct Position<IntervalTree<TValue, TSpec> const>
{
    typedef typename Position<TValue const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

/**
 * @mfn IntervalTree#Value
 * @brief The value type of the <tt>IntervalTree</tt>.
 */
// TODO(holtgrew): Redefined in legacy module.
// template <typename TValue, typename TSpec>
// struct Value<IntervalTree<TValue, TSpec> >
// {
//     typedef TValue Type;
// };

// template <typename TValue, typename TSpec>
// struct Value<IntervalTree<TValue, TSpec> const>
// {
//     typedef TValue Type;
// };

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // INCLUDE_SEQAN_INTERVAL_TREE_BASE_INTERVAL_TREE_H_
