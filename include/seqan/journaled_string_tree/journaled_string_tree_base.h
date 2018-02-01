// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// Basic defintions and forwards used globally for this module.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Traversal Tags
// ----------------------------------------------------------------------------

struct ForwardTraversal_;
typedef Tag<ForwardTraversal_> ForwardTraversal;

struct BidirectionalTraversal_;
typedef Tag<BidirectionalTraversal_> BidirectionalTraversal;

// ----------------------------------------------------------------------------
// Tag JSTree
// ----------------------------------------------------------------------------

template <typename TSpec = ForwardTraversal>
struct JSTree;

// ----------------------------------------------------------------------------
// Tag JstBufferMember
// ----------------------------------------------------------------------------

struct JstBufferMember_;
typedef Tag<JstBufferMember_> JstBufferMember;

// ----------------------------------------------------------------------------
// Tag JstSourceMember
// ----------------------------------------------------------------------------

struct JstSourceMember_;
typedef Tag<JstSourceMember_> JstSourceMember;

// ----------------------------------------------------------------------------
// Tag JstDeltaMapMember
// ----------------------------------------------------------------------------

struct JstDeltaMapMember_;
typedef Tag<JstDeltaMapMember_> JstDeltaMapMember;

// ----------------------------------------------------------------------------
// Class DefaultJstConfig
// ----------------------------------------------------------------------------

/*!
 * @class DefaultJstConfig
 * @headerfile <seqan/journaled_string_tree.h>
 * @brief The default Journaled-String-Tree configuration object.
 *
 * @signature template <typename TSequence>
 *            class DefaultJstConfig<TDeltaStore>
 *
 * @tparam TSequence Type of the underlying base sequence.
 *
 * Defines the follwoing types used to configure the @link DeltaMap @endlink:
 * @code{.cpp}
 * typedef typename Size<TSequence>::Type  TDeltaPos;
 * typedef typename Value<TSequence>::Type TSnpValue;
 * typedef String<TSnpValue>               TInsValue;
 * typedef typename Size<TSequence>::Type  TDelValue;
 * typedef Pair<TDelValue, TInsValue>      TSVValue;
 * @endcode
 */
template <typename TSequence>
struct DefaultJstConfig
{
    typedef typename Size<TSequence>::Type  TDeltaPos;  // Position type for delta entry.
    typedef typename Value<TSequence>::Type TSnpValue;  // Value type of the SNPs.
    typedef String<TSnpValue>               TInsValue;  // Value type of insertions.
    typedef typename Size<TSequence>::Type  TDelValue;  // Value type of deletions.
    typedef Pair<TDelValue, TInsValue>      TSVValue;   // Value type of structural variants (combination of deletion and insertion).
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Traverser
// ----------------------------------------------------------------------------

template <typename TContainer, typename TObserverList = ObserverList<> >
struct Traverser;

// ----------------------------------------------------------------------------
// Metafunction ContextIterator
// ----------------------------------------------------------------------------

template <typename TObj>
struct ContextIterator;

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_
