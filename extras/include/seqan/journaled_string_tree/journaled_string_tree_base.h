// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Tag ContextPositionLeft
// ----------------------------------------------------------------------------

struct ContextPositionLeft_;
typedef Tag<ContextPositionLeft_> ContextPositionLeft;


// ----------------------------------------------------------------------------
// Tag ContextPositionRight
// ----------------------------------------------------------------------------

struct ContextPositionRight_;
typedef Tag<ContextPositionRight_> ContextPositionRight;

// ----------------------------------------------------------------------------
// Struct JstTraverserSpec
// ----------------------------------------------------------------------------

template <typename TContextPosition = ContextPositionLeft, typename RequireFullContext = True>
struct JstTraverserSpec;

// ----------------------------------------------------------------------------
// Tag StringTreeDefault
// ----------------------------------------------------------------------------

struct DefaultStringTree_;
typedef Tag<DefaultStringTree_> StringTreeDefault;

// ----------------------------------------------------------------------------
// Tag StringTreeSparse
// ----------------------------------------------------------------------------

struct SparseStringTree_;
typedef Tag<SparseStringTree_> StringTreeSparse;

// ----------------------------------------------------------------------------
// Class JournaledStringTree
// ----------------------------------------------------------------------------

template <typename TDeltaMap, typename TSpec = StringTreeDefault>
class JournaledStringTree;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction GetState
// ----------------------------------------------------------------------------

template <typename T>
struct GetState
{
    typedef Nothing Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetJstTraverser
// ----------------------------------------------------------------------------

template <typename T>
struct GetJstTraverser;

// ----------------------------------------------------------------------------
// Metafunction ContextIteratorPosition
// ----------------------------------------------------------------------------

template <typename T>
struct ContextIteratorPosition
{
    typedef ContextPositionLeft Type;
};

// ----------------------------------------------------------------------------
// Metafunction RequireFullContext
// ----------------------------------------------------------------------------

template <typename T>
struct RequireFullContext : True{};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_
