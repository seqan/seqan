// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2014, Knut Reinert, FU Berlin
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
// Tag EnhancedDeltaMap
// ----------------------------------------------------------------------------

// TODO(rmaerker): Docu!
struct EnhancedDeltaMap_;
typedef Tag<EnhancedDeltaMap_> EnhancedDeltaMap;

// ----------------------------------------------------------------------------
// Class JstConfig
// ----------------------------------------------------------------------------

template <typename TReference>
struct DefaultJstConfig
{
    typedef typename Value<TReference>::Type TSnpValue;  // Value type of the SNPs.
    typedef String<TSnpValue>                TInsValue;  // Value type of insertions.
    typedef typename Size<TReference>::Type  TDelValue;  // Value type of deletions.
    typedef Pair<TDelValue, TInsValue>       TSVValue;   // Value type of structural variants (combination of deletion and insertion).
};

// ----------------------------------------------------------------------------
// Class JournaledStringTree
// ----------------------------------------------------------------------------

template <typename TReference, typename TSpec = Default, typename TConfig = DefaultJstConfig<TReference> >cp 
class JournaledStringTree;

// TODO(rmaerker): Should be in a different commit. Likely to change.
//// ----------------------------------------------------------------------------
//// Tag ContextPositionLeft
//// ----------------------------------------------------------------------------
//
//struct ContextPositionLeft_;
//typedef Tag<ContextPositionLeft_> ContextPositionLeft;
//
//
//// ----------------------------------------------------------------------------
//// Tag ContextPositionRight
//// ----------------------------------------------------------------------------
//
//struct ContextPositionRight_;
//typedef Tag<ContextPositionRight_> ContextPositionRight;
//
//// TODO(rmaerker): Need later to distinguish between extension (local search) and global search.
////// ----------------------------------------------------------------------------
////// Tag TraversePrefix
////// ----------------------------------------------------------------------------
////
////struct TraversePrefix_;
////typedef Tag<TraversePrefix_> TraversePrefix;
////
////// ----------------------------------------------------------------------------
////// Tag TraverseInfix
////// ----------------------------------------------------------------------------
////
////struct TraverseInfix_;
////typedef Tag<TraverseInfix_> TraverseInfix;
//
//// ----------------------------------------------------------------------------
//// Struct JstTraverserConfig
//// ----------------------------------------------------------------------------
//
//template <typename TContextPosition = ContextPositionLeft, typename TRequireFullContext = True>
//struct JstTraverserConfig;
//
// ============================================================================
// Metafunctions
// ============================================================================

//// ----------------------------------------------------------------------------
//// Metafunction GetState
//// ----------------------------------------------------------------------------
//
//template <typename T>
//struct GetState
//{
//    typedef Nothing Type;
//};
//
//// ----------------------------------------------------------------------------
//// Metafunction GetJstTraverser
//// ----------------------------------------------------------------------------
//
//template <typename T>
//struct GetJstTraverser;
//
//// ----------------------------------------------------------------------------
//// Metafunction ContextIteratorPosition
//// ----------------------------------------------------------------------------
//
//template <typename T>
//struct ContextIteratorPosition
//{
//    typedef ContextPositionLeft Type;
//};
//
//// ----------------------------------------------------------------------------
//// Metafunction RequireFullContext
//// ----------------------------------------------------------------------------
//
//template <typename T>
//struct RequireFullContext : True{};
//
//// ----------------------------------------------------------------------------
//// Metafunction ConfigureJstTraversal
//// ----------------------------------------------------------------------------
//
//template <typename TAlgo, typename TType, typename TDirection>
//struct ConfigureJstTraversal;

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_JOURNALED_STRING_TREE_BASE_H_
