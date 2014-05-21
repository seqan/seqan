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
// Facade header for the journaled string tree.
// ==========================================================================

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_H_

// ===========================================================================
// Prerequisites.
// ===========================================================================

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/sequence_journaled.h>
#include <seqan/journaled_set.h>
#include <seqan/parallel.h>

// ===========================================================================
// Basics.
// ===========================================================================

#include <seqan/journaled_string_tree/delta_store.h>
#include <seqan/journaled_string_tree/delta_coverage_store.h>
#include <seqan/journaled_string_tree/delta_map.h>
#include <seqan/journaled_string_tree/delta_map_iterator.h>

// ===========================================================================
// Journaled String Tree.
// ===========================================================================

#include <seqan/journaled_string_tree/journaled_string_tree_base.h>
#include <seqan/journaled_string_tree/journaled_string_tree_util.h>
#include <seqan/journaled_string_tree/journaled_string_tree_adapter.h>
#include <seqan/journaled_string_tree/journaled_string_tree_impl.h>

// ===========================================================================
// Journaled String Tree Traversal.
// ===========================================================================

#include <seqan/journaled_string_tree/journaled_string_tree_traversal_concept.h>
#include <seqan/journaled_string_tree/journaled_string_tree_traversal_branch_stack.h>
#include <seqan/journaled_string_tree/journaled_string_tree_traversal.h>

#endif // EXTRAS_INCLUDE_SEQAN_JOURNALED_STRING_TREE_H_
