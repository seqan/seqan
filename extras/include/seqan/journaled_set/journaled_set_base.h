// ==========================================================================
//                               journaled_set
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

#ifndef EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_
#define EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Spec JournaledSet
// ----------------------------------------------------------------------------

/**
.Spec.Journaled Set:
..summary:A string set storing the string as members. It can store a global reference sequence to which
all members can be journaled if they are of type @Spec.Journaled String@.
..cat:Sequences
..signature:StringSet<TString, Owner<JournaledSet> >
..param.TString:The string type.
...type:Class.String
...type:Spec.Journaled String
..remarks:The strings are internally stored in a $String<TString>$ object and the character position type is a
@Class.Pair@ $(seqNo,seqOfs)$ where seqNo identifies the string within the stringset and seqOfs identifies the position within this string.
..remarks:The global reference is of type $Host<TString>$. Only strings of type @Spec.Journaled String@ or $Host<$@Spec.Journaled String@>$
can be used for the advanced functionality supported by this string set.
..include:seqan/journal_set.h
 */

struct JournaledSet_;
typedef Tag<JournaledSet_> JournaledSet;


// ----------------------------------------------------------------------------
// Class JournalTraceBuffer
// ----------------------------------------------------------------------------

template <typename TString>
class JournalTraceBuffer;

/**
.Tag.Join Strategies
..cat:Alignments
..summary:Tags used for selecting journaling strategies when joining a @Spec.Journaled String@ to a global reference sequence.
..tag
...JournaledManhatten:Constructs @Spec.Journaled String@ based on Manhatten distance.
....remarks:This strategy is very fast on the cost of memory.
...JournaledCompact:Computes an optimal alignment to construct @Spec.Journaled String@.
....remarks:This strategy is slow but depending on the scoring function minimizes the memory requirements for the computed @Spec.Journaled String@ (default).
..see:Function.join
..include:seqan/align.h
*/

// ----------------------------------------------------------------------------
// Tag JournaledManhatten
// ----------------------------------------------------------------------------

struct JournaledManhatten_;
typedef Tag<JournaledManhatten_> JournaledManhatten;

// ----------------------------------------------------------------------------
// Tag JournaledCompact
// ----------------------------------------------------------------------------

struct JournaledCompact_;
typedef Tag<JournaledCompact_> JournaledCompact;


// ----------------------------------------------------------------------------
// Spec GlobalAlign
// ----------------------------------------------------------------------------

/**
.Spec.GlobalAlign:
..summary:Selects a global alignment method to join a @Spec.Journaled String@ to a global reference sequence.
..cat:Sequences
..general:Class.JoinConfig
..signature:GlobalAlign<TStrategyTag>
..param.TStrategyTag:The strategy used to compute the journal.
...type:Tag.Join Strategies.tag.JournaledManhatten
...type:Tag.Join Strategies.tag.GlobalChain
..default:Tag.Join Strategies.tag.JournaledManhatten
..remarks:If @Tag.Join Strategies.tag.JournaledManhatten@ is selected, then the resulting @Spec.Journaled String@ consists
of one insertion node covering the complete joined sequence.
..see:Spec.GlobalChain
..include:seqan/journal_set.h
 */

template <typename TSpec = JournaledManhatten>
struct GlobalAlign{};

// ----------------------------------------------------------------------------
// Spec GlobalChain
// ----------------------------------------------------------------------------

/**
.Spec.GlobalChain:
..summary:Selects an anchor based method to join a @Spec.Journaled String@ to a global reference sequence.
..cat:Sequences
..general:Class.JoinConfig
..signature:GlobalChain<TStrategyTag>
..param.TStrategyTag:The strategy used to compute the journal.
...type:Tag.Join Strategies.tag.JournaledManhatten
...type:Tag.Join Strategies.tag.JournaledCompact
..default:Tag.Join Strategies.tag.JournaledManhatten
..remarks:The @Tag.Join Strategies.tag.JournaledManhatten@ strategy fills the gaps between the anchors with a single insertion node,
while the corresponding part of the reference sequence is deleted.
..see:Spec.GlobalAlign
..include:seqan/journal_set.h
 */

template <typename TSpec = JournaledManhatten>
struct GlobalChain{};

// ----------------------------------------------------------------------------
// Class JoinConfig
// ----------------------------------------------------------------------------

/**
.Class.JoinConfig:
..summary:Specifies the strategy and all necessary parameters used to journal a sequence to a reference sequence.
..cat:Sequences
..signature:JoinConfig<TMethod>
..param.TMethod:The method type.
...type:Spec.GlobalAlign
...type:Spec.GlobalChain
..default:Spec.JournaledManhatten
..remarks:SeqAn offers two general methods to compute the journal. The first one uses a @Function.globalAlignment@ function
and the second one uses an anchor based approach.
..include:seqan/journal_set.h
 */

template <typename TSpec = GlobalAlign<JournaledManhatten> >
struct JoinConfig{};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


}  // namespace seqan

#endif  // EXTRAS_INCLUDE_SEQAN_JOURNALED_SET_JOURNALED_SET_BASE_H_
