// ==========================================================================
//                                   ANISE
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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Haplotype reconstruction / repeat separation globalization using clique
// enumeration.
// ==========================================================================

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_ENUM_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_ENUM_H_

#include <iosfwd>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>
#include <vector>

#include "asm/frag_store.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

struct ReadSeparatorOptions;
struct FeatureReadSet;
class FeatureMap;
class LocalVariationStore;
class CliqueSet;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Function buildAtomicReadSet()
// ----------------------------------------------------------------------------

// Builds a FeatureReadSet where each Read represents the corresponding atomic read from the store.

void buildAtomicReadSet(FeatureReadSet & out,
                        TFragmentStore const & fragStore,
                        LocalVariationStore const & localVarStore);

// ----------------------------------------------------------------------------
// Function buildFeatureReadSet()
// ----------------------------------------------------------------------------

// Build an initial read set from a FragmentStore and a LocalVariationStore into out.  This set is appropriate for
// clique enumeration, meaning that read pairs on the same contig are already joined into one Read.

void buildFeatureReadSet(FeatureReadSet & out,
                         TFragmentStore const & fragStore,
                         LocalVariationStore const & localVarStore);

// ----------------------------------------------------------------------------
// Function mergeToSuperReads()
// ----------------------------------------------------------------------------

// Merge the reads in each clique to give rise to superreads.
//
// The reads in in must be sorted by id.

void mergeToSuperReads(FeatureReadSet & out,
                       CliqueSet const & cliques,
                       FeatureReadSet const & in,
                       bool logging = false);

// ----------------------------------------------------------------------------
// Function performCliqueEnumeration()
// ----------------------------------------------------------------------------

// Process FeatureReadSet into CliqueSet.
//
// Returns true if reads in out can now be joined into super reads and iteration can continue.

bool performCliqueEnumeration(CliqueSet & out,
                              FeatureReadSet const & in,
                              ReadSeparatorOptions const & options);

// ----------------------------------------------------------------------------
// Function mergeNonConflictingReadsOnContig()
// ----------------------------------------------------------------------------

// For each contig, pick sets of non-conflicting ("independent sets") reads to merge to larger reads.

void mergeNonConflictingReadsOnContig(FeatureReadSet & out,
                                      FeatureReadSet const & in,
                                      FeatureReadSet const & initialReadSet,
                                      ReadSeparatorOptions const & options);

// ----------------------------------------------------------------------------
// Function split()
// ----------------------------------------------------------------------------

// Split the given TFragmentStore and BamAlignmentRecord objects using the FeatureReadSet.  Split reads will be created with
// the given prefix.

void split(TFragmentStore & outStore,
           FeatureReadSet const & readSet,
           TFragmentStore const & inStore,
           ReadSeparatorOptions const & options);

// TODO(holtgrew): Check whether the assumption is correct that we get a "clean" separation over contig borders with pairs using a SEQAN_CHECK().

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_CLIQUE_ENUM_H_
