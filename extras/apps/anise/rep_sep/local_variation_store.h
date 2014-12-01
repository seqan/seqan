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

#ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_LOCAL_VARIATION_STORE_H_
#define SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_LOCAL_VARIATION_STORE_H_

#include <seqan/basic.h>

#include "asm/frag_store.h"

namespace rep_sep {

// ============================================================================
// Forwards
// ============================================================================

struct ReadSeparatorOptions;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class LocalVariationStore
{
public:
    // The alphabet used in the reads themselves.
    typedef seqan::Dna5 TAlphabet;
    // The alphabet used for the consensus, includes gaps.
    typedef seqan::ModifiedAlphabet<TAlphabet, seqan::ModExpand<'-'> > TConsensusAlphabet;
    // Type of the profile, consensus, and for the occurences in the columns.
    typedef seqan::ProfileChar<TAlphabet> TProfileChar;
    typedef seqan::String<TProfileChar> TProfile;
    typedef seqan::String<TConsensusAlphabet> TConsensus;
    typedef seqan::Triple<unsigned, TConsensusAlphabet, char> TColumnEntry;  // (read id, value, phred qual)

    // Contigs and positions of the columns (contigID, pos);
    seqan::String<std::pair<int, int> > positions;
    // The consensus character for the columns.
    TConsensus consensus;
    // The profile in the selected columns, better resolution than consensus but worse than columns.
    TProfile profile;
    // The entries for each column.
    seqan::StringSet<seqan::String<TColumnEntry> > columns;
    // Reads covered by the column.
    seqan::StringSet<seqan::String<unsigned> > coveringReads;
    // Mapping between the compressed variation store's read ids to the uncompressed ones.
    std::map<unsigned, seqan::String<unsigned> > compressionMap;

    LocalVariationStore() = default;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

void clear(LocalVariationStore & store);

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

inline unsigned length(LocalVariationStore const & store)
{
    return length(store.positions);
}

// ----------------------------------------------------------------------------
// Function fill()
// ----------------------------------------------------------------------------

// Fill local variation store with deviating columns.

void fill(LocalVariationStore & varStore,
          TFragmentStore const & fragStore,
          ReadSeparatorOptions const & options);

// ----------------------------------------------------------------------------
// Function filter()
// ----------------------------------------------------------------------------

// Filter the given local variation store entries to a list of (sorted) columns.

void filter(LocalVariationStore & varStore,
            seqan::String<unsigned> const & columnIds);

}  // namespace rep_sep

#endif  // #ifndef SANDBOX_HERBARIUM_APPS_ANISE_REP_SEP_LOCAL_VARIATION_STORE_H_
