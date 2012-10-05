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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_INDEXER_H_
#define SEQAN_EXTRAS_MASAI_INDEXER_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include "store.h"
#include "index.h"

using namespace seqan;

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TGenomeIndex>
struct Indexer
{
    TFragmentStore & store;
    TGenome             genome;
    TGenomeIndex        genomeIndex;

    Indexer(TFragmentStore & store) :
        store(store)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TGenomeIndex, typename TString>
bool loadGenome(Indexer<TGenomeIndex> & indexer, TString const & genomeFile)
{
    if (!loadContigs(indexer.store, genomeFile))
        return 0;

    reserve(indexer.genome, length(indexer.store.contigStore));
    for (unsigned contigId = 0; contigId < length(indexer.store.contigStore); ++contigId)
    {
        shrinkToFit(indexer.store.contigStore[contigId].seq);
        _reverseContig(indexer, indexer.store.contigStore[contigId].seq);
        appendValue(indexer.genome, indexer.store.contigStore[contigId].seq);
    }

    return length(indexer.store.contigStore);
}

template <typename TGenomeIndex, typename TString>
void _reverseContig(Indexer<TGenomeIndex> &, TString &)
{}

template <typename TString>
void _reverseContig(Indexer<TGenomeFM> &, TString & contig)
{
    reverse(contig);
}

template <typename TGenomeIndex, typename TString>
bool loadGenomeIndex(Indexer<TGenomeIndex> & indexer, TString const & genomeIndexFile)
{
    indexer.genomeIndex = TGenomeIndex(indexer.genome);
    return open(indexer.genomeIndex, toCString(genomeIndexFile));
}

template <typename TGenomeIndex, typename TString>
bool dumpIndexedGenome(Indexer<TGenomeIndex> & indexer, TString const & genomeIndexFile)
{
    return save(indexer.genomeIndex, toCString(genomeIndexFile));
}

template <typename TGenomeIndex>
void indexGenome(Indexer<TGenomeIndex> & indexer)
{
    indexer.genomeIndex = TGenomeIndex(indexer.genome);

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TGenomeIndex, TopDown<> >::Type it(indexer.genomeIndex);
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_INDEXER_H_
