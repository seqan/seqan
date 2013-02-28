//// ==========================================================================
////                 SeqAn - The Library for Sequence Analysis
//// ==========================================================================
//// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
//// All rights reserved.
////
//// Redistribution and use in source and binary forms, with or without
//// modification, are permitted provided that the following conditions are met:
////
////     * Redistributions of source code must retain the above copyright
////       notice, this list of conditions and the following disclaimer.
////     * Redistributions in binary form must reproduce the above copyright
////       notice, this list of conditions and the following disclaimer in the
////       documentation and/or other materials provided with the distribution.
////     * Neither the name of Knut Reinert or the FU Berlin nor the names of
////       its contributors may be used to endorse or promote products derived
////       from this software without specific prior written permission.
////
//// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
//// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
//// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//// DAMAGE.
////
//// ==========================================================================
//// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
//// ==========================================================================
//// TODO(esiragusa): Add file description.
//// ==========================================================================
//
//#ifndef SEQAN_EXTRAS_MASAI_INDEXER_H_
//#define SEQAN_EXTRAS_MASAI_INDEXER_H_
//
//#include <seqan/basic.h>
//#include <seqan/sequence.h>
//#include <seqan/file.h>
//
//#include "index.h"
//
//using namespace seqan;
//
//// ============================================================================
//// Forwards
//// ============================================================================
//
//// ============================================================================
//// Tags, Classes, Enums
//// ============================================================================
//
//// ----------------------------------------------------------------------------
//// Class GenomeIndex
//// ----------------------------------------------------------------------------
//
//template <typename TGenome, typename TIndex, typename TSpec = void>
//struct GenomeIndex
//{
//    TGenome             & genome;
//    TIndex              index;
//
//    GenomeIndex(TGenome & genome) :
//        genome(genome)
//    {}
//};
//
//// ----------------------------------------------------------------------------
//// Class ReadsIndex
//// ----------------------------------------------------------------------------
//
//template <typename TReads, typename TIndex, typename TSpec = void>
//struct ReadsIndex
//{
//    TReads              const & reads;
//    TIndex              index;
//
//    ReadsIndex(TReads const & reads) :
//        reads(reads)
//    {}
//};
//
//// ============================================================================
//// Metafunctions
//// ============================================================================
//
//// ============================================================================
//// Functions
//// ============================================================================
//
//
//// ----------------------------------------------------------------------------
//// Function load()                                                [GenomeIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TGenome, typename TIndex, typename TSpec, typename TString>
//bool load(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex, TString const & genomeIndexFile)
//{
//    genomeIndex.index = TIndex(genomeIndex.genome.contigs);
//
//    return open(genomeIndex.index, toCString(genomeIndexFile));
//}
//
//// ----------------------------------------------------------------------------
//// Function build()                                               [GenomeIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TGenome, typename TIndex, typename TSpec>
//void build(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex)
//{
//    genomeIndex.index = TIndex(genomeIndex.genome.contigs);
//
//    // Iterator instantiation calls automatic index construction.
//    typename Iterator<TIndex, TopDown<> >::Type it(genomeIndex.index);
//}
//
//template <typename TGenome, typename TSpec>
//void build(GenomeIndex<TGenome, TGenomeFM, TSpec> & genomeIndex)
//{
//    // IndexFM is built on the reversed genome.
//    reverse(genomeIndex.genome);
//
//    genomeIndex.index = TGenomeFM(genomeIndex.genome.contigs);
//
//    // Iterator instantiation calls automatic index construction.
//    typename Iterator<TGenomeFM, TopDown<> >::Type it(genomeIndex.index);
//
//    reverse(genomeIndex.genome);
//}
//
//// ----------------------------------------------------------------------------
//// Function dump()                                                [GenomeIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TGenome, typename TIndex, typename TSpec, typename TString>
//bool dump(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex, TString const & genomeIndexFile)
//{
//    return save(genomeIndex.index, toCString(genomeIndexFile));
//}
//
//// ----------------------------------------------------------------------------
//// Function clear()                                               [GenomeIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TGenome, typename TSpec>
//void clear(GenomeIndex<TGenome, TGenomeFM, TSpec> & genomeIndex)
//{
//    clear(genomeIndex.index);
//}
//// ----------------------------------------------------------------------------
//// Function build()                                                [ReadsIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TReads, typename TIndex, typename TSpec>
//void build(ReadsIndex<TReads, TIndex, TSpec> & readsIndex)
//{
//}
//
//// ----------------------------------------------------------------------------
//// Function clear()                                                [ReadsIndex]
//// ----------------------------------------------------------------------------
//
//template <typename TReads, typename TIndex, typename TSpec>
//void clear(ReadsIndex<TReads, TIndex, TSpec> & readsIndex)
//{
//}
//
//#endif  // #ifndef SEQAN_EXTRAS_MASAI_INDEXER_H_
