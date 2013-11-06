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
// This file contains index type definitions.
// ==========================================================================

#ifndef SEQAN_EXTRAS_MASAI_INDEX_H_
#define SEQAN_EXTRAS_MASAI_INDEX_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

#include "index/index_qgram_stretched.h"
#include "index/find_backtracking_stretched.h"

#include "store.h"

using namespace seqan;

// ============================================================================
// Types
// ============================================================================

// ----------------------------------------------------------------------------
// Shape Type
// ----------------------------------------------------------------------------

typedef UngappedShape<10>                               TShape;

// ============================================================================
// Contigs Index Fibres
// ============================================================================

// ----------------------------------------------------------------------------
// Contigs Suffix Array Value Type
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TContigs>
{
    typedef Pair<unsigned char, unsigned int, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// Contigs Enhanced Suffix Array Fibres
// ----------------------------------------------------------------------------

typedef Index<TContigs, IndexEsa<> >                     TGenomeEsa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP)
typedef DefaultIndexStringSpec<TGenomeEsa>::Type        TGenomeEsaStringSpec;
#else
typedef MMap<>                                          TGenomeEsaStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeEsa, FibreSA>
{
    typedef String<SAValue<TContigs>::Type, TGenomeEsaStringSpec>    Type;
};

template <>
struct Fibre<TGenomeEsa, FibreLcp>
{
    typedef String<unsigned int, TGenomeEsaStringSpec>   Type;
};

template <>
struct Fibre<TGenomeEsa, FibreChildtab>
{
    typedef String<unsigned int, TGenomeEsaStringSpec>   Type;
};
}

// ----------------------------------------------------------------------------
// Contigs Suffix Array Fibres
// ----------------------------------------------------------------------------

typedef Index<TContigs, IndexSa<> >                      TGenomeSa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP)
typedef DefaultIndexStringSpec<TGenomeSa>::Type         TGenomeSaStringSpec;
#else
typedef MMap<>                                          TGenomeSaStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeSa, FibreSA>
{
    typedef String<SAValue<TContigs>::Type, TGenomeSaStringSpec>     Type;
};
}

// ----------------------------------------------------------------------------
// Contigs QGram Index with Bucket Refinement Fibres
// ----------------------------------------------------------------------------

typedef IndexQGram<TShape>                              TQGramBaseIndex;
typedef IndexQGram<TShape, BucketRefinement>            TQGramBucketRefinementIndex;

typedef Index<TContigs, TQGramBaseIndex>                 TGenomeBaseQGram;
typedef Index<TContigs, IndexSa<InfixSegment> >          TGenomeInfixSa;
typedef Index<TContigs, TQGramBucketRefinementIndex>     TGenomeQGram;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP)
typedef DefaultIndexStringSpec<TGenomeBaseQGram>::Type  TGenomeQGramStringSpec;
#else
typedef MMap<>                                          TGenomeQGramStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeBaseQGram, FibreDir>
{
    typedef String<unsigned int, TGenomeQGramStringSpec>   Type;
};
}

// ----------------------------------------------------------------------------
// Contigs FM Index Fibres
// ----------------------------------------------------------------------------

typedef Index<TContigs, FMIndex<WT<>, CompressText> >    TGenomeFM;
//typedef Index<TContigs, FMIndex<BMS<>, CompressText> >    TGenomeFM;

// ----------------------------------------------------------------------------
// Contigs Uncompressed FM Index Fibres
// ----------------------------------------------------------------------------

//#if defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
//typedef DefaultIndexStringSpec<TGenomeFM>::Type         TGenomeFMStringSpec;
//#else
//typedef External<ExternalConfigLarge<File<>,8192,2> >   TGenomeFMStringSpec;
//#endif
//
//namespace seqan {
//template <>
//struct Fibre<TGenomeFM, FibreSA> {
//    typedef String<SAValue<TGenomeFM>::Type, TGenomeFMStringSpec> Type;
//};
//}

// ============================================================================
// Reads Index Fibres
// ============================================================================

// ----------------------------------------------------------------------------
// Reads Suffix Array Value Type
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TReadSeqStore>
{
    typedef Pair<unsigned int, unsigned short, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// Reads Wotd Fibres
// ----------------------------------------------------------------------------

typedef Index<TReadSeqStore, IndexWotd<> >              TReadsWotd;

//namespace seqan {
//template <>
//struct Fibre<TReadsWotd, FibreDir>
//{
//    typedef String<unsigned int, DefaultIndexStringSpec<TReadsWotd>::Type>   Type;
//};
//}

// ----------------------------------------------------------------------------
// Reads QGram Index with Bucket Refinement Fibres
// ----------------------------------------------------------------------------

//typedef Index<TReadSeqStore, TQGramBaseIndex>                 TReadsBaseQGram;
//typedef Index<TReadSeqStore, TQGramBucketRefinementIndex>     TReadsQGram;
//
//namespace seqan
//{
//template <>
//struct Fibre<TReadsBaseQGram, FibreDir>
//{
//    typedef String<unsigned int, DefaultIndexStringSpec<TReadsBaseQGram>::Type >   Type;
//};
//}

// ----------------------------------------------------------------------------
// Reads QGram Index Fibres
// ----------------------------------------------------------------------------

typedef Index<TReadSeqStore, TQGramBaseIndex>               TReadsQGram;

namespace seqan {
template <>
struct Fibre<TReadsQGram, FibreDir>
{
    typedef String<unsigned int, DefaultIndexStringSpec<TReadsQGram>::Type>    Type;
};
}

// ----------------------------------------------------------------------------
// Reads QGram Index with Wotd Subtree Fibres
// ----------------------------------------------------------------------------

//typedef Index<TReadSeqStore, IndexWotd< Subtree<> > >   TReadsWotdSubtree;
//
//namespace seqan
//{
//    template <>
//    struct Fibre<TReadsWotdSubtree, FibreSA>
//    {
//        typedef Segment< Fibre<TReadsQGram, FibreSA>::Type, InfixSegment>  Type;
//    };
//}

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class GenomeIndex
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndex, typename TSpec = void>
struct GenomeIndex
{
    TGenome             & genome;
    TIndex              index;

    GenomeIndex(TGenome & genome) :
        genome(genome)
    {}
};

// ----------------------------------------------------------------------------
// Class ReadsIndex
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec = void>
struct ReadsIndex
{
    TReads              const & reads;
    TIndex              index;

    ReadsIndex(TReads const & reads) :
        reads(reads)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _indexCreate()                              [Uncompressed FM Index]
// ----------------------------------------------------------------------------

//namespace seqan {
//template <>
//inline bool _indexCreate(TGenomeFM & index, TGenome & text)
//{
//    if (empty(text)) return false;
//
//    // Create the suffix array table.
//    indexCreate(index, FibreSA());
//
//    // Create the lf table.
//    _indexCreateLfTables(index, text, indexSA(index));
//
//    return true;
//}
//}

// ----------------------------------------------------------------------------
// Function open()                                      [Uncompressed FM Index]
// ----------------------------------------------------------------------------

//namespace seqan {
//template <>
//inline bool open(TGenomeFM & index, const char * fileName, int openMode)
//{
//    String<char> name;
//
//    String<Pair<unsigned, Size<TGenomeFM>::Type> > infoString;
//
//    name = fileName;    append(name, ".txt");
//    if (!open(getFibre(index, FibreText()), toCString(name), openMode)) return false;
//
//    name = fileName;    append(name, ".sa");
//    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;
//
//    name = fileName;    append(name, ".lf");
//    if (!open(getFibre(index, FibreLfTable()), toCString(name), openMode)) return false;
//
//    name = fileName;    append(name, ".fma");
//    if (!open(infoString, toCString(name), openMode)) return false;
//
//    index.compressionFactor = infoString[0].i1;
//    index.n = infoString[0].i2;
//    //getFibre(index, FibreSA()).lfTable = & getFibre(index, FibreLfTable());
//
//    return true;
//}
//}

// ----------------------------------------------------------------------------
// Function _getNodeByChar() for Dna5Q                  [Iter<FMIndex, VSTree>]
// ----------------------------------------------------------------------------

namespace seqan {
    template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
    inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
                               typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type const & vDesc,
                               Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> & _range,
                               Dna5Q c)
    {
        typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
        typedef typename Value<TIndex>::Type                        TAlphabet;
        typedef typename ValueSize<TAlphabet>::Type                 TAlphabetSize;
        typedef typename Size<TIndex>::Type                         TSize;

        typedef typename Fibre<TIndex, FibreLfTable>::Type          TLfTable;
        typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

        if (__MASK_DNA5Q_LT[ordValue(c)] >= ValueSize<TAlphabet>::VALUE) return false;

        TIndex const & _index = container(it);
        TPrefixSumTable const & pst = getFibre(getFibre(_index, FibreLfTable()), FibrePrefixSumTable());

        TAlphabetSize cPosition = getCharacterPosition(pst, c);

        if (_isRoot(vDesc))
        {
            _range.i1 = getPrefixSum(pst, cPosition);
            _range.i2 = getPrefixSum(pst, cPosition + 1);
        }
        else
        {
            TSize prefixSum = getPrefixSum(pst, cPosition);
            _range.i1 = prefixSum + countOccurrences(_index.lfTable.occTable, c, vDesc.range.i1 - 1);
            _range.i2 = prefixSum + countOccurrences(_index.lfTable.occTable, c, vDesc.range.i2 - 1);
        }

        return _range.i1 + 1 <= _range.i2;
    }
}

// ----------------------------------------------------------------------------
// Function countOccurrences() for Dna5Q                          [WaveletTree]
// ----------------------------------------------------------------------------

//namespace seqan {
//template <typename TText, typename TSpec, typename TPos>
//inline unsigned countOccurrences(WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const & tree,
//                                 Dna5Q const character,
//                                 TPos const pos)
//{
//    typedef typename Value<TText>::Type             TAlphabet;
//
//    if (__MASK_DNA5Q_LT[ordValue(character)] >= ValueSize<TAlphabet>::VALUE) return 0;
//
//    unsigned occ = _countOccurrences(tree, character, pos);
//    if (ordEqual(getDollarSubstitute(tree), character))
//        occ -= getRank(getFibre(tree, FibreDollarPosition()), pos);
//
//    return occ;
//}
//}

// ----------------------------------------------------------------------------
// Function load()                                                [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndex, typename TSpec, typename TString>
bool load(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex, TString const & genomeIndexFile)
{
    genomeIndex.index = TIndex(genomeIndex.genome.contigs);

    return open(genomeIndex.index, toCString(genomeIndexFile));
}

// ----------------------------------------------------------------------------
// Function build()                                               [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndex, typename TSpec>
void build(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex)
{
    genomeIndex.index = TIndex(genomeIndex.genome.contigs);

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TIndex, TopDown<> >::Type it(genomeIndex.index);
}

template <typename TGenome, typename TSpec>
void build(GenomeIndex<TGenome, TGenomeFM, TSpec> & genomeIndex)
{
    // IndexFM is built on the reversed genome.
    reverse(genomeIndex.genome);

    genomeIndex.index = TGenomeFM(genomeIndex.genome.contigs);

    // Iterator instantiation calls automatic index construction.
    typename Iterator<TGenomeFM, TopDown<> >::Type it(genomeIndex.index);

    // NOTE(esiragusa): This only disables a warning.
    goRoot(it);

    reverse(genomeIndex.genome);
}

// ----------------------------------------------------------------------------
// Function dump()                                                [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TIndex, typename TSpec, typename TString>
bool dump(GenomeIndex<TGenome, TIndex, TSpec> & genomeIndex, TString const & genomeIndexFile)
{
    return save(genomeIndex.index, toCString(genomeIndexFile));
}

// ----------------------------------------------------------------------------
// Function clear()                                               [GenomeIndex]
// ----------------------------------------------------------------------------

template <typename TGenome, typename TSpec>
void clear(GenomeIndex<TGenome, TGenomeFM, TSpec> & genomeIndex)
{
    clear(genomeIndex.index);
}
// ----------------------------------------------------------------------------
// Function build()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec>
void build(ReadsIndex<TReads, TIndex, TSpec> & /* readsIndex */)
{
}

// ----------------------------------------------------------------------------
// Function clear()                                                [ReadsIndex]
// ----------------------------------------------------------------------------

template <typename TReads, typename TIndex, typename TSpec>
void clear(ReadsIndex<TReads, TIndex, TSpec> & /* readsIndex */)
{
}

#endif  // #ifndef SEQAN_EXTRAS_MASAI_INDEX_H_
