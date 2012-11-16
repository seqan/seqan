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
#include <seqan/file.h>
#include <seqan/index.h>
#include <seqan/index_fm.h>

#include "index/index_qgram_stretched.h"
#include "index/find_backtracking_stretched.h"

using namespace seqan;


// ============================================================================
// QGram Shape Definitions
// ============================================================================

typedef UngappedShape<10>                               TShape;


// ============================================================================
// Genome Index Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Genome Suffix Array Value Type Definition
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TGenome>
{
    typedef Pair<unsigned char, unsigned int, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// Genome Enhanced Suffix Array Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, IndexEsa<> >                     TGenomeEsa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
typedef DefaultIndexStringSpec<TGenomeEsa>::Type        TGenomeEsaStringSpec;
#else
typedef MMap<>                                          TGenomeEsaStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeEsa, FibreSA>
{
    typedef String<SAValue<TGenome>::Type, TGenomeEsaStringSpec>    Type;
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
// Genome Suffix Array Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, IndexSa<> >                      TGenomeSa;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
typedef DefaultIndexStringSpec<TGenomeSa>::Type         TGenomeSaStringSpec;
#else
typedef MMap<>                                          TGenomeSaStringSpec;
#endif

namespace seqan {
template <>
struct Fibre<TGenomeSa, FibreSA>
{
    typedef String<SAValue<TGenome>::Type, TGenomeSaStringSpec>     Type;
};
}

// ----------------------------------------------------------------------------
// Genome QGram Index with Bucket Refinement Type Definitions
// ----------------------------------------------------------------------------

typedef IndexQGram<TShape>                              TQGramBaseIndex;
typedef IndexQGram<TShape, BucketRefinement>            TQGramBucketRefinementIndex;

typedef Index<TGenome, TQGramBaseIndex>                 TGenomeBaseQGram;
typedef Index<TGenome, IndexSa<InfixSegment> >          TGenomeInfixSa;
typedef Index<TGenome, TQGramBucketRefinementIndex>     TGenomeQGram;

#if defined(SEQAN_EXTRAS_MASAI_DISABLE_MMAP) || defined(SEQAN_EXTRAS_MASAI_INDEXER_H_)
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

template <>
struct Fibre<TGenomeInfixSa, FibreSA>
{
    typedef Segment<Fibre<TGenomeBaseQGram, FibreSA>::Type const, InfixSegment>     Type;
};
}

// ----------------------------------------------------------------------------
// Genome FM Index Type Definitions
// ----------------------------------------------------------------------------

typedef Index<TGenome, FMIndex<WT<FmiDollarSubstituted<MultiDollar<void> > >, CompressText> > TGenomeFM;

// ----------------------------------------------------------------------------
// Genome FM Index Uncompressed Suffix Array Definitions
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
//
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
//
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
// _getNodeByChar() function for Dna5Q
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
// countOccurrences() function for Dna5Q
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

// ============================================================================
// Reads Index Type Definitions
// ============================================================================

// ----------------------------------------------------------------------------
// Reads Suffix Array Value Type Definition
// ----------------------------------------------------------------------------

namespace seqan {
template <>
struct SAValue<TReadSeqStore>
{
    typedef Pair<unsigned int, unsigned short, Pack> Type;
};
}

// ----------------------------------------------------------------------------
// Reads Wotd Type Definitions
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
// Reads QGram Index with Bucket Refinement Type Definitions
// ----------------------------------------------------------------------------

//typedef Index<TReadSeqStore, TQGramBaseIndex>                 TReadsBaseQGram;
//typedef Index<TReadSeqStore, IndexSa<InfixSegment> >          TReadsInfixSa;
//typedef Index<TReadSeqStore, TQGramBucketRefinementIndex>     TReadsQGram;
//
//namespace seqan
//{
//    template <>
//    struct Fibre<TReadsBaseQGram, FibreDir>
//    {
//        typedef String<unsigned int, DefaultIndexStringSpec<TReadsBaseQGram>::Type >   Type;
//    };
//
//    template <>
//    struct Fibre<TReadsInfixSa, FibreSA>
//    {
//        typedef Segment<Fibre<TReadsBaseQGram, FibreSA>::Type const, InfixSegment>           Type;
//    };
//}

// ----------------------------------------------------------------------------
// Reads QGram Index Type Definitions
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
// Reads QGram Index with Wotd Subtree Type Definitions
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


#endif  // #ifndef SEQAN_EXTRAS_MASAI_INDEX_H_
