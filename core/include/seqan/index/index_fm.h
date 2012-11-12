// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_H
#define INDEX_FM_H

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

template <typename TChar, typename TSpec>
class PrefixSumTable;

template <typename TSpec>
struct FmiDollarSubstitutedDefault_;

// WT = WaveletTree
/**
.Tag.WT
..summary:Tag that specifies the @Spec.FMIndex@ to use a wavelet tree as the occurrence table.
..cat:Index
*/
template <typename TSpec = void>
class WT;

template <typename TOccSpec = WT<>, typename TSpec = void>
class FMIndex;

struct FinderFMIndex_;
typedef Tag<FinderFMIndex_> FinderFMIndex;

template <typename TText, typename TOccSpec, typename TSpec>
struct DefaultFinder<Index<TText, FMIndex<TOccSpec, TSpec> > >
{
	typedef FinderFMIndex Type;
};

// FM index fibres
/**
.Tag.FM Index Fibres
..summary:Tag to select a specific fibre of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibrePrefixSumTable:The prefix sum table of the index.
..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.
..tag.FibreLfTable:The lf table.
..tag.FibreSaLfTable:The lf table as well as the compressed suffix array.
...remarks:This tag can only be used with the functions @Function.indexRequire@ or @Function.indexSupplied@.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index_fm.h
*/

struct FibrePrefixSumTable_;
struct FibreSA_;
struct FibreText_;
struct FibreLfTable_;
struct FibreSaLfTable_;
struct CompressText_;

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

typedef Tag<FibrePrefixSumTable_> const     FmiPrefixSumTable;  typedef FmiPrefixSumTable   FibrePrefixSumTable;
typedef Tag<FibreSA_> const                 FmiCompressedSA;    typedef FmiCompressedSA     FibreSA;
typedef Tag<FibreText_> const               FmiText;            typedef FmiText             FibreText;
typedef Tag<FibreLfTable_> const            FmiLfTable;         typedef FmiLfTable          FibreLfTable;
typedef Tag<FibreSaLfTable_> const          FmiSaLfTable;       typedef FmiSaLfTable        FibreSaLfTable;
typedef Tag<CompressText_> const            CompressText;

// ==========================================================================
// Metafunctions
// ==========================================================================
// Wavelet tree based FM index fibres
template <typename TText, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<WT<TWaveletTreeSpec>, TSpec> >, FibreOccTable>
{
	typedef WaveletTree<TText, FmiDollarSubstituted<SingleDollar<void> > > Type;
};

template <typename TText, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<WT<TWaveletTreeSpec>, TSpec> > const, FibreOccTable>
{
    typedef WaveletTree<TText, FmiDollarSubstituted<SingleDollar<void> > > const Type;
};
    
template <typename TText, typename TStringSetSpec, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TStringSetSpec>, FMIndex<WT<TWaveletTreeSpec>, TSpec > >, FibreOccTable>
{
     typedef WaveletTree<TText, FmiDollarSubstituted<MultiDollar<void> > > Type;
};

template <typename TText, typename TStringSetSpec, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TStringSetSpec>, FMIndex<WT<TWaveletTreeSpec>, TSpec > > const, FibreOccTable>
{
    typedef WaveletTree<TText, FmiDollarSubstituted<MultiDollar<void> > > const Type;
};
// ==========================================================================

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreOccTable>::Type        TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable>::Type  TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                        Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreLfTable>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreOccTable>::Type          TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable>::Type    TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                                Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type  TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                             TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                   Type;
};
    
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > const>::Type TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                                 TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                       Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
	typedef typename SAValue<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type                TSAValue_;
	typedef SparseString<String<TSAValue_>, void>                                           TSparseString_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>::Type     TLfTable_;
	typedef CompressedSA<TSparseString_, TLfTable_, void>                                   Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreSA> 
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA>::Type const Type;
};

template <typename TText, typename TSetSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> > const, FibreSA>
{
    typedef Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> >    TNonConstIndex_;
    typedef typename Fibre<TNonConstIndex_, FibreSA>::Type const            Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Tag.CompressText
..summary:Tag to select a FM index variant that can be used such that it is 
not necessary to store the text after index construction. This index is very
space efficient.
*/

// TODO (singer): type of TOccSpec
/**
.Spec.FMIndex:
..summary:An index based on the Burrows-Wheeler transform.
..cat:Index
..general:Class.Index
..signature:Index<TText, FMIndex<TOccSpec, TSpec> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..param.TOccSpec:Occurrence table specialisation. 
...type:Tag.WT
...default.Class:WaveletTree<FmiDollarSubstituted>
..param.TSpec:FM index specialisation.
...type:Tag.CompressText
...default:void
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
class Index<TText, FMIndex<TOccSpec, TSpec> >
{
    public:
    Holder<typename Fibre<Index, FibreText>::Type>  text;
	typename Fibre<Index, FibreLfTable>::Type       lfTable;
	typename Fibre<Index, FibreSA>::Type            compressedSA;
	typename Size<TText>::Type                      n;
    unsigned                                        compressionFactor; 

	Index() :
		n(0)
	{}

	Index(TText & text, unsigned compressionFactor = 10) :
	    text(text),
		n(_computeBwtLength(text)),
		compressionFactor(compressionFactor)
	{}

	inline bool operator==(const Index & b) const
    {
        return lfTable == b.lfTable &&
               compressedSA == b.compressedSA &&
               n == b.n &&
               compressionFactor == b.compressionFactor;
    }
};

// ==========================================================================
// Functions
// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
inline void clear(Index<TText, FMIndex<TOccSpec, TSpec> > & index)
{
    clear(getFibre(index, FibreLfTable()));
    clear(getFibre(index, FibreSA()));
}

// ==========================================================================
// This function computes the length of the bwt string.
template <typename TText>
unsigned _computeBwtLength(TText const & text)
{
    return(length(text) + 1);
}

// This function computes the length of the bwt string.
template <typename TText, typename TSetSpec>
unsigned _computeBwtLength(StringSet<TText, TSetSpec> const & text)
{
    return(lengthSum(text) + countSequences(text));
}

// ==========================================================================
// This function computes the BWT of a text. Note that the dollar sign is substituted and its position stored.
// The function is tested implicitly in lfTableLfMapping() in test_index_fm.h
template <typename TBwt, typename TDollarPosition, typename TText, typename TSA, typename TDollarSub>
inline void _createBwTable(
		TBwt & bwt,
		TDollarPosition & dollarPos,
		TText const & text,
		TSA const & SA,
		TDollarSub const dollarSub)
{
	//typedefs
	typedef typename GetValue<TSA>::Type    TSAValue;
	typedef typename Size<TSA>::Type        TSize;

	//little helpers
	TSize n = length(text);
	bwt[0] = back(text);
	for (TSize i = 0; i < n; i++)
	{
		TSAValue sa = getValue(SA, i);
		if (sa != 0)
			bwt[i + 1] = getValue(text, sa - 1);
		else
		{
			bwt[i + 1] = dollarSub;
			dollarPos = i + 1;
		}
	}
}

// This function computes the BWT of a text. Note that the dollar sign is substituted and its position stored.
// The function is tested implicitly in lfTableLfMapping() in test_index_fm.h
template <typename TBwt, typename TDollarPosition, typename TText, typename TSetSpec, typename TSA, typename TDollarSub>
inline void _createBwTable(
		TBwt & bwt,
		TDollarPosition & dollarPos,
		StringSet<TText, TSetSpec> const & text,
		TSA const & SA,
		TDollarSub const dollarSub)
{
    typedef typename Value<TSA>::Type   TSAValue;
    typedef typename Size<TSA>::Type    TSize;

    // little Helpers
    TSize n = countSequences(text);
    TSize nn = lengthSum(text);

    resize(dollarPos, n + nn);
    
    // fill the dollar positions (they are all at the beginning of the bwt)
    for (TSize i = 0; i < n; ++i)
    {
        bwt[i] = back(text[n - 1 - i]);
    }

    // compute the rest of the bwt
    for (TSize i = 0; i < nn; ++i)
    {
        TSAValue sa;    // = SA[i];
        posLocalize(sa, getValue(SA, i), stringSetLimits(text));
        if (sa.i2 != 0)
            bwt[i + n] = text[sa.i1][sa.i2 - 1];
        else
        {
            bwt[i + n] = dollarSub;
            setBit(dollarPos, i + n, 1);
        }
    }

    // update the auxiliary rank support bit string information.
    _updateRanks(dollarPos);
}

// ==========================================================================
// This function determines the '$' substitute.
// The character with the smallest number of occurrences greater 0 is chosen.
template <typename TPrefixSumTable, typename TChar>
inline void _determineDollarSubstitute(TPrefixSumTable const & pst,
		TChar & sub)
{
    typedef typename RemoveConst<TPrefixSumTable>::Type TNonConstPrefixSumTable;
	typedef typename Value<TNonConstPrefixSumTable>::Type TValue;

	TValue min = getPrefixSum(pst, length(pst) - 1);
	unsigned pos = length(pst) - 1;
	for (unsigned i = 0; i < length(pst) - 1; ++i)
	{
		unsigned diff = pst[i + 1] - pst[i];
		if (diff != 0 && diff < min)
		{
			min = diff;
			pos = i;
		}
	}
	sub = getCharacter(pst, pos);
}

// ==========================================================================
// This function checks whether the index is empty.
template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
inline bool empty(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index)
{
    return empty(getFibre(index, FibreLfTable()))
            && empty(getFibre(index, FibreSA()));
}

// ==========================================================================
// This function is used by the finder interface. It initializes the range of the finder.
template <typename TText, typename TPattern, typename TIndexSpec, typename TFMISpeedEnhancement>
inline void
_findFirstIndex(Finder<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FinderFMIndex> & finder,
		TPattern const & pattern,
		FinderFMIndex const &)
{
	typedef Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >    TIndex;
	typedef typename Fibre<TIndex, FibreSA>::Type                       TCompressedSA;
	typedef typename Iterator<TCompressedSA const>::Type                TIterator;

	TIndex & index = haystack(finder);

	indexRequire(index, FibreSaLfTable());
    setContainer(finder.range.i1, getFibre(container(finder), FibreSA()));
    setContainer(finder.range.i2, getFibre(container(finder), FibreSA()));

	_range(index, pattern, finder.range);
}

// ==========================================================================
/**
.Function.FMIndex#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Spec.FMIndex
..cat:Index
..param.container:The container holding the fibre.
...type:Spec.FMIndex
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.FM Index Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
*/
template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreLfTable >::Type &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > & index, FibreLfTable /*tag*/)
{
	return index.lfTable;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreLfTable >::Type const &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index, FibreLfTable /*tag*/)
{
	return index.lfTable;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreSA >::Type &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > & index, FibreSA /*tag*/)
{
	return index.compressedSA;
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreSA >::Type const &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index, FibreSA /*tag*/)
{
	return index.compressedSA;
}

// ==========================================================================

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TPos, typename TSize>
inline typename SAValue<Index<TText, FMIndex<TOccSpec, TIndexSpec > > const>::Type
toSuffixPosition(Index<TText, FMIndex<TOccSpec, TIndexSpec > > const & index, TPos i, TSize offset)
{
    SEQAN_ASSERT_GEQ(suffixLength(i, index), offset);
    setSeqOffset(i, suffixLength(i, index) - offset);
    return i;
}

// ==========================================================================
// This function determines the number of the different characters in the text.
template <typename TText, typename TSetSpec, typename TFreq>
inline void _getFrequencies(TFreq & freq,
						    StringSet<TText, TSetSpec> const & text)
{
	typedef typename Value<TText>::Type TChar;
	typedef typename Value<TFreq>::Type TFreqValue;
	resize(freq, ValueSize<TChar>::VALUE, 0);

    typedef typename Size<TText>::Type TSize;
	for (TSize i = 0; i < length(text); ++i)
	    for (TSize j = 0; j < length(text[i]); ++j)
		    ++freq[getCharacterPosition(freq, text[i][j])];
}

// This function determines the number of the different characters in the text.
template <typename TText, typename TFreq>
inline void _getFrequencies(TFreq & freq,
		   	   	   	   	   TText const & text)
{
	typedef typename Value<TText>::Type TChar;
	typedef typename Value<TFreq>::Type TFreqValue;

	resize(freq, ValueSize<TChar>::VALUE, 0);

    typedef typename Size<TText>::Type TSize;
	for (TSize i = 0; i < length(text); ++i)
		++freq[getCharacterPosition(freq, text[i])];
}

// ==========================================================================
template <typename TText, typename TWaveletTreeSpec, typename TPrefixSumTable, typename TText2, typename TChar, typename TPos>
inline bool createOccurrenceTable(LfTable<WaveletTree<TText, TWaveletTreeSpec>, TPrefixSumTable> & lfTable, TText2 & text, TChar dollarSub, TPos const & dollarPos)
{
	waveletTreeCreate(lfTable, text, dollarSub, dollarPos);
	return true;
}

// ==========================================================================
// This function computes the full and compressed suffix array. 
// Note, in contrast to indexCreate(index, FibreSA()) the full suffix array is also computed.
template <typename TText, typename TIndexSpec, typename TSpec, typename TSA>
inline bool _indexCreateSA(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TSA & fullSa, TText const & text)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> > TIndex;
	typedef typename Fibre<TIndex, FibreSA>::Type TCompressedSA;
	typedef typename Fibre<TCompressedSA, FibreSparseString>::Type TSparseString;
	typedef typename Fibre<TSparseString, FibreIndicatorString>::Type TIndicatorString;
	typedef typename Size<TCompressedSA>::Type TSize;

    // TODO(singer): If there is a lfTable we do not need the Skew7
	// create the fulle sa
	resize(fullSa, length(text));
	createSuffixArray(fullSa,
			text,
			Skew7());

	// create the compressed sa
	TCompressedSA & compressedSA = getFibre(index, FibreSA());
    setLfTable(compressedSA, getFibre(index, FibreLfTable()));

    unsigned numDollar = countSequences(text);
    compressedSaCreate(compressedSA, fullSa, index.compressionFactor, numDollar); 

	return true;
}

// ==========================================================================
// This function creates all table of the lf table given a text and a suffix array.
template <typename TIndexSpec, typename TSpec, typename TText, typename TSA>
inline bool _indexCreateLfTables(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TText & text, TSA const & sa)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> >		        TIndex;
	typedef typename Fibre<TIndex, FibreLfTable>::Type              TLfTable;
	typedef typename Fibre<TLfTable, FibreOccTable>::Type           TOccTable;
	typedef typename Fibre<TOccTable, FibreDollarPosition>::Type    TDollarPosition;
	typedef typename Value<TIndex>::Type						    TAlphabet;

	createPrefixSumTable(index.lfTable.prefixSumTable, text);

	TAlphabet dollarSub(0);
	_determineDollarSubstitute(index.lfTable.prefixSumTable, dollarSub);

	String<TAlphabet> bwt;
	resize(bwt, index.n);
	TDollarPosition dollarPos = 0;
	_createBwTable(bwt, dollarPos, text, sa, dollarSub);

    createOccurrenceTable(index.lfTable, bwt, dollarSub, dollarPos);

	_insertDollar(index.lfTable.prefixSumTable, countSequences(text));

	return true;
}

// ==========================================================================
/**
.Function.FMIndex#indexCreate
..summary:Creates a specific @Metafunction.Fibre@.
..param.fibreTag:
...type:Tag.FM Index Fibres.tag.FibreSaLfTable
..remarks:If you call this function on the compressed text version of the FM index
you will get an error message: "Logic error. It is not possible to create this index without a text."
*/

// This function creates the index.
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool _indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec > > & index,
		TText & text)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> > TIndex;
	typedef typename Fibre<TIndex, FibreSA>::Type TCompressedSA;
	typedef typename Fibre<TCompressedSA, FibreSparseString>::Type TSparseString;
	typedef typename Fibre<TSparseString, FibreValueString>::Type TValueString;

    if (empty(text))
        return false;

    TValueString fullSa;

	// create the compressed SA
	_indexCreateSA(index, fullSa, text);

	// create the lf table
	_indexCreateLfTables(index, text, fullSa);

	return true;
}


template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, FibreSaLfTable const)
{
    return _indexCreate(index, getFibre(index, FibreText()));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index)
{
    return indexCreate(index, FibreSaLfTable());
}

// ==========================================================================
/**
.Function.indexSupplied:
..summary:Returns whether a specific @Metafunction.Fibre@ is present.
..param.fibreTag:
...type:Tag.FM Index Fibres
*/
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > & index, FibreSaLfTable const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLfTable())));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexSupplied(Index<TText, FMIndex<TIndexSpec, TSpec > > const & index, FibreSaLfTable const)
{
    return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreLfTable())));
}

// ==========================================================================
// This function computes a range in the suffix array who's entries point to location
// in the text where the pattern occurs. 
template <typename TText, typename TOccSpec, typename TSpec, typename TPattern, typename TIter, typename TPairSpec>
inline void _range(const Index<TText, FMIndex<TOccSpec, TSpec> > & index,
		const TPattern & pattern,
		Pair<TIter, TPairSpec> & range)
{
    typedef Index<TText, FMIndex<TOccSpec, TSpec> >             TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename ValueSize<TAlphabet>::Type                 TAlphabetSize;
    typedef typename Size<TIndex>::Type                         TSize;
 	typedef typename Value<TPattern>::Type                      TChar;

    if (empty(pattern))
    {
	    setPosition(range.i1, countSequences(index));
	    setPosition(range.i2, index.n);
    }

	TSize i = length(pattern) - 1;
	TChar letter = pattern[i];

    // initilization
	TAlphabetSize letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
	TSize sp = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
	TSize ep = getPrefixSum(index.lfTable.prefixSumTable, letterPosition + 1) - 1;
	
	// the search as proposed by Ferragina and Manzini
	while ((sp <= ep) && (i > 0))
	{
	    --i;
		letter = pattern[i];
		letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
		TSize prefixSum = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
		sp = prefixSum + countOccurrences(index.lfTable.occTable, letter, sp - 1);
		ep = prefixSum + countOccurrences(index.lfTable.occTable, letter, ep) - 1;
	}

    setPosition(range.i1, sp);
	setPosition(range.i2, ep + 1);
}

// ==========================================================================
/**
.Function.open
..param.object:
...type:Spec.FMIndex
*/
// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

    String<Pair<unsigned, typename Size<TText>::Type> > infoString;

    name = fileName;    append(name, ".txt");
    if (!open(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!open(getFibre(index, FibreLfTable()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!open(infoString, toCString(name), openMode)) return false;
    
    index.compressionFactor = infoString[0].i1;
    index.n = infoString[0].i2;
    getFibre(index, FibreSA()).lfTable = & getFibre(index, FibreLfTable());

    return true;
}

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(Index<TText, FMIndex<TOccSpec, TSpec> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

// ==========================================================================
// This function can be used to save an index on disk.
/**
.Function.save:
..class:Class.Index
..summary:Saves an object.
..cat:Index
..signature:save(object, fileName[, openMode])
..param.object:The object to be saved.
...type:Spec.FMIndex
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName, int openMode)
{
    String<char> name;

    String<Pair<unsigned, typename Size<TText>::Type> > infoString;
    appendValue(infoString, Pair<unsigned, typename Size<TText>::Type>(index.compressionFactor, index.n));

    name = fileName;    append(name, ".txt");
    if (!save(getFibre(index, FibreText()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, FibreSA()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".lf");
    if (!save(getFibre(index, FibreLfTable()), toCString(name), openMode)) return false;

    name = fileName;    append(name, ".fma");
    if (!save(infoString, toCString(name), openMode)) return false;

    return true;
}

// This function can be used to save an index on disk.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(Index<TText, FMIndex<TOccSpec, TSpec> > const & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

}
#endif // INDEX_FM_H
