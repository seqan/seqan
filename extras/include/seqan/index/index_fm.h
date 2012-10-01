// ==========================================================================
//                                  FMIndex
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

#ifndef INDEX_FM_H_
#define INDEX_FM_H_

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
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.FMIndex@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a FM index.
..cat:Index

..tag.FibrePrefixSumTable:The prefix sum table of the index.
..tag.FibreSA:The compressed suffix array of the text.
..tag.FibreText:The original text of the index.

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

template <typename TText, typename TStringSetSpec, typename TWaveletTreeSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TStringSetSpec>, FMIndex<WT<TWaveletTreeSpec>, TSpec > >, FibreOccTable>
{
     typedef WaveletTree<TText, FmiDollarSubstituted<MultiDollar<void> > > Type;
};
// ==========================================================================

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable const>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreOccTable>::Type        TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable>::Type  TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                        Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreLfTable const>
{
    typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreOccTable>::Type          TOccTable_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable>::Type    TPrefixSumTable_;
	typedef LfTable<TOccTable_, TPrefixSumTable_>   	                                                Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibrePrefixSumTable const>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type  TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                             TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                   Type;
};
    
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibrePrefixSumTable const>
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > const>::Type TChar_;
    typedef typename MakeUnsigned<TChar_>::Type                                 TUChar_;
	typedef PrefixSumTable<TUChar_, void>                                       Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreSA const>
{
	typedef unsigned int                                                                    TSAValue_;
	typedef SparseString<String<TSAValue_>, void>                                           TSparseString_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>::Type     TLfTable_;
	typedef CompressedSA<TSparseString_, TLfTable_, void>                                   Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreSA const>
{
	typedef unsigned int                                                                TSAValue_;
	typedef SparseString<String<TSAValue_>, void>                                       TSparseString_;
	typedef typename Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreLfTable>::Type TLfTable_;
	typedef CompressedSA<TSparseString_, TLfTable_, void> const                         Type;
};

template <typename TText, typename TSetSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> >, FibreSA>
{
    // TODO(singer): Note that we use a pair of unsigned values here at the moment instead of using SAValue to force using 32 bit integers. Is this what we ultimately want?
	typedef Pair<unsigned, unsigned>                                                                            TSAValue_;
	typedef SparseString<String<TSAValue_>, void>                                                               TSparseString_;
	typedef typename Fibre<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> >, FibreLfTable>::Type    TLfTable_;
	typedef CompressedSA<TSparseString_, TLfTable_, void>                                                       Type;
};

template <typename TText, typename TSetSpec, typename TOccSpec, typename TSpec>
struct Fibre<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> > const, FibreSA>
{
    typedef Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TSpec> >    TNonConstIndex_;
    typedef typename Fibre<TNonConstIndex_, FibreSA>::Type const            Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> >, FibreText>
{
	typedef TText Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Fibre<Index<TText, FMIndex<TOccSpec, TSpec> > const, FibreText const>
{
	typedef TText const Type;
};

// ==========================================================================
template <typename TText, typename TOccSpec, typename TSpec>
struct Size<Index<TText, FMIndex<TOccSpec, TSpec> > >
{
	typedef typename Size<TText>::Type Type;
};

template <typename TText, typename TOccSpec, typename TSpec>
struct Size<Index<TText, FMIndex<TOccSpec, TSpec> > const>
{
	typedef typename Size<TText>::Type const Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Tag.CompressText
..summary:Tag to select a FM index variant that can be used such that it is 
not necessary to store the text after index construction. 
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
..param.TOccSpec:FM index specialisation. 
...type:Tag.WT
...default.Class:WaveletTree<FmiDollarSubstituted>
..param.TSpec:Speed optimization method.
...type:Tag.CompressText
...default:void
..include:seqan/index.h
*/

//template <typename TText, typename TOccSpec , typename TSpec>
//class Index<TText, FMIndex<TOccSpec = typename FmiDollarSubstitutedDefault_<TText>::Type, TSpec > >;

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
		n(computeBwtLength_(text)),
		compressionFactor(compressionFactor)
	{
	    if (IsSameType<TSpec, CompressText>::VALUE)
    		indexCreate_(*this, text);
	}

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
unsigned computeBwtLength_(TText const & text)
{
    return(length(text) + 1);
}

// This function computes the length of the bwt string.
template <typename TText, typename TSetSpec>
unsigned computeBwtLength_(StringSet<TText, TSetSpec> const & text)
{
    return(lengthSum(text) + countSequences(text));
}

// ==========================================================================
// This function computes the BWT of a text. Note that the dollar sign is substituted and its position stored.
// The function is tested implicitly in lfTableLfMapping() in test_index_fm.h
template <typename TBwt, typename TDollarPosition, typename TText, typename TSA, typename TDollarSub>
inline void createBwTable_(
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
inline void createBwTable_(
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
    updateRanks_(dollarPos);
}

// ==========================================================================
// This function determines the '$' substitute.
// The character with the smallest number of occurrences greater 0 is chosen.
template <typename TPrefixSumTable, typename TChar>
inline void determineDollarSubstitute_(TPrefixSumTable const & pst,
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
// This function is used by the finder interface. It initializes the range of the findet.
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

	range_(index, pattern, finder.range);
}

// ==========================================================================
/**
.Function.getFibre:
..param.container:
...type:Spec.FMIndex
..param.fibreTag:
...type:Tag.FM Index Fibres
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

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreText >::Type & 
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > & index, FibreText /*tag*/)
{
	return value(index.text); 
}

template <typename TText, typename TIndexSpec, typename TFMISpeedEnhancement>
typename Fibre<Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> >, FibreText >::Type const &
getFibre(Index<TText, FMIndex<TIndexSpec, TFMISpeedEnhancement> > const & index, FibreText /*tag*/)
{
	return value(index.text); 
}

template <typename TText, typename TIndexSpec, typename TSpecSpec>
Nothing
getFibre(Index<TText, FMIndex<TIndexSpec, CompressText> > & /*tag*/, FibreText /*tag*/)
{
	return Nothing(); 
}

template <typename TText, typename TIndexSpec, typename TSpecSpec>
Nothing
getFibre(Index<TText, FMIndex<TIndexSpec, CompressText> > const & /*tag*/, FibreText /*tag*/)
{
	return Nothing(); 
}

// ==========================================================================
template <typename TText, typename TSetSpec, typename TFreq>
inline void getFrequencies_(TFreq & freq,
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

/*
 * This function determines the number of the different characters in the text.
 * In order to do so 'getNumCharsImpl' is called
 */
template <typename TText, typename TFreq>
inline void getFrequencies_(TFreq & freq,
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
inline bool indexCreateSA_(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TSA & fullSa, TText const & text)
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
inline bool indexCreateLfTables_(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, TText & text, TSA & sa)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> >		    TIndex;
	typedef typename Fibre<TIndex, FibreLfTable>::Type              TLfTable;
	typedef typename Fibre<TLfTable, FibreOccTable>::Type           TOccTable;
	typedef typename Fibre<TOccTable, FibreDollarPosition>::Type    TDollarPosition;
	typedef typename Value<TIndex>::Type						    TAlphabet;

	createPrefixSumTable(index.lfTable.prefixSumTable, text);

	TAlphabet dollarSub(0);
	determineDollarSubstitute_(index.lfTable.prefixSumTable, dollarSub);

	String<TAlphabet> bwt;
	resize(bwt, index.n);
	TDollarPosition dollarPos = 0;
	createBwTable_(bwt, dollarPos, text, sa, dollarSub);

	clear(sa);
    
    createOccurrenceTable(index.lfTable, bwt, dollarSub, dollarPos);
	clear(bwt);

	insertDollar_(index.lfTable.prefixSumTable, countSequences(text));

	return true;
}

// ==========================================================================
// This function creates the index.
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate_(Index<TText, FMIndex<TIndexSpec, TSpec > > & index,
		TText & text)
{
	typedef Index<TText, FMIndex<TIndexSpec, TSpec> > TIndex;
	typedef typename Fibre<TIndex, FibreSA>::Type TCompressedSA;
	typedef typename Fibre<TCompressedSA, FibreSparseString>::Type TSparseString;
	typedef typename Fibre<TSparseString, FibreValueString>::Type TValueString;

    if(empty(text))
        return false;

	// create the compressed SA
	TValueString fullSa;
	indexCreateSA_(index, fullSa, text);
	
	// create the lf table
	indexCreateLfTables_(index, text, fullSa);

	return true;
}


template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, FibreSaLfTable const)
{
    return indexCreate_(index, getFibre(index, FibreText()));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index)
{
    return indexCreate(index, FibreSaLfTable());
}
    
// template <typename TText, typename TIndexSpec, typename TSpec, typename TCompression>
// inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, TSpec> > & index, FibreSaLfTable const, TCompression compressionFactor)
// {
//     return indexCreate_(index, getFibre(index, FibreText()), compressionFactor);
// }

// TODO (singer): this needs to be documented. In addition, indexCreate should be extended
// to the different fibres.
template <typename TText, typename TIndexSpec>
inline bool indexCreate(Index<TText, FMIndex<TIndexSpec, CompressText> > & /*tag*/, FibreSaLfTable const)
{
    SEQAN_FAIL("Logic error. It is not possible to create this index without a text");
    return false;
}

// ==========================================================================
/**
.Function.indexSupplied:
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
template <typename TText, typename TPattern, typename TIter, typename TPairSpec, typename TOccSpec, typename TSpec>
inline void range_(const Index<TText, FMIndex<TOccSpec, TSpec> > & index,
		const TPattern & pattern,
		Pair<TIter, TPairSpec> & range)
{
    typedef unsigned TPos;

    if (empty(pattern))
    {
	    setPosition(range.i1, countSequences(index));
	    setPosition(range.i2, index.n);
    }

 	typedef typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type TCharValue;
 	typedef typename MakeUnsigned<TCharValue>::Type TUchar;
 	typedef typename Size<TPattern>::Type TSize;

	TSize i = length(pattern) - 1;
	TUchar letter = pattern[i];

	unsigned letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
	TPos sp = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
	TPos ep = getPrefixSum(index.lfTable.prefixSumTable, letterPosition + 1) - 1;
	
	while ((sp <= ep) && (i > 0))
	{
	    --i;
		letter = pattern[i];
		letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
		unsigned long prefixSum = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
		sp = prefixSum + getOccurrences(index.lfTable.occTable, letter, sp - 1);
		ep = prefixSum + getOccurrences(index.lfTable.occTable, letter, ep) - 1;
	}

    setPosition(range.i1, sp);
	setPosition(range.i2, ep + 1);
}

// ==========================================================================
// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(
    Index<TText, FMIndex<TOccSpec, TSpec> > & index,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<Pair<unsigned, typename Size<TText>::Type> > infoString;
    name = fileName;    append(name, ".txt");   open(getFibre(index, FibreText()), toCString(name), openMode);
    name = fileName;    append(name, ".sa");    open(getFibre(index, FibreSA()), toCString(name), openMode);
    name = fileName;    append(name, ".lf");    open(getFibre(index, FibreLfTable()), toCString(name), openMode);
    name = fileName;    append(name, ".fma");   open(infoString, toCString(name), openMode);
    
    index.compressionFactor = infoString[0].i1;
    index.n = infoString[0].i2;
    getFibre(index, FibreSA()).lfTable = & getFibre(index, FibreLfTable());
    return true;
}

// This function can be used to open a previously saved index.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool open(
    Index<TText, FMIndex<TOccSpec, TSpec> > & index,
    const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}

// ==========================================================================
// This function can be used to save an index on disk.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(
    Index<TText, FMIndex<TOccSpec, TSpec> > const & index,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<Pair<unsigned, typename Size<TText>::Type> > infoString;
    appendValue(infoString, Pair<unsigned, typename Size<TText>::Type>(index.compressionFactor, index.n));
    name = fileName;    append(name, ".txt");   save(getFibre(index, FibreText()), toCString(name), openMode);
    name = fileName;    append(name, ".sa");    save(getFibre(index, FibreSA()), toCString(name), openMode);
    name = fileName;    append(name, ".lf");    save(getFibre(index, FibreLfTable()), toCString(name), openMode);
    name = fileName;    append(name, ".fma");   save(infoString, toCString(name), openMode);
    return true;
}

// This function can be used to save an index on disk.
template <typename TText, typename TOccSpec, typename TSpec>
inline bool save(
    Index<TText, FMIndex<TOccSpec, TSpec> > const & index,
    const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TText, FMIndex<TOccSpec, TSpec> > >::VALUE);
}



// // This function computes a range in the suffix array who's entries point to location
// // in the text where the pattern occurs. 
// template <typename TText, typename TPattern, typename TPos, typename TOccSpec, typename TSpec>
// inline void range_(Index<TText, FMIndex<TOccSpec, TextVerification<TSpec> > > & index,
// 		const TPattern & pattern,
// 		Pair<TPos, TPos> & range)
// {
// 	typedef typename Alphabet<TPattern>::Type TCharValue;
// 	typedef typename MakeUnsigned<TCharValue>::Type TUChar;
// 	typedef typename Size<TPattern>::Type TSize;
// 
// 	int i = length(pattern) - 1;
// 	TUChar letter = pattern[i];
// 
// 	unsigned letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
// 	TPos sp = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
// 	TPos ep = getPrefixSum(index.lfTable.prefixSumTable, letterPosition + 1) - 1;
// 	--i;
// 
// 	while ((sp <= ep) && (i >= 0))
// 	{
// 		letter = pattern[i];
// 		letterPosition = getCharacterPosition(index.lfTable.prefixSumTable, letter);
// 		TPos prefixSum = getPrefixSum(index.lfTable.prefixSumTable, letterPosition);
// 		sp = prefixSum + getOccurrences(index.lfTable.occTable, letter, sp - 1);
// 		ep = prefixSum + getOccurrences(index.lfTable.occTable, letter, ep) - 1;
// 		--i;
// 	}
// 
// 	if(sp > ep || i == -1)
//     {
//     	range.i1 = sp;
// 	    range.i2 = ep;
// 	    return;
//     }
// 
//     bool match = true;
//     //bool posEntryStored = entryStored(index.compressedSA, sp);
// //     if(i >= 0)
// // 	{
// // 	    TUChar expectedChararacter = pattern[i];
// //         do
// //         {
// //             std::cerr << sp << " " << i << " " << getCharacter(getFibre(getFibre(index, FibreLfTable()), FibreOccTable()), sp) << " " << expectedChararacter << std::endl; 
// //             match = (expectedChararacter == getCharacter(getFibre(getFibre(index, FibreLfTable()), FibreOccTable()), sp));
// //             if (!match)
// //             {
// //                 range.i1 = ep + 1;
// // 	            range.i2 = ep;
// // 	            return;
// //             }
// // 
// //             posEntryStored = entryStored(index.compressedSA, sp);
// //             
// //             if(i == 0)
// //             {
// //                 index.compressedSA.sparseString.valueString[0] = index.compressedSA[sp] - i - 1;
// //                 setBit(index.compressedSA.sparseString.indicatorString, 0, 1); 
// //                 range.i1 = 0;
// //                 range.i2 = 0;
// //                 return;
// //             }
// // 
// //             if(posEntryStored)
// //                 break;
// //             
// //             sp = lfMapping(getFibre(index, FibreLfTable()), sp);
// //             --i;
// //             expectedChararacter = pattern[i];
// //         } while(match && !posEntryStored);
// //     }
// // 	else
// //     {
// // 		sp = lfMapping(getFibre(index, FibreLfTable()), sp);
// //     }
//     
// //     for(; !posEntryStored; --i)
//     for(;true; --i)
//     {
//         std::cerr << sp << " " << i << " " << getCharacter(getFibre(getFibre(index, FibreLfTable()), FibreOccTable()), sp) << " " << pattern[i] << std::endl; 
//         match = pattern[i] == getCharacter(getFibre(getFibre(index, FibreLfTable()), FibreOccTable()), sp);
// 
//         if (!match)
//         {
//             range.i1 = ep + 1;
//             range.i2 = ep;
//             return;
//         }
// 
//         //posEntryStored = entryStored(index.compressedSA, sp);
//         
//         if(i == 0)
//         {
//             std::cerr << "index.compressedAS: " << index.compressedSA[sp - 1] << " " << index.compressedSA[sp + 1] << " "<< index.compressedSA[sp] << " "<< getDollarPosition(index.lfTable.occTable) << std::endl;
//             index.compressedSA.sparseString.valueString[getDollarPosition(index.lfTable.occTable)] = index.compressedSA[sp] - 1;
//             setBit(index.compressedSA.sparseString.indicatorString, getDollarPosition(index.lfTable.occTable), 1); 
//             range.i1 = getDollarPosition(index.lfTable.occTable);
//             range.i2 = getDollarPosition(index.lfTable.occTable);
//             return;
//         }
// 
//         //posEntryStored = entryStored(index.compressedSA, sp);
//         sp = lfMapping(getFibre(index, FibreLfTable()), sp);
//     }
// 
//     std::cerr << "ii: " << i << " sp: " << sp << " ep: " << ep << std::endl;    
// 
//     if(match)
//     {
//        unsigned posInText = index.compressedSA[sp];
//        std::cerr << posInText << std::endl;
//        if (infix(pattern, 0, i + 1) == infix(getFibre(index, FibreText()), posInText - i - 1, posInText))
//        {
//            index.compressedSA.sparseString.valueString[0] = posInText - i - 1;
//            setBit(index.compressedSA.sparseString.indicatorString, 0, 1); 
//            sp = 0;
//            ep = sp;
//        }
//        else 
//        {
//            sp = ep + 1;
//        }
//     }
//     
// 	range.i1 = sp;
// 	range.i2 = ep;
//     
//     return;
// }

}
#endif
