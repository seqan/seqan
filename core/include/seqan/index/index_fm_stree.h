// ==========================================================================
//                 seqan - the library for sequence analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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

#ifndef INDEX_FM_STREE_H_
#define INDEX_FM_STREE_H_

namespace seqan {
// ==========================================================================
// Classes
// ==========================================================================

template <typename TAlphabet, typename TSize>
struct VertexFmi
{
    Pair<TSize> range;
    TSize       repLen;
    TAlphabet   lastChar;

    VertexFmi() :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexFmi(MinimalCtor) :
        range(0, 0),
        repLen(0),
        lastChar(0)
    {}

    VertexFmi(Pair<TSize> newCurrentRange, TSize newRepLen, TAlphabet newChar) :
        range(newCurrentRange),
        repLen(newRepLen),
        lastChar(newChar)
    {}

    VertexFmi(VertexFmi const & other) :
        range(other.range),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}

    inline VertexFmi &
    operator = (VertexFmi const & _origin)
    {
        range = _origin.range;
        repLen = _origin.repLen;
        lastChar = _origin.lastChar;
        return *this;
    }
};

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec>
struct VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >
{
    typedef Index<TText,FMIndex<TOccSpec, TIndexSpec> >     TIndex;
	typedef typename Value<TIndex>::Type                    TAlphabet;
    typedef typename Size<TIndex>::Type                     TSize;

    typedef VertexFmi<TAlphabet, TSize> Type;
};

template <typename TAlphabet, typename TSize>
struct HistoryStackFmi_
{
    Pair<TSize> range;		// current SA interval of hits
    TAlphabet   lastChar;

    HistoryStackFmi_() {}

    template <typename TAlphabet_, typename TSize_>
    HistoryStackFmi_(TAlphabet_ const _lastChar, Pair<TSize_> const &_range):
        range(_range),
        lastChar(_lastChar)
    {}

    inline HistoryStackFmi_ const &
    operator=(HistoryStackFmi_ const & _origin)
    {
        range = _origin.range;
        lastChar = _origin.lastChar;
    }
};

template <typename TText, typename TOccSpec, typename TSpec, typename TIterSpec>
struct HistoryStackEntry_<Iter<Index<TText, FMIndex<TOccSpec, TSpec> >,
                               VSTree< TopDown< ParentLinks<TIterSpec> > > > >
{
    typedef HistoryStackFmi_<typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type,
                             typename Size<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type>     Type;
};


// ============================================================================
// Metafunctions
// ============================================================================
    
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
struct Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, TopDown<TSpec> >
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef Iter<TIndex, VSTree< TopDown<TSpec> > >         Type;
};

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
struct Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> > const, TopDown<TSpec> >
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> > const  TIndex;
    typedef Iter<TIndex, VSTree< TopDown<TSpec> > >             Type;
};

// ============================================================================

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type Type;
};


// ============================================================================
// Functions
// ============================================================================

// ==========================================================================
///.Function.begin.param.type:Spec.FMIndex
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, TSpec>::Type TIter;

	TIter it(index);
	value(it).range.i1 = index.lfTable.prefixSumTable[0];

	return it;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TIndexSpec> > const & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TIndexSpec> > const, TSpec>::Type TIter;

	TIter it(index);
	value(it).range.i1 = index.lfTable.prefixSumTable[0];

	return it;
}

// ==========================================================================
///.Function.clear.param.type:Spec.FMIndex
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline void clear(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TSpec> > &it)
{
    _clear(it);
}

// ==========================================================================
template <typename TAlphabet, typename TSize>
inline bool _isRoot(VertexFmi<TAlphabet, TSize> const &value)
{
    return _isSizeInval(value.range.i2);
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec>
void _indexRequireTopDownIteration(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index) 
{
    indexRequire(index, FibreSaLfTable());
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 + 1 >= value(it).range.i2 &&
                value(it).range.i1 == getDollarPosition(container(it).lfTable.occTable));
}

template <typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 + 1 >= value(it).range.i2 &&
                dollarPosition(getFibre(getFibre(container(it), FibreLfTable()), FibreOccTable()), value(it).range.i1));
}

// ==========================================================================
// This function the node of the corresponding character
//template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TValue >
//    inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TSpec> > & it, TValue c)

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline bool _goDownChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                        TChar const & c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename MakeUnsigned<TAlphabet>::Type              TUChar;
    typedef typename Size<TIndex>::Type                         TSize;


    if (isLeaf(it)) return false;

    TIndex const & _index = container(it);

    TSize sp, ep;

    // TODO(esiragusa): Remove cast to TAlphabet
    unsigned cPosition = getCharacterPosition(_index.lfTable.prefixSumTable, (TAlphabet)c);

    if (isRoot(it))
    {
        sp = getPrefixSum(_index.lfTable.prefixSumTable, cPosition);
        ep = getPrefixSum(_index.lfTable.prefixSumTable, cPosition + 1);
    }
    else
    {
        TSize prefixSum = getPrefixSum(_index.lfTable.prefixSumTable, cPosition);
        sp = prefixSum + countOccurrences(_index.lfTable.occTable, (TAlphabet)c, value(it).range.i1 - 1);
        ep = prefixSum + countOccurrences(_index.lfTable.occTable, (TAlphabet)c, value(it).range.i2 - 1);
    }

    if (sp + 1 > ep)
        return false;

    _historyPush(it);

    value(it).range.i1 = sp;
    value(it).range.i2 = ep;
    value(it).lastChar = (TAlphabet)c;
    value(it).repLen++;

    return true;
}

// ==========================================================================

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename Value<TIndex>::Type                    TAlphabet;

    if (isLeaf(it)) return false;

    for (unsigned char c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
        if (_goDownChar(it, c)) return true;

    return false;
}

// ==========================================================================

// TODO(esiragusa): Reimplement _goDownString()
// NOTE(esiragusa): Pushing vertex descriptors on the history for each char is costly.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                          TString const & string,
                          TSize & lcp)
{
    typedef Index<TText, FMIndex<TOccSpec, CompressText> >      TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename Iterator<TString const, Standard>::Type    TStringIter;

    lcp = 0;

    for (TStringIter stringIt = begin(string, Standard()); stringIt != end(string, Standard()); ++stringIt, ++lcp)
        if (!_goDownChar(it, value(stringIt)))
            return false;

    return true;
}

// ==========================================================================

// TODO(esiragusa): Implement this.
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                     VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goRight(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                     VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;

    typedef typename Fibre<TIndex, FibreLfTable>::Type          TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

    TPrefixSumTable & pst = getFibre(getFibre(container(it), FibreLfTable()), FibrePrefixSumTable());

    TAlphabet lastChar = value(it).lastChar;

    // TODO(esiragusa): Rewrite it without calling _goUp()
    if (_goUp(it))
    {
	    for (unsigned c = getCharacterPosition(pst, lastChar) + 1; c < getAlphabetSize(pst); ++c)
	    	if (_goDownChar(it, getCharacter(pst, c)))
		    	return true;
        _goDownChar(it, lastChar);
	}
	
	return false;
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
	if (!isRoot(it))
	{
        value(it).range = it._parentDesc.range;
        value(it).lastChar = it._parentDesc.lastChar;
        --value(it).repLen;
        return true;
    }
    
    return false;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it)
{
	if (!isRoot(it))
	{
        value(it).range = back(it.history).range;
        value(it).lastChar = back(it.history).lastChar;
        --value(it).repLen;
        pop(it.history);
        return true;
    }
    
    return false;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, class TSpec>
inline typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type
nodeUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown< ParentLinks<TSpec> > > > const & it)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename VertexDescriptor<TIndex>::Type         TVertexDescriptor;

    if (!empty(it.history))
        return TVertexDescriptor(back(it.history).range, value(it).repLen - 1, back(it.history).lastChar);
    else
        return value(it);
}

// ==========================================================================

// NOTE(esiragusa): this implementation of getOccurrences() might be very inefficient
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
String<typename SAValue<TText>::Type>
getOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
    String<typename SAValue<TText>::Type> bla;
    for (unsigned i = value(it).range.i1; i <= value(it).range.i2; ++i)
        appendValue(bla, container(it).compressedSA[i]);
    return bla;
}

// NOTE(esiragusa): this implementation of getOccurrences() might be very inefficient
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
String<typename SAValue<TText>::Type>
getOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
    String<typename SAValue<TText>::Type> bla;
    for (unsigned i = value(it).range.i1; i <= value(it).range.i2; ++i)
        appendValue(bla, container(it).compressedSA[i]);
    return bla;
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > &it)
{
    it._parentDesc = value(it);
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline void
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > &it)
{
    typedef Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;

    typename HistoryStackEntry_<TIter>::Type h;

    h.range = value(it).range;
    h.lastChar = value(it).lastChar;

    appendValue(it.history, h);
}

// ==========================================================================
template <typename TIndex, typename TAlphabet, typename TSize>
inline typename Size<TIndex>::Type repLength(TIndex const &, VertexFmi<TAlphabet, TSize> const & vDesc)
{
	return vDesc.repLen;
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline typename Value<Index<TText, FMIndex<TOccSpec, TIndexSpec > > >::Type
parentEdgeFirstChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

}
#endif  // INDEX_FM_STREE_H_
