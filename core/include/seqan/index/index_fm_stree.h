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

#include <algorithm>

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
        range(0, -1),
        repLen(0),
        lastChar(0)
    {}

    VertexFmi(MinimalCtor) :
        range(0, -1),
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
    TAlphabet c;
    Pair<TSize> range;		// current SA interval of hits
    
    HistoryStackFmi_() {}

    template <typename TAlphabet_, typename TSize_>
    HistoryStackFmi_(TAlphabet_ const _c, Pair<TSize_> const &_range): 
        c(_c),
        range(_range) 
    {}

    inline HistoryStackFmi_ const &
    operator=(HistoryStackFmi_ const & _origin)
    {
        c = _origin.c;
        range = _origin.range;
    }
};

template < typename TText, typename TOccSpec, typename TSpec, typename TIterSpec >
struct HistoryStackEntry_< Iter< Index<TText, FMIndex<TOccSpec, TSpec> >, VSTree< TopDown< ParentLinks<TIterSpec> > > > >
{
    typedef HistoryStackFmi_<typename Value<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type, typename Size<Index<TText, FMIndex<TOccSpec, TSpec> > >::Type>	Type;
};

// ==========================================================================
/**
.Spec.TopDown Iterator:
..param.TContainer:
...type:Spec.FMIndex
*/
template <typename TText, typename TOccSpec, typename TIterSpec>
class Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TIterSpec> > >
{
public:
   
    typedef Index<TText, FMIndex<TOccSpec, CompressText> >    TIndex;
    typedef typename VertexDescriptor<TIndex>::Type	        TVertexDesc;
    //TODO(singer): What is this
    typedef Iter									        iterator;
    
    TIndex const	*index;		    // container of all necessary tables
    TVertexDesc		vDesc;		    // current interval in suffix array and
    TText           representative; // the representative

    // pseudo history stack (to go up at most one node)
    TVertexDesc		_parentDesc;

    Iter() {}
    
    Iter(TIndex &_index):
        index(&_index)
    {
        _indexRequireTopDownIteration(_index);
        _initRepresentative(*this);
        goRoot(*this);
    }
};

// ==========================================================================
/**
.Spec.TopDownHistory Iterator:
..param.TContainer:
...type:Spec.FMIndex
*/
template <typename TText, typename TOccSpec, typename TIterSpec>
class Iter< Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree< TopDown< ParentLinks<TIterSpec> > > >:
    public Iter< Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree< TopDown<> > >
{
public:

    typedef Index<TText, FMIndex<TOccSpec, CompressText> >  TIndex;
    typedef Iter< TIndex, VSTree< TopDown<> > >		        TBase;
    typedef	typename HistoryStackEntry_<Iter>::Type         TStackEntry;
    typedef String<TStackEntry, Block<> >			        TStack;
    typedef Iter									        iterator;

    TStack			history;	// contains all previously visited intervals (allows to go up)

    Iter() {}
    
    Iter(TIndex &_index):
        TBase(_index) {}

    Iter(TIndex &_index, MinimalCtor):
        TBase(_index, MinimalCtor()) {}

    Iter(Iter const &_origin):
        TBase((TBase const &)_origin),
        history(_origin.history) {}

    inline Iter const &
    operator = (Iter const &_origin)
    {
        *(TBase*)(this) = _origin;
        history = _origin.history;
        return *this;
    }
};

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec>
struct Iterator<Index<TText, FMIndex<TOccSpec, TSpeedSpec> >, TopDown<TSpec> >
{
    typedef Index<TText, FMIndex<TOccSpec, TSpeedSpec> > TIndex;
    typedef Iter<TIndex, VSTree< TopDown<TSpec> > > Type;
};

template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec>
struct Iterator<Index<TText, FMIndex<TOccSpec, TSpeedSpec> > const, TopDown<TSpec> >
{
    typedef Index<TText, FMIndex<TOccSpec, TSpeedSpec> > const TIndex;
    typedef Iter<TIndex, VSTree< TopDown<TSpec> > > Type;
};

// ==========================================================================
///.Function.begin.param.type:Spec.FMIndex
template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec >
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TSpeedSpec> >, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TSpeedSpec> > & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TSpeedSpec> >, TSpec>::Type TIter;

	TIter it(index);
	value(it).range.i1 = index.lfTable.prefixSumTable[0];

	return it;
}

template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec >
inline
typename Iterator<Index<TText,FMIndex<TOccSpec, TSpeedSpec> > const, TSpec>::Type
begin(Index<TText, FMIndex<TOccSpec, TSpeedSpec> > const & index, TSpec const /*Tag*/)
{
	typedef typename Iterator<Index<TText, FMIndex<TOccSpec, TSpeedSpec> > const, TSpec>::Type TIter;

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
    clear(_getRepresentative(it));
}

// ==========================================================================
template < typename TAlphabet, typename TSize >
inline bool _isRoot(VertexFmi<TAlphabet, TSize> const &value)
{
    return !(value.repLen);
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec>
void _indexRequireTopDownIteration(Index<TText, FMIndex<TOccSpec, TIndexSpec> > & index) 
{
    indexRequire(index, FibreSaLfTable());
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder >
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder >
inline bool _isLeaf(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 == value(it).range.i2 &&
                value(it).range.i1 == getDollarPosition(container(it).lfTable.occTable));
}

template < typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder >
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template < typename TText, typename TSetSpec, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder >
inline bool _isLeaf(Iter<Index<StringSet<TText, TSetSpec>, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    return (value(it).range.i1 == value(it).range.i2 &&
                dollarPosition(getFibre(getFibre(container(it), FibreLfTable()), FibreOccTable()), value(it).range.i1));
}

// ==========================================================================
// This function the node of the corresponding character
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TValue >
inline bool
_getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TSpec> > &it,
               TValue c)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> > TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename MakeUnsigned<TAlphabet>::Type TUChar;
    typedef typename Size<TIndex>::Type                       TSize;


    if (isLeaf(it)) return false;

    TIndex const & _index = container(it);

    TSize sp, ep;

    unsigned cPosition = getCharacterPosition(_index.lfTable.prefixSumTable, (TAlphabet)c);
    if(value(it).repLen == 0)
    {
        sp = getPrefixSum(_index.lfTable.prefixSumTable, cPosition);
        ep = getPrefixSum(_index.lfTable.prefixSumTable, cPosition + 1) - 1;
    }
    else
    {
        TSize prefixSum = getPrefixSum(_index.lfTable.prefixSumTable, cPosition);
        sp = prefixSum + countOccurrences(_index.lfTable.occTable, (TAlphabet)c, value(it).range.i1 - 1);
        ep = prefixSum + countOccurrences(_index.lfTable.occTable, (TAlphabet)c, value(it).range.i2) - 1;
    }

    if(sp > ep)
        return false;

    _historyPush(it, (TAlphabet)c);
    it._parentDesc = value(it);
    value(it).range.i1 = sp;
    value(it).range.i2 = ep;
    value(it).lastChar = (TAlphabet)c;

    return true;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec >
inline void _initRepresentative(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TSpec> > &) {}

template < typename TText, typename TStringSetSpec, typename TOccSpec, typename TSpec >
inline void _initRepresentative(Iter<Index<StringSet<TText, TStringSetSpec>, FMIndex<TOccSpec, CompressText> >, VSTree<TSpec> > &it)
{
    appendValue(it.representative, TText(), Exact());
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec >
inline TText &
_getRepresentative(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TSpec> > &it)
{
    return it.representative;
}

template < typename TText, typename TStringSetSpec, typename TOccSpec, typename TSpec >
inline TText &
_getRepresentative(Iter<Index<StringSet<TText, TStringSetSpec>, FMIndex<TOccSpec, CompressText> >, VSTree<TSpec> > &it)
{
    return back(it.representative);
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec, typename TChar >
void _storeCharacter(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > /*tag*/,
                     TChar const /*tag*/) {}

template < typename TText, typename TOccSpec, typename TSpec, typename TChar >
void _storeCharacter(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TSpec> > > &it,
                     TChar const c)
{
    appendValue(_getRepresentative(it), c);
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder >
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec,TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template < typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec, typename TDfsOrder >
inline bool _goDown(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >    TIndex;
    typedef typename Value<TIndex>::Type                     TAlphabet;

    if (isLeaf(it)) return false;

    for (unsigned char c = 0; c < ValueSize<TAlphabet>::VALUE; ++c)
        if (_getNodeByChar(it, c))
        {
            _storeCharacter(it, c);
            return true;
        }

    return false;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec, typename TChar>
inline bool _goDownChar(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TSpec> > > &it,
                        TChar const & character)
{
    typedef Index<TText, FMIndex<TOccSpec, CompressText> > TIndex;
    typedef typename Value<TIndex>::Type                         TAlphabet;

    if (isLeaf(it)) return false;

    if (_getNodeByChar(it, character))
    {
        _storeCharacter(it, character);
        return true;
    }
    return false;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline bool _goDownChar(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                        TChar const & character)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> > TIndex;
    typedef typename Value<TIndex>::Type                         TAlphabet;
    
    if (isLeaf(it)) return false;

    if (_getNodeByChar(it, character))
    {
        return true;
    }

    return false;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString>
inline bool _goDownString(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
                        TString const & string)
{
    typedef Index<TText, FMIndex<TOccSpec, CompressText> >      TIndex;
    typedef typename Value<TIndex>::Type                        TAlphabet;
    typedef typename Iterator<TString const, Standard>::Type    TStringIter;

    for (TStringIter stringIt = begin(string, Standard()); stringIt != end(string, Standard()); ++stringIt)
        if (!goDown(it, value(stringIt)))
            return false;

    return true;
}

// ==========================================================================
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TObject >
inline bool 
_goDownObject(
    Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown<TSpec> > > &it, 
    TObject const &obj,
    False)
{
    return _goDownChar(it, obj);
}

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TObject >
inline bool 
_goDownObject(
    Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown<TSpec> > > &it, 
    TObject const &obj,
    True)
{
    return _goDownString(it, obj);
}

// ==========================================================================
// public interface for goDown(it, ...)
template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TObject >
inline bool
goDown(
    Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree< TopDown<TSpec> > > &it, 
    TObject const &obj) 
{
    return _goDownObject(it, obj, typename IsSequence<TObject>::Type());
}


// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpeedSpec, class TSpec >
inline void goRoot(Iter<Index<TText, FMIndex<TOccSpec, TSpeedSpec> >, VSTree<TSpec> > &it) 
{
	_historyClear(it);
	clear(it);			
	if (!empty(indexSA(container(it))))
	{
		value(it).range.i1 = countSequences(TText());
	    value(it).range.i2 = container(it).n - 1;
		it._parentDesc.range = value(it).range;
	}
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec, typename TDfsOrder >
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TSpeedSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _goRight(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template < typename TText, typename TOccSpec, typename TSpeedSpec, typename TSpec, typename TDfsOrder >
inline bool _goRight(Iter<Index<TText, FMIndex<TOccSpec, TSpeedSpec> >, VSTree<TopDown<TSpec> > > &it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, FMIndex<TOccSpec, TSpeedSpec> >    TIndex;
    typedef typename Value<TIndex>::Type                  TAlphabet;

    typedef typename Fibre<TIndex, FibreLfTable>::Type TLfTable;
    typedef typename Fibre<TLfTable, FibrePrefixSumTable>::Type TPrefixSumTable;

    TPrefixSumTable pst = getFibre(getFibre(container(it), FibreLfTable()), FibrePrefixSumTable());

    if(_goUp(it))
    {
	    for(unsigned i = getCharacterPosition(pst, value(it).lastChar) + 1; i < getAlphabetSize(pst); ++i)
	    	if(goDown(it, getCharacter(pst, i)))
		    	return true;
        goDown(it, value(it).lastChar);
	}
	
	return false;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec, typename TIndexSpec>
void _eraseCharacter(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > /*tag*/) 
{} 

template < typename TText, typename TOccSpec, typename TSpec>
void _eraseCharacter(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TSpec> > > &it)
{
    resize(_getRepresentative(it), value(it).repLen - 1);
} 

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
	if(!isRoot(it))
	{
	    _eraseCharacter(it);
        value(it).range = it._parentDesc.range;
        --value(it).repLen;
        return true;
    }
    return false;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
bool _goUp(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it)
{
	if(!isRoot(it))
	{
	    _eraseCharacter(it);
        value(it).range = back(it.history).range;
        value(it).lastChar = back(it.history).c; 
        pop(it.history);
        --value(it).repLen;
        return true;
    }
    return false;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type
countOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
	return value(it).range.i2 - value(it).range.i1 + 1;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type
countOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
	return value(it).range.i2 - value(it).range.i1 + 1;
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
String<typename SAValue<TText>::Type>
getOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
    String<typename SAValue<TText>::Type> bla;
    for (unsigned i = value(it).range.i1; i <= value(it).range.i2; ++i)
        appendValue(bla, container(it).compressedSA[i]);
    return bla;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
String<typename SAValue<TText>::Type>
getOccurrences(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
    String<typename SAValue<TText>::Type> bla;
    for (unsigned i = value(it).range.i1; i <= value(it).range.i2; ++i)
        appendValue(bla, container(it).compressedSA[i]);
    return bla;
}


template < typename TText, typename TOccSpec, typename TIndexSpec, typename TDesc >
inline Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type>
range(Index<TText, FMIndex<TOccSpec, TIndexSpec> > const & /*tag*/, TDesc const &desc)
{
    return Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type> (desc.range.i1, desc.range.i2 + 1);
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec>
inline Pair<typename Size<Index<TText, FMIndex<TOccSpec, TIndexSpec> > >::Type>
range(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TSpec> > const &it)
{
    return  Pair<typename Size<Index<TText, TSpec> >::Type>(value(it).range.i1, value(it).range.i2 + 1);
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline void 
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<TSpec> > > &it, TChar const & /*dummy*/) 
{
    it._parentDesc = value(it);
    ++value(it).repLen;
}

template < typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TChar>
inline void 
_historyPush(Iter<Index<TText, FMIndex<TOccSpec, TIndexSpec > >, VSTree<TopDown<ParentLinks<TSpec> > > > &it, TChar const & _c) 
{
    resize(it.history, length(it.history) + 1);
    it.history[length(it.history) - 1].range = value(it).range;
    it.history[length(it.history) - 1].c = _c;
    ++value(it).repLen;
}

// ==========================================================================
template < typename TIndex, typename TAlphabet, typename TSize>
inline unsigned
repLength(TIndex const & /*tag*/, VertexFmi<TAlphabet, TSize> const &vDesc)
{
	return vDesc.repLen;
}

// ==========================================================================
/**
.Function.pathLabel:
..summary:Returns a substring representing the path from root to $iterator$ node.
..cat:Index
..signature:representative(iterator)
..class:Spec.VSTree Iterator
..param.iterator:An iterator of a virtual string tree.
...type:Spec.VSTree Iterator
..returns:An @Spec.InfixSegment@ of the text of an index (see @Tag.ESA Index Fibres.EsaText@).
If $iterator$'s container type is $TIndex$ the return type is $Infix<Fibre<TIndex, EsaText>::Type const>::Type$.
..include:seqan/index.h
*/
template < typename TText, typename TOccSpec, typename TCompression, typename TSpec>
inline ModifiedString<typename Fibre<Index<TText, FMIndex<TOccSpec, TCompression> >, FibreText>::Type, ModReverse>
pathLabel(Iter<Index<TText, FMIndex<TOccSpec, TCompression> >, VSTree<TopDown<TSpec> > > &it)
{
    return ModifiedString<typename Fibre<Index<TText, FMIndex<TOccSpec, TCompression> >, FibreText>::Type, ModReverse >(representative(it));
}

template < typename TText, typename TOccSpec, typename TSpec>
inline typename Infix< typename Fibre<Index<TText, FMIndex<TOccSpec, CompressText> >, FibreText>::Type const >::Type
pathLabel(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TSpec> > > &it)
{
    return infixWithLength(it.representative, 0, repLength(it));
}

// ==========================================================================
template < typename TText, typename TOccSpec, typename TSpec>
inline typename Infix< typename Fibre<Index<TText, FMIndex<TOccSpec, CompressText> >, FibreText>::Type const >::Type
representative(Iter<Index<TText, FMIndex<TOccSpec, CompressText> >, VSTree<TopDown<TSpec> > > &it)
{
    return infixWithLength(_getRepresentative(it), 0, repLength(it));
}

}
#endif  // INDEX_FM_STREE_H_
