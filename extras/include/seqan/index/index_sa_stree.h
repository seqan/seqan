// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_EXTRAS_INDEX_SA_STREE_H_
#define SEQAN_EXTRAS_INDEX_SA_STREE_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(esiragusa): Add IndexSa fibres.
//typedef FibreText         SaText;
//typedef FibreRawText      SaRawText;
//typedef FibreSA           SaSA;
//typedef FibreRawSA        SaRawSA;

template <typename TSpec = void>
struct IndexSa {};

/**
.Spec.IndexSa:
..summary:An index based on a suffix array.
..cat:Index
..general:Class.Index
..signature:Index<TText, IndexSa<> >
..param.TText:The text type.
...type:Class.String
...type:Class.StringSet
..include:seqan/index_extras.h
*/

template <typename TText, typename TSpec>
class Index<TText, IndexSa<TSpec> >
{
public:
    Holder<typename Fibre<Index, EsaText>::Type>    text;
    typename Fibre<Index, EsaSA>::Type              sa;
    typename Cargo<Index>::Type                     cargo;  // user-defined cargo

    Index() {}

    Index(Index & other) :
        text(other.text),
        sa(other.sa)
    {}

    Index(Index const & other) :
        text(other.text),
        sa(other.sa)
    {}

    template <typename TText_>
    Index(TText_ & _text) :
        text(_text)
    {}

    template <typename TText_>
    Index(TText_ const & _text) :
        text(_text)
    {}
};

template <typename TSize, typename TAlphabet>
struct VertexSA : public VertexEsa<TSize>
{
    typedef VertexEsa<TSize>                        TBase;

    TSize       repLen;
    TAlphabet   lastChar;

    VertexSA() :
        TBase(),
        repLen(0),
        lastChar(0)
    {}

    VertexSA(MinimalCtor) :
        TBase(MinimalCtor()),
        repLen(0),
        lastChar(0)
    {}

    VertexSA(VertexSA const & other) :
        TBase(other),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}
};

template <typename TSize, typename TAlphabet>
struct HistoryStackSA_
{
    Pair<TSize> range;
    TAlphabet   lastChar;

    HistoryStackSA_() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TSpec>
struct VertexDescriptor<Index<TText, IndexSa<TSpec> > >
{
    typedef Index<TText, IndexSa<TSpec> >           TIndex;
    typedef typename Size<TIndex>::Type             TSize;
    typedef typename Value<TIndex>::Type            TAlphabet;

    typedef VertexSA<TSize, TAlphabet>              Type;
};

template <typename TText, typename TIndexSpec, typename TSpec>
struct HistoryStackEntry_<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > >
{
private:
    typedef Index<TText, IndexSa<TIndexSpec> >      TIndex;
    typedef typename Size<TIndex>::Type             TSize;
    typedef typename Value<TIndex>::Type            TAlphabet;

public:
    typedef HistoryStackSA_<TSize, TAlphabet>       Type;
};

template <typename TText, typename TIndexSpec, typename TSpec>
struct EdgeLabel<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > >
{
    typedef typename Value<Index<TText, IndexSa<TIndexSpec> > >::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TIndexSpec>
void _indexRequireTopDownIteration(Index<TText, IndexSa<TIndexSpec> > & index)
{
    indexRequire(index, EsaSA());
}

template <typename TText>
void _indexRequireTopDownIteration(Index<TText, IndexSa<InfixSegment> > &)
{
    // The SA fibre must be provided by calling setHost explicitely.
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename SAValue<Index<TText, IndexSa<TIndexSpec> > >::Type
_lastOccurrence(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    if (_isSizeInval(value(it).range.i2))
        return back(indexSA(container(it)));
    else
        return saAt(value(it).range.i2 - 1, container(it));
}


// is this the first child of the parent node?
template <typename TIndex, typename TSpec>
inline bool
isFirstChild(Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    if (!empty(it.history))
    {
        typename Size<TIndex>::Type parentLeft = back(it.history).range.i1;
        if (value(it).range.i1 == parentLeft)
            return true;

        TIndex const & index = container(it);
        return suffixLength(saAt(value(it).range.i1 - 1, index), index) + 1 == value(it).repLen;
    }
    else
        return value(it).range.i1 == 0;

}

// is this the last child of the parent node?
template <typename TIndex, typename TSpec>
inline bool
isLastChild(Iter<TIndex, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).range.i2 == value(it).parentRight;
}

template <typename TIndex, typename TSpec>
inline bool
isLastChild(Iter<TIndex, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    if (!empty(it.history))
        return value(it).range.i2 == back(it.history).range.i2;
    else
        return _isSizeInval(value(it).range.i2);

}

// go down and skip non-branching nodes
// this function emulates suffix tree traversal
template <typename TIndex, typename TSpec>
inline bool
goDownSkipSingletons(Iter<TIndex, VSTree< TopDown<TSpec> > > &it)
{
    do {
        if (!goDown(it))
            return false;
    } while (isLastChild(it));
    return true;
}

// is this a leaf? (hide empty $-edges)
template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                      TIndex;
    typedef typename SAValue<TIndex>::Type                          TOcc;

//    if (_isLeaf(value(it))) return true;

    TIndex const & index = container(it);

    typename Size<TIndex>::Type lcp = repLength(it);

    // if the last suffix in the interval is larger than the lcp,
    // not all outgoing edges are empty (uses lex. sorting)
    TOcc oc = _lastOccurrence(it);
    return getSeqOffset(oc, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(oc, stringSetLimits(index)), index);
}

template <typename TIndex, typename TSize, typename TAlphabet>
inline typename Size<TIndex>::Type
repLength(TIndex const &, VertexSA<TSize, TAlphabet> const & vDesc)
{
    return vDesc.repLen;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Value<Index<TText, IndexSa<TIndexSpec> > >::Type
parentEdgeFirstChar(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).lastChar;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline void goRoot(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > & it)
{
    _historyClear(it);
    clear(it);
    if (!empty(indexSA(container(it))))
        _setSizeInval(value(it).range.i2);
    value(it).repLen = 0;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    // TODO(esiragusa): goDown including empty $-edges
    return _goDown(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    //std::cout << "goDown" << std::endl;
#endif

    // TODO(esiragusa): use HideEmptyEdges()
    if (_isLeaf(it, HideEmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

    // TODO(esiragusa): check nodeHullPredicate

#ifdef SEQAN_DEBUG
    //std::cout << "parent: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    Pair<typename Size<TIndex>::Type> saRange = range(it);

    // TODO(esiragusa): remove this check.
    if (saRange.i1 >= saRange.i2) return false;

    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // TODO(esiragusa): remove this check and ++saRange.i1 in loop.
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }

    // Get first and last characters in interval.
    TAlphabet cLeft = textAt(posAdd(saAt(saRange.i1, index), value(it).repLen), index);
    TAlphabet cRight = textAt(posAdd(saAt(saRange.i2 - 1, index), value(it).repLen), index);

#ifdef SEQAN_DEBUG
    //std::cout << "cLeft: " << cLeft << std::endl;
    //std::cout << "cRight: " << cRight << std::endl;
#endif

    // Save vertex descriptor.
    _historyPush(it);

    // Update left range.
    value(it).range.i1 = saRange.i1;

    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);

        TSAIterator upperBound = _upperBoundSA(text, node, cLeft, value(it).repLen);

        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }
    
    // Update child repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

#ifdef SEQAN_DEBUG
    //std::cout << "child: " <<  value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goRight(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    //std::cout << "goRight" << std::endl;
#endif

    if (isRoot(it))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);

#ifdef SEQAN_DEBUG
    //std::cout << "current: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    Pair<typename Size<TIndex>::Type> saRange;
    saRange.i1 = value(it).range.i2;
    saRange.i2 = (_isSizeInval(value(it).parentRight)) ? length(sa) : value(it).parentRight;

    if (saRange.i1 >= saRange.i2) return false;

    // Change repLen to parent repLen.
    value(it).repLen--;

    // TODO(esiragusa): don't check for empty edges (do it in goDown)
    // Skip $-edges.
    while (suffixLength(saAt(saRange.i1, index), index) <= value(it).repLen)
    {
        // Interval contains only $-edges.
        if (++saRange.i1 >= saRange.i2)
            return false;
    }

    // Get first and last characters in interval.
    TAlphabet cLeft = textAt(posAdd(saAt(saRange.i1, index), value(it).repLen), index);
    TAlphabet cRight = textAt(posAdd(saAt(saRange.i2 - 1, index), value(it).repLen), index);

    SEQAN_ASSERT_NEQ(ordValue(cLeft), ordValue(value(it).lastChar));

#ifdef SEQAN_DEBUG
    //std::cout << "cLeft: " << cLeft << std::endl;
    //std::cout << "cRight: " << cRight << std::endl;
#endif

    // Update left range.
    value(it).range.i1 = saRange.i1;

    // Update right range.
    // NOTE(esiragusa): I should use ordLess(cLeft, cRight) but masai redefines it for Ns.
    if (ordValue(cLeft) != ordValue(cRight))
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRange.i1;
        TSASize saLen = saRange.i2 - saRange.i1;
        TSearchTreeIterator node(saBegin, saLen);

        TText const & text = indexText(index);
        TSAIterator upperBound = _upperBoundSA(text, node, cLeft, value(it).repLen);

        value(it).range.i2 = upperBound - begin(sa, Standard());
    }
    else
    {
        value(it).range.i2 = saRange.i2;
    }

    // Update repLen, lastChar.
    value(it).repLen++;
    value(it).lastChar = cLeft;

#ifdef SEQAN_DEBUG
    //std::cout << "sibling: " <<  value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TValue>
inline bool _goDownChar(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, TValue c)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TIndex>::Type                     TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    //typedef typename Value<TIndex>::Type                    TAlphabet;

    if (_isLeaf(it, HideEmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

#ifdef SEQAN_DEBUG
    //std::cout << "parent: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);

    Pair<TSAIterator> range = _equalRangeSA(text, node, c, value(it).repLen);

    if (range.i1 >= range.i2)
        return false;

    // Save vertex descriptor.
    _historyPush(it);
    
    // Update range, lastChar and repLen.
    value(it).range.i1 = range.i1 - begin(sa, Standard());
    value(it).range.i2 = range.i2 - begin(sa, Standard());
    value(it).lastChar = c;
    value(it).repLen++;
    
#ifdef SEQAN_DEBUG
    //std::cout << "child: " <<  value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool _goDownString(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                          TString const & pattern, TSize & lcp)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TSA const>::Type                  TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;

    if (_isLeaf(it, HideEmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

#ifdef SEQAN_DEBUG
    //std::cout << "parent: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);
    Pair<TSAIterator> range = _equalRangeSA(text, node, pattern, value(it).repLen);

    if (range.i1 >= range.i2)
        return false;

    // Save vertex descriptor.
    _historyPush(it);
    
    // Update range, lastChar and repLen.
    value(it).range.i1 = range.i1 - begin(sa, Standard());
    value(it).range.i2 = range.i2 - begin(sa, Standard());
    value(it).lastChar = back(pattern);
    value(it).repLen += length(pattern);

#ifdef SEQAN_DEBUG
    //std::cout << "child: " <<  value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    lcp = length(pattern);

    return true;
}

// return vertex descriptor of parent's node
template <typename TText, typename TIndexSpec, typename TSpec>
inline typename VertexDescriptor<Index<TText, IndexSa<TIndexSpec> > >::Type
nodeUp(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > const & it)
{
    typedef Index<TText, IndexSa<TIndexSpec> > TIndex;
    
    if (!empty(it.history))
    {
        typename VertexDescriptor<TIndex>::Type desc;
        desc.range = back(it.history).range;
        desc.lastChar = back(it.history).lastChar;
        desc.repLen = value(it).repLen - 1;
        if (length(it.history) >= 2)
            desc.parentRight = topPrev(it.history).range.i2;
        else
            desc.parentRight = value(it).parentRight;
        return desc;
    } else
        return value(it);
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool _goUp(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
#ifdef SEQAN_DEBUG
    //std::cout << "goUp" << std::endl;
#endif

    if (!empty(it.history))
    {
        value(it).range = back(it.history).range;
        value(it).lastChar = back(it.history).lastChar;
        value(it).repLen--;
        pop(it.history);
        if (!empty(it.history))
            value(it).parentRight = back(it.history).range.i2;
        return true;
    }
    return false;
}

// ============================================================================

template <typename TText, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
    value(it).parentRight = value(it).range.i2;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;
    typename HistoryStackEntry_<TIter>::Type h;
    h.range = value(it).range;
    h.lastChar = value(it).lastChar;

    value(it).parentRight = value(it).range.i2;
    appendValue(it.history, h);
}

// ============================================================================

template <typename TText, typename TSpec>
inline void clear(Index<TText, IndexSa<TSpec> > & index)
{
    clear(getFibre(index, EsaSA()));
}

template <typename TObject, typename TSpec>
inline bool open(Index<TObject, IndexSa<TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;    append(name, ".txt");
    if ((!open(getFibre(index, EsaText()), toCString(name), openMode)) &&
        (!open(getFibre(index, EsaText()), fileName, openMode))) return false;

    name = fileName;    append(name, ".sa");
    if (!open(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TSpec>
inline bool open(Index<TObject, IndexSa<TSpec> > & index, const char * fileName)
{
    return open(index, fileName, DefaultOpenMode<Index<TObject, IndexSa<TSpec> > >::VALUE);
}

template <typename TObject, typename TSpec>
inline bool save(Index<TObject, IndexSa<TSpec> > & index, const char * fileName, int openMode)
{
    String<char> name;
    
    name = fileName;    append(name, ".txt");
    if ((!save(getFibre(index, EsaText()), toCString(name), openMode)) &&
        (!save(getFibre(index, EsaText()), fileName, openMode))) return false;

    name = fileName;    append(name, ".sa");
    if (!save(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

    return true;
}

template <typename TObject, typename TSpec>
inline bool save(Index<TObject, IndexSa<TSpec> > & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TObject, IndexSa<TSpec> > >::VALUE);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INDEX_SA_STREE_H_
