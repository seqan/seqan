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
 ..include:seqan/index.h
 */

template <typename TText, typename TSpec>
class Index<TText, IndexSa<TSpec> >
{
public:
    Holder<typename Fibre<Index, EsaText>::Type>    text;
    typename Fibre<Index, EsaSA>::Type              sa;

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
struct VertexSA :
    public VertexEsa<TSize>
{
    typedef VertexEsa<TSize>                        TBase;
    typedef Array<ValueSize<TAlphabet>::VALUE>      TDirSpec;
    typedef String<TSize, TDirSpec>                 TDir;

    TDir        dir;
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
        dir(other.dir),
        repLen(other.repLen),
        lastChar(other.lastChar)
    {}
};

template <typename TSize, typename TAlphabet>
struct HistoryStackSA_
{
    typedef Array<ValueSize<TAlphabet>::VALUE>      TDirSpec;
    typedef String<TSize, TDirSpec>                 TDir;

    TDir        dir;
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

// ============================================================================

template <typename TText>
struct Fibre<Index<TText, IndexSa<InfixSegment> >, FibreSA>
{
    typedef Segment<typename Fibre<Index<TText, IndexSa<> >, FibreSA>::Type const, InfixSegment>  Type;
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

// is this a leaf? (including empty $-edges)
template <typename TText, typename TIndexSpec, class TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, False> const)
{
    return _isLeaf(value(it)) && _isLeaf(it, VSTreeIteratorTraits<TDfsOrder, True>());
}

// is this a leaf? (hide empty $-edges)
template <typename TText, typename TIndexSpec, class TSpec, typename TDfsOrder>
inline bool _isLeaf(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >                                  TIndex;
    typedef typename Infix<typename Fibre<TIndex, EsaSA>::Type const>::Type   TOccs;
    typedef typename Iterator<TOccs, Standard>::Type                            TIter;

    TIndex const & index = container(it);

    typename Size<TIndex>::Type lcp = repLength(it);

    // if the last suffix in the interval is larger than the lcp,
    // not all outgoing edges are empty (uses lex. sorting)
    TOccs occs = getOccurrences(it);
    TIter oc = begin(occs, Standard()) + length(occs) - 1;
    return getSeqOffset(*oc, stringSetLimits(index)) + lcp == sequenceLength(getSeqNo(*oc, stringSetLimits(index)), index);
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, IndexSa<TIndexSpec> > >::Type
repLength(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return value(it).repLen;
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
    typedef typename Size<TSA const>::Type                  TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    std::cout << "goDown" << std::endl;
#endif

    if (_isLeaf(it, EmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

    // TODO(esiragusa): check nodeHullPredicate

    TSASize saRangeBegin = value(it).range.i1;
    TSASize saRangeEnd = isRoot(it) ? length(sa) : value(it).range.i2;

    typename VertexDescriptor<TIndex>::Type childDesc = value(it);
    resize(childDesc.dir, ValueSize<TAlphabet>::VALUE);
    arrayFill(begin(childDesc.dir, Standard()), end(childDesc.dir, Standard()), saRangeBegin);

    // Skip $-edges.
    while (suffixLength(saAt(saRangeBegin, index), index) <= value(it).repLen)
    {
        // Interval contains only $-edges.
        if (++saRangeBegin >= saRangeEnd)
            return false;
    }

    childDesc.range.i1 = saRangeBegin;

    // Get first and last characters in interval.
    TAlphabet cLeft = textAt(posAdd(saAt(saRangeBegin, index), value(it).repLen), index);
    TAlphabet cRight = textAt(posAdd(saAt(saRangeEnd - 1, index), value(it).repLen), index);

#ifdef SEQAN_DEBUG
    std::cout << "cLeft: " << cLeft << std::endl;
    std::cout << "cRight: " << cRight << std::endl;
#endif

    // Fill child dir.
    // TODO(esiragusa): Remove workaround for alphabets with quality values
    for (typename ValueSize<TAlphabet>::Type c = ordValue(cLeft); c < ordValue(cRight); ++c)
    {
        TSAIterator saBegin = begin(sa, Standard()) + saRangeBegin;
        TSASize saLen = saRangeEnd - saRangeBegin;
        TSearchTreeIterator node(saBegin, saLen);
        TAlphabet edgeLabel = c;
        TSAIterator upperBound = _upperBoundSA(text, node, edgeLabel, value(it).repLen);

        childDesc.dir[ordValue(c)] = upperBound - begin(sa, Standard());
        saRangeBegin = childDesc.dir[ordValue(c)];
    }
//    childDesc.dir[ordValue(cRight)] = saRangeEnd;
    arrayFill(begin(childDesc.dir, Standard()) + ordValue(cRight), end(childDesc.dir, Standard()), saRangeEnd);

#ifdef SEQAN_DEBUG
    for (typename ValueSize<TAlphabet>::Type i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
        std::cout << value(childDesc.dir, i) << std::endl;
    std::cout << std::endl;
#endif

    // Update child repLen, lastChar and range.
    childDesc.repLen++;
    childDesc.lastChar = cLeft;
    childDesc.range.i2 = childDesc.dir[ordValue(childDesc.lastChar)];

    _historyPush(it);
    value(it) = childDesc;

    return true;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goRight(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it,
                     VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Value<TIndex>::Type                    TAlphabet;

#ifdef SEQAN_DEBUG
    std::cout << "goRight" << std::endl;
#endif

    if (isRoot(it))
        return false;

#ifdef SEQAN_DEBUG
    std::cout << length(value(it).dir) << std::endl;
    for (typename ValueSize<TAlphabet>::Type i = 0; i < ValueSize<TAlphabet>::VALUE; ++i)
        std::cout << value(value(it).dir, i) << std::endl;
    std::cout << std::endl;
#endif

    // TODO(esiragusa): Remove workaround for alphabets with quality values
    for (typename ValueSize<TAlphabet>::Type lastChar = ordValue(value(it).lastChar) + 1; lastChar < ValueSize<TAlphabet>::VALUE; ++lastChar)
    {
        value(it).range.i1 = value(it).range.i2;
        value(it).range.i2 = value(it).dir[ordValue(value(it).lastChar)];
        value(it).lastChar = lastChar;

        if (value(it).range.i1 < value(it).range.i2)
            return true;
    }

    return false;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TValue>
inline bool _getNodeByChar(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TSpec> > const & it,
                           TValue c, typename VertexDescriptor<Index<TText, IndexSa<TIndexSpec> > >::Type & childDesc)
{
    typedef Index<TText, IndexSa<TIndexSpec> >              TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type           TSA;
    typedef typename Size<TSA const>::Type                  TSASize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIterator;
    typedef SearchTreeIterator<TSA const, SortedList>       TSearchTreeIterator;
    typedef typename Value<TIndex>::Type                    TAlphabet;

    if (_isLeaf(it, EmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

    childDesc = value(it);
    resize(childDesc.dir, ValueSize<TAlphabet>::VALUE);
    arrayFill(begin(childDesc.dir, Standard()), end(childDesc.dir, Standard()), MaxValue<TSASize>::VALUE);

#ifdef SEQAN_DEBUG
    std::cout << "parent: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);

    Pair<TSAIterator> range = _equalRangeSA(text, node, c, value(it).repLen);

    if (range.i1 >= range.i2)
        return false;

    childDesc.range.i1 = range.i1 - begin(sa, Standard());
    childDesc.range.i2 = range.i2 - begin(sa, Standard());

#ifdef SEQAN_DEBUG
    std::cout << "child: " <<  childDesc.range.i1 << " " << childDesc.range.i2 << std::endl;
#endif

    childDesc.repLen++;

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

    if (_isLeaf(it, EmptyEdges()))
        return false;

    TIndex const & index = container(it);
    TSA const & sa = indexSA(index);
    TText const & text = indexText(index);

    typename VertexDescriptor<TIndex>::Type childDesc = value(it);

#ifdef SEQAN_DEBUG
    std::cout << "parent: " << value(it).range.i1 << " " << value(it).range.i2 << std::endl;
#endif

    TSAIterator saBegin = begin(sa, Standard()) + value(it).range.i1;
    TSASize saLen = isRoot(it) ? length(sa) : value(it).range.i2 - value(it).range.i1;
    TSearchTreeIterator node(saBegin, saLen);
    Pair<TSAIterator> range = _equalRangeSA(text, node, pattern, value(it).repLen);

    if (range.i1 >= range.i2)
        return false;

    childDesc.range.i1 = range.i1 - begin(sa, Standard());
    childDesc.range.i2 = range.i2 - begin(sa, Standard());

#ifdef SEQAN_DEBUG
    std::cout << "child: " <<  childDesc.range.i1 << " " << childDesc.range.i2 << std::endl;
#endif

    childDesc.repLen += length(pattern);

    _historyPush(it);
    value(it) = childDesc;

    lcp = childDesc.repLen;

    return true;
}

//template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
//inline bool _goDownString(Iter< Index<TText, IndexSa<TIndexSpec> >, VSTree< TopDown<TSpec> > > &node,
//                          TString const &pattern,
//                          TSize &lcp)
//{
//    typedef typename Iterator<TString const, Standard>::Type	PatternIterator;
//
//    lcp = 0;
//
//    PatternIterator patternBegin = begin(pattern, Standard());
//    PatternIterator patternEnd = end(pattern, Standard());
//
//    // TODO(esiragusa): Specialize to compute the whole hash at once
//    for (PatternIterator patternIt = patternBegin; patternIt != patternEnd; ++patternIt)
//    {
//        if (!_goDownChar(node, *patternIt)) return false;
//        ++lcp;
//    }
//
//    return true;
//}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool _goUp(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
#ifdef SEQAN_DEBUG
    std::cout << "goUp" << std::endl;
#endif

    if (!empty(it.history))
    {
        value(it).range = back(it.history).range;
        value(it).lastChar = back(it.history).lastChar;
        value(it).dir = back(it.history).dir;
        // TODO(esiragusa): Restore repLen from history
        value(it).repLen--;
        pop(it.history);
        return true;
    }
    return false;
}

// ============================================================================

template <typename TText, class TIndexSpec, class TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    it._parentDesc = value(it);
    value(it).parentRight = value(it).range.i2;
}

template <typename TText, class TIndexSpec, class TSpec>
inline void _historyPush(Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    typedef Iter<Index<TText, IndexSa<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > TIter;
    typename HistoryStackEntry_<TIter>::Type h;
    h.dir = value(it).dir;
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

    bool result = true;
    if ((!open(getFibre(index, EsaText()), toCString(name), openMode)) &&
        (!open(getFibre(index, EsaText()), fileName, openMode)))
        result = false;

    name = fileName;    append(name, ".sa");    open(getFibre(index, EsaSA()), toCString(name), openMode);
    return result;
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
        (!save(getFibre(index, EsaText()), fileName, openMode)))
        return false;

    name = fileName;    append(name, ".sa");    save(getFibre(index, EsaSA()), toCString(name), openMode);
    return true;
}

template <typename TObject, typename TSpec>
inline bool save(Index<TObject, IndexSa<TSpec> > & index, const char * fileName)
{
    return save(index, fileName, DefaultOpenMode<Index<TObject, IndexSa<TSpec> > >::VALUE);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INDEX_SA_STREE_H_
