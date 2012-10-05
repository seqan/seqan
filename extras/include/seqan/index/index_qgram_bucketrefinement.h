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

#ifndef SEQAN_EXTRAS_INDEX_QGRAM_BUCKETREFINEMENT_H_
#define SEQAN_EXTRAS_INDEX_QGRAM_BUCKETREFINEMENT_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
 .Spec.BucketRefinement
 ..summary:An index based on a refined array of sorted q-grams.
 ..cat:Index
 ..general:Spec.IndexQGram
 ..signature:Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >
 ..param.TText:The text type.
 ...type:Class.String
 ..param.TShapeSpec:The @Class.Shape@ specialization type.
 ...note:This can be either a $TSpec$ argument (e.g. $SimpleShape$) or a complete @Class.Shape@ class (e.g. Shape<Dna, SimpleShape>).
 ..remarks:This index refines q-gram buckets by sorting and uses binary search to locate them.
 ..include:seqan/index.h
 */

struct BucketRefinement_;
typedef Tag<BucketRefinement_> BucketRefinement;
    
template <typename TText, typename TShapeSpec>
class Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >:
    public Index<TText, IndexQGram<TShapeSpec> >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBase;
    typedef Index<TText, IndexSa<InfixSegment> >     TIndexSa;

    TIndexSa    _indexSa;

    Index() :
        TBase(),
        _indexSa(indexText(*this))
    {
        _setHost(*this);
    }

    Index(Index & other) :
        TBase(static_cast<TBase &>(other)),
        _indexSa(other._indexSA)
    {
        _setHost(*this);
    }

    Index(Index const & other) :
        TBase(static_cast<TBase const &>(other)),
        _indexSa(other._indexSa)
    {
        _setHost(*this);
    }

    template <typename TText_>
    Index(TText_ & _text) :
        TBase(_text),
        _indexSa(indexText(*this))
    {
        _setHost(*this);
    }

    template <typename TText_>
    Index(TText_ const & _text) :
        TBase(_text),
        _indexSa(indexText(*this))
    {
        _setHost(*this);
    }

    template <typename TText_, typename TShape_>
    Index(TText_ & _text, TShape_ const & _shape) :
        TBase(_text, _shape),
        _indexSa(indexText(*this))
    {
        _setHost(*this);
    }

    template <typename TText_, typename TShape_>
    Index(TText_ const & _text, TShape_ const & _shape) :
        TBase(_text, _shape),
        _indexSa(indexText(*this))
    {
        _setHost(*this);
    }

};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
class Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >         TIndex;

    typedef Index<TText, IndexQGram<TShapeSpec> >                           TBase;
    typedef typename Iterator<TBase, TopDown<TSpec> >::Type                 TBaseIterator;

    typedef Index<TText, IndexSa<InfixSegment> >                            TIndexSa;
    typedef typename Iterator<TIndexSa, TopDown<TSpec> >::Type              TIndexSaIterator;

    TBaseIterator       _topIterator;
    TIndexSaIterator    _bottomIterator;

// ============================================================================

    Iter() {}

    Iter(TIndex & _index) :
        _topIterator(),
        _bottomIterator(_index._indexSa)
    {
        _indexRequireTopDownIteration(_index);
        _topIterator = TBaseIterator(static_cast<TBase &>(_index));
        goRoot(_topIterator);
        goRoot(_bottomIterator);
    }

    Iter(Iter const & _origin) :
        _topIterator(_origin._topIterator),
        _bottomIterator(_origin._bottomIterator)
    {}

// ============================================================================

    inline Iter const &
    operator=(Iter const & _origin)
    {
        _bottomIterator = _origin._bottomIterator;
        _topIterator = _origin._topIterator;
        return *this;
    }

};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
class Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >                 TIndex;

    typedef Index<TText, IndexQGram<TShapeSpec> >                                   TBase;
    typedef typename Iterator<TBase, TopDown<ParentLinks<TSpec> > >::Type          TBaseIterator;

    typedef Index<TText, IndexSa<InfixSegment> >                                    TIndexSa;
    typedef typename Iterator<TIndexSa, TopDown<ParentLinks<TSpec> > >::Type       TIndexSaIterator;

    TBaseIterator       _topIterator;
    TIndexSaIterator    _bottomIterator;

// ============================================================================

    Iter() {}

    Iter(TIndex & _index) :
        _topIterator(),
        _bottomIterator(_index._indexSa)
    {
        _indexRequireTopDownIteration(_index);
        _topIterator = TBaseIterator(static_cast<TBase &>(_index));
        goRoot(_topIterator);
        goRoot(_bottomIterator);
    }

    Iter(Iter const & _origin) :
        _topIterator(_origin._topIterator),
        _bottomIterator(_origin._bottomIterator)
    {}

// ============================================================================

    inline Iter const &
    operator=(Iter const & _origin)
    {
        _bottomIterator = _origin._bottomIterator;
        _topIterator = _origin._topIterator;
        return *this;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreText>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreText>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreRawText>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreRawText>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreSA>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreRawSA>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreRawSA>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreDir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreDir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSADir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreSADir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreShape>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreShape>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreCounts>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreCounts>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreCountsDir>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreCountsDir>
{};

template <typename TText, typename TShapeSpec>
struct Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreBucketMap>:
    public Fibre<Index<TText, IndexQGram<TShapeSpec> >, FibreBucketMap>
{};

// ============================================================================

template <typename TText, typename TShapeSpec>
struct DefaultIndexCreator<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>:
    public DefaultIndexCreator<Index<TText, IndexSa<InfixSegment> >, FibreSA>
{};

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
struct Iterator<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >
{
    typedef Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >     Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TShapeSpec>
inline bool indexCreate(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, FibreSADir, Default const)
{
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >                 TIndex;
    typedef Index<TText, IndexQGram<TShapeSpec> >                                   TBase;

//    indexCreate(static_cast<TBase &>(index), FibreSADir(), Default());
//    _refineQGramIndex(indexSA(index), indexDir(index), indexText(index), weight(indexShape(index)), lengthSum(indexText(index)));

    // Create QGram directory.
    resize(indexDir(index), _fullDirLength(index), Exact());
    createQGramIndexDirOnly(indexDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));

    // Create full SA.
    indexCreate(index, FibreSA());

    // Remove too short suffixes from SA.
    _pruneSA(index);

    // Update indexSA host.
    _setHost(index);

    return true;
}

template <typename TText, typename TShapeSpec>
void _setHost(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index)
{
    // Make SA index fibre point to QGramSA index fibre.
    setHost(indexSA(index._indexSa), indexSA(index));
    setBeginPosition(indexSA(index._indexSa), 0);
    setEndPosition(indexSA(index._indexSa), length(indexSA(index)));
}

template <typename TText, typename TShapeSpec>
void _pruneSA(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index)
{
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >     TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA, Standard>::Type                      TSAIterator;

    TSA & sa = indexSA(index);

    TSAIterator saBegin = begin(sa, Standard());
    TSAIterator saEnd = end(sa, Standard());

    TSAIterator saOld = saBegin;
    TSAIterator saNew = saOld;

    while (saOld != saEnd)
    {
        if (suffixLength(*saOld, index) < weight(indexShape(index)))
            ++saOld;
        else
        {
            *saNew = *saOld;
            ++saOld;
            ++saNew;
        }
    }

    resize(sa, saNew - saBegin, Exact());
}

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool _atTop(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return isRoot(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline void _implantSa(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    typedef Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >     TIndex;
    typedef typename Value<TIndex>::Type                                TAlphabet;

    _historyPush(it._bottomIterator);

    resize(value(it._bottomIterator).dir, ValueSize<TAlphabet>::VALUE);
    arrayFill(begin(value(it._bottomIterator).dir, Standard()),
              end(value(it._bottomIterator).dir, Standard()),
              value(it._topIterator).range.i1);

    value(it._bottomIterator).repLen = value(it._topIterator).repLen;
    value(it._bottomIterator).range = value(it._topIterator).range;
    value(it._bottomIterator).lastChar = value(it._topIterator).lastChar;
}

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
inline Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > const &
container(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    SEQAN_ASSERT(_atTop(it));
    return container(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > &
container(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    SEQAN_ASSERT(_atTop(it));
    return container(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename VertexDescriptor<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type &
value(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > >&it)
{
    SEQAN_ASSERT(_atTop(it));
    return value(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename VertexDescriptor<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type const &
value(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    SEQAN_ASSERT(_atTop(it));
    return value(it._topIterator);
}

// ============================================================================

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool isRoot(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return isRoot(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool isLeaf(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    // TODO(esiragusa): Fix isLeaf: last qgram is a top leaf but there is no subtree attached
    return isLeaf(it._topIterator) && isLeaf(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline void goRoot(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    goRoot(it._bottomIterator);
    goRoot(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool goDown(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    if (_atTop(it))
    {
        if (goDown(it._topIterator))
            return true;

        _implantSa(it);
    }
    return goDown(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TObject>
inline bool goDown(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it,
                   TObject const & obj)
{
    unsigned lcp = 0;
    return goDown(it, obj, lcp);
}

template <typename TText, typename TShapeSpec, typename TSpec, typename TString, typename TSize>
inline bool goDown(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it,
                   TString const & pattern, TSize & lcp)
{
    if (_atTop(it))
    {
        if (goDown(it._topIterator, pattern, lcp))
            return true;

        if (!isLeaf(it._topIterator))
            return false;

        _implantSa(it);
    }

    return goDown(it._bottomIterator, suffix(pattern, lcp));
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool goRight(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > & it)
{
    return _atTop(it) ? goRight(it._topIterator) : goRight(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline bool goUp(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    if (!_atTop(it))
    {
        goUp(it._bottomIterator);
        if (!isRoot(it._bottomIterator))
            return true;
    }

    return goUp(it._topIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
countOccurrences(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? countOccurrences(it._topIterator) : countOccurrences(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename SAValue<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
getOccurrence(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? getOccurrence(it._topIterator) : getOccurrence(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreSA>::Type const>::Type
getOccurrences(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? getOccurrences(it._topIterator) : getOccurrences(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type
repLength(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? repLength(it._topIterator) : repLength(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline typename Infix<typename Fibre<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, FibreText>::Type const>::Type
representative(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? representative(it._topIterator) : representative(it._bottomIterator);
}

template <typename TText, typename TShapeSpec, typename TSpec>
inline Pair<typename Size<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > >::Type>
range(Iter<Index<TText, IndexQGram<TShapeSpec, BucketRefinement> >, VSTree<TopDown<TSpec> > > const & it)
{
    return _atTop(it) ? range(it._topIterator) : range(it._bottomIterator);
}

template <typename TText, typename TShapeSpec>
inline bool open(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, const char * fileName)
{
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBase;

    if (open(static_cast<TBase &>(index), fileName))
    {
        _setHost(index);
        return true;
    }

    return false;
}

template <typename TText, typename TShapeSpec>
inline bool save(Index<TText, IndexQGram<TShapeSpec, BucketRefinement> > & index, const char * fileName)
{
    typedef Index<TText, IndexQGram<TShapeSpec> >    TBase;

    return save(static_cast<TBase &>(index), fileName);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_INDEX_QGRAM_BUCKETREFINEMENT_H_
