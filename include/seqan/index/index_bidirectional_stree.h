// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Author: Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_BIDIRECTIONAL_STREE_H_
#define INDEX_BIDIRECTIONAL_STREE_H_

namespace seqan {

struct BidirectionalDirection {};

struct BidirectionalFwd_ : BidirectionalDirection {};
struct BidirectionalRev_ : BidirectionalDirection {};
typedef Tag<BidirectionalFwd_> Fwd;
typedef Tag<BidirectionalRev_> Rev;

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class BidirectionalIndex-Iter
// ----------------------------------------------------------------------------

// TopDown-Iterator for bidirectional indices
template <typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Index<TText, BidirectionalIndex<TIndexSpec> > TBiIndex;

    typedef typename RevTextFibre<TText>::Type          TRevText;
    typedef Index<TText, TIndexSpec>                    TFwdIndex;
    typedef Index<TRevText, TIndexSpec>                 TRevIndex;
    typedef Iter<TFwdIndex, VSTree<TopDown<TSpec> > >   TFwdIndexIter;
    typedef Iter<TRevIndex, VSTree<TopDown<TSpec> > >   TRevIndexIter;

    TFwdIndexIter    fwdIter;
    TRevIndexIter    revIter;

//____________________________________________________________________________

    Iter():
        fwdIter(),
        revIter() {}

    Iter(TBiIndex &_index):
        fwdIter(_index.fwd),
        revIter(_index.rev)
    {}

    template <typename TSpec2>
    Iter(Iter<TBiIndex, VSTree<TopDown<TSpec2> > > const &_origin):
        fwdIter(_origin.fwdIter),
        revIter(_origin.revIter)
    {}
};

// TopDown-Iterator with history stack for bidirectional indices
template <typename TText, typename TIndexSpec, typename TSpec>
class Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:
    typedef Index<TText, BidirectionalIndex<TIndexSpec> >  TBiIndex;

    typedef typename RevTextFibre<TText>::Type                      TRevText;
    typedef Index<TText, TIndexSpec>                                TFwdIndex;
    typedef Index<TRevText, TIndexSpec>                             TRevIndex;
    typedef Iter<TFwdIndex, VSTree<TopDown<ParentLinks<TSpec> > > > TFwdIndexIter;
    typedef Iter<TRevIndex, VSTree<TopDown<ParentLinks<TSpec> > > > TRevIndexIter;

    TFwdIndexIter    fwdIter;
    TRevIndexIter    revIter;

//____________________________________________________________________________

    Iter():
        fwdIter(),
        revIter() {}

    Iter(TBiIndex &_index):
        fwdIter(_index.fwd),
        revIter(_index.rev)
    {}

    template <typename TSpec2>
    Iter(Iter<TBiIndex, VSTree<TopDown<ParentLinks<TSpec2> > > > const &_origin):
        fwdIter(_origin.fwdIter),
        revIter(_origin.revIter)
    {}
};

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TIndexSpec, typename TSpec>
inline Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & _iter(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, Fwd)
{
    return it.fwdIter;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > const & _iter(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Fwd)
{
    return it.fwdIter;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline Iter<Index<typename RevTextFibre<TText>::Type, TIndexSpec>, VSTree<TopDown<TSpec> > > & _iter(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, Rev)
{
    return it.revIter;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline Iter<Index<typename RevTextFibre<TText>::Type, TIndexSpec>, VSTree<TopDown<TSpec> > > const & _iter(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Rev)
{
    return it.revIter;
}

// ----------------------------------------------------------------------------
// Function goDown() directional interface for unidirectional indexes[Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & it,
       Fwd const &)
{
    return goDown(it);
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & it,
       TObject const & obj,
       Fwd const &)
{
    return goDown(it, obj);
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & it,
       TString const & pattern,
       TSize & lcp,
       Fwd const &)
{
    return _goDownString(it, pattern, lcp);
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > &,
       Rev const &)
{
    SEQAN_ASSERT_MSG(false, "ERROR: Cannot goDown(it, Rev) on uni-directional index");
    return false;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & it,
       TObject const &,
       Rev const &)
{
    return goDown(it, Rev()); // fail above
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline bool
goDown(Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > & it,
       TString const &,
       TSize &,
       Rev const &)
{
    return goDown(it, Rev()); // fail above
}

// ----------------------------------------------------------------------------
// Function goDown()                                                 [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool goDown(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, Fwd const &)
{
    if (goDown(_iter(it, Fwd())))
    {
        update(it, Fwd());
        return true;
    }
    return false;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool goDown(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it, Rev const &)
{
    if (goDown(_iter(it, Rev())))
    {
        update(it, Rev());
        return true;
    }
    return false;
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool goDown(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > & it)
{
    return goDown(it, Fwd());
}

//------------------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject, typename TDirection>
inline bool
_goDownObject(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    TObject const &obj,
    False,
    TDirection)
{
    if (_goDownChar(_iter(it, TDirection()), obj))
    {
        update(it, TDirection());
        return true;
    }
    return false;
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject, typename TDirection>
inline bool
_goDownObject(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    TObject const &obj,
    True,
    TDirection)
{
    typename Size<TIndexSpec>::Type dummy;
    return _goDownString(it, obj, dummy, TDirection());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject>
inline std::enable_if_t<!(std::is_same<TObject, Fwd>::value || std::is_same<TObject, Rev>::value),  bool>
goDown(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    TObject const &obj)
{
    return goDown(it, obj, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject>
inline bool
goDown(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    TObject const &obj,
    Rev const &)
{
    return _goDownObject(it, obj, typename IsSequence<TObject>::Type(), Rev());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TObject>
inline bool
goDown(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it,
    TObject const &obj,
    Fwd const &)
{
    return _goDownObject(it, obj, typename IsSequence<TObject>::Type(), Fwd());
}

//------------------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize>
inline std::enable_if_t<std::is_integral<TSize>::value,  bool>
goDown(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
    TString const &pattern,
    TSize const & lcp)
{
    return _goDownString(it, pattern, lcp, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TString, typename TSize, typename TDirection>
inline std::enable_if_t<std::is_integral<TSize>::value,  bool>
goDown(
    Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree< TopDown<TSpec> > > &it,
    TString const &pattern,
    TSize &lcp,
    Tag<TDirection>)
{
    return _goDownString(it, pattern, lcp, Tag<TDirection>());
}

// ----------------------------------------------------------------------------
// Function repLength()                                              [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, BidirectionalIndex<TIndexSpec> > >::Type
repLength(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const &it)
{
    return repLength(container(it.fwdIter), value(it.fwdIter));
}

// ----------------------------------------------------------------------------
// Function goRight()                                                [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool goRight(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it)
{
    return goRight(it, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDirection>
inline bool goRight(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > &it, Tag<TDirection>)
{
    if (goRight(_iter(it, Tag<TDirection>())))
    {
        update(it, Tag<TDirection>());
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function parentEdgeLabel()                                        [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename EdgeLabel<Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return parentEdgeLabel(it, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDirection>
inline typename EdgeLabel<Iter<Index<TText, TIndexSpec>, VSTree<TopDown<TSpec> > > >::Type
parentEdgeLabel(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Tag<TDirection>)
{
    return parentEdgeLabel(_iter(it, Tag<TDirection>()));
}

// ----------------------------------------------------------------------------
// Function getOccurrences()                                         [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Infix< typename Fibre<Index<TText, TIndexSpec>, FibreSA>::Type const >::Type
getOccurrences(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return getOccurrences(it, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Infix< typename Fibre<Index<TText, TIndexSpec>, FibreSA>::Type const >::Type
getOccurrences(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Fwd)
{
    return getOccurrences(_iter(it, Fwd()));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline typename Infix< typename Fibre<Index<typename RevTextFibre<TText>::Type, TIndexSpec>, FibreSA>::Type const >::Type
getOccurrences(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Rev)
{
    return getOccurrences(_iter(it, Rev()));
}

// ----------------------------------------------------------------------------
// Function isLeaf()                                                 [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool isLeaf(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return isLeaf(it, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDirection>
inline bool isLeaf(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Tag<TDirection>)
{
    return isLeaf(_iter(it, Tag<TDirection>()));
}

// ----------------------------------------------------------------------------
// Function isRoot()                                                 [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec>
inline bool isRoot(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return isRoot(it, Fwd());
}

template <typename TText, typename TIndexSpec, typename TSpec, typename TDirection>
inline bool isRoot(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it, Tag<TDirection>)
{
    return isRoot(_iter(it, Tag<TDirection>()));
}

// ----------------------------------------------------------------------------
// Function goUp()                                                   [Iterator]
// ----------------------------------------------------------------------------

// go up one edge
template <typename TText, typename TIndexSpec, typename TSpec>
inline bool goUp(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<ParentLinks<TSpec> > > > &it)
{
    return goUp(it.fwdIter) && goUp(it.revIter);
}

// ----------------------------------------------------------------------------
// Function goRoot()                                                 [Iterator]
// ----------------------------------------------------------------------------

template < typename TText, typename TIndexSpec, typename TSpec>
inline void goRoot(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TSpec> > &it)
{
    goRoot(it.fwdIter);
    goRoot(it.revIter);
}

// ----------------------------------------------------------------------------
// Function countOccurrences()                                       [Iterator]
// ----------------------------------------------------------------------------

template < typename TText, typename TIndexSpec, typename TSpec >
inline typename Size<Index<TText, BidirectionalIndex<TIndexSpec> > >::Type
countOccurrences(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TSpec> > const &it)
{
    return countOccurrences(it.fwdIter);
}

// ----------------------------------------------------------------------------
// Function representative()                                         [Iterator]
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TSpec, typename TDirection>
inline auto
representative(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it,
               Tag<TDirection> const &)
{
    return representative(_iter(it, Tag<TDirection>()));
}

template <typename TText, typename TIndexSpec, typename TSpec>
inline auto
representative(Iter<Index<TText, BidirectionalIndex<TIndexSpec> >, VSTree<TopDown<TSpec> > > const & it)
{
    return representative(it, Fwd());
}

}

#endif /* INDEX_BIDIRECTIONAL_STREE_H_ */
