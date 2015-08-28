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

#ifndef INDEX_BIFM_STREE_H_
#define INDEX_BIFM_STREE_H_

namespace seqan {

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Class BidirectionalFMIndex-Iter
// ----------------------------------------------------------------------------

// TopDown-Iterator for bidirectional FM index
template <typename TText, typename TOccSpec, typename TLengthSum, typename TBidirectional, typename TSpec>
class Iter<Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, TBidirectional> > >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, TBidirectional> > >  TBiIndex;

    typedef typename RevTextFibre<TText>::Type                                                      TRevText;

    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >     TFwdIndex;
    typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >  TRevIndex;

    typedef Iter<TFwdIndex, VSTree<TopDown<TSpec> > >                                                   TFwdIndexIter;
    typedef Iter<TRevIndex, VSTree<TopDown<TSpec> > >                                                   TRevIndexIter;

    TFwdIndexIter    fwdIter;
    TRevIndexIter    bwdIter;

//____________________________________________________________________________

    Iter():
        fwdIter(),
        bwdIter() {}

    Iter(TBiIndex &_index):
        fwdIter(*(&_index.fwd)),
        bwdIter(*(&_index.rev))
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

    Iter(TBiIndex &_index, MinimalCtor):
        fwdIter(*(&_index.fwd), MinimalCtor()),
        bwdIter(*(&_index.rev), MinimalCtor())
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

    template <typename TSpec2>
    Iter(Iter<TBiIndex, VSTree<TopDown<TSpec2> > > const &_origin):
        fwdIter(*(&_origin.fwdIter)),
        bwdIter(*(&_origin.bwdIter))
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

//____________________________________________________________________________

    template <typename TSpec2>
    inline Iter const &
    operator = (Iter<TBiIndex, VSTree<TopDown<TSpec2> > > const &_origin)
    {
        fwdIter = *(&_origin.fwdIter);
        bwdIter = *(&_origin.bwdIter);
        return *this;
    }
};

// TopDown-Iterator with history stack for bidirectional FM index
template <typename TText, typename TOccSpec, typename TLengthSum, typename TBidirectional, typename TSpec>
class Iter<Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, TBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >
{
public:
    typedef Index<TText, BidirectionalFMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, TBidirectional> > >  TBiIndex;

    typedef typename RevTextFibre<TText>::Type                                                            TRevText;

    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >           TFwdIndex;
    typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >        TRevIndex;

    typedef Iter<TFwdIndex, VSTree<TopDown<ParentLinks<TSpec> > > >                                           TFwdIndexIter;
    typedef Iter<TRevIndex, VSTree<TopDown<ParentLinks<TSpec> > > >                                           TRevIndexIter;

    TFwdIndexIter    fwdIter;
    TRevIndexIter    bwdIter;

//____________________________________________________________________________

    Iter():
        fwdIter(),
        bwdIter() {}

    Iter(TBiIndex &_index):
        fwdIter(*(&_index.fwd)),
        bwdIter(*(&_index.rev))
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

    Iter(TBiIndex &_index, MinimalCtor):
        fwdIter(*(&_index.fwd), MinimalCtor()),
        bwdIter(*(&_index.rev), MinimalCtor())
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

    template <typename TSpec2>
    Iter(Iter<TBiIndex, VSTree<TopDown<ParentLinks<TSpec2> > > > const &_origin):
        fwdIter(*(&_origin.fwdIter)),
        bwdIter(*(&_origin.bwdIter))
    {
        fwdIter.setRevIter(bwdIter);
        bwdIter.setRevIter(fwdIter);
    }

//____________________________________________________________________________

    template <typename TSpec2>
    inline Iter const &
    operator = (Iter<TBiIndex, VSTree<TopDown<ParentLinks<TSpec2> > > > const &_origin)
    {
        fwdIter = *(&_origin.fwdIter);
        bwdIter = *(&_origin.bwdIter);
        return *this;
    }
};

// Specialization of undirectional FM index iterator for use in the bidirectional FM index iterator
template <typename TText, typename TOccSpec, typename TLengthSum, typename TSpec>
class Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<TSpec> > >
{
public:
    typedef Iter    iterator;

    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >       TFwdIndex;
    typedef typename RevTextFibre<TText>::Type                                                        TRevText;
    typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >    TRevIndex;

    typedef Iter<TFwdIndex, VSTree< TopDown<TSpec> > >    TFwdIndexIter;
    typedef Iter<TRevIndex, VSTree< TopDown<TSpec> > >    TRevIndexIter;

    typedef typename VertexDescriptor<TFwdIndex>::Type    TVertexDesc;

    TRevIndexIter     *revIter;    // container of all necessary tables

    TFwdIndex const    *index;        // container of all necessary tables
    TVertexDesc        vDesc;        // current interval in suffix array and
                                // right border of parent interval (needed in goRight)

    // pseudo history stack (to go up at most one node)
    TVertexDesc        _parentDesc;

//____________________________________________________________________________

    Iter() : index() {}

    Iter(TFwdIndex &_index):
        revIter(0),
        index(&_index)
    {
        _indexRequireTopDownIteration(_index);
        goRoot(*this);
    }

    Iter(TFwdIndex &_index, MinimalCtor):
        index(&_index),
        vDesc(MinimalCtor()),
        _parentDesc(MinimalCtor()) {}

    // NOTE(esiragusa): _parentDesc is unitialized
    Iter(TFwdIndex &_index, TVertexDesc const &_vDesc):
        index(&_index),
        vDesc(_vDesc)
    {
        _indexRequireTopDownIteration(_index);
    }

    template <typename TSpec2>
    Iter(Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin):
        index(&container(_origin)),
        vDesc(value(_origin)),
        _parentDesc(nodeUp(_origin)) {}

//____________________________________________________________________________

    template <typename TSpec2>
    inline Iter const &
    operator = (Iter<TFwdIndex, VSTree<TopDown<TSpec2> > > const &_origin)
    {
        revIter = 0;
        index = &container(_origin);
        vDesc = value(_origin);
        _parentDesc = nodeUp(_origin);
        return *this;
    }

    void setRevIter(TRevIndexIter &_revIter)
    {
        revIter = &_revIter;
    }
};

// Specialization of undirectional FM index iterator with history stack for use in the bidirectional FM index iterator with history stack
template <typename TText, typename TOccSpec, typename TLengthSum, typename TSpec>
class Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > >:
    public Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<TSpec> > >
{
public:

    typedef typename RevTextFibre<TText>::Type TRevText;
    typedef Index<TRevText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > TRevIndex;
    typedef Iter<TRevIndex, VSTree< TopDown<ParentLinks<TSpec> > > > TRevIndexIter;

    TRevIndexIter     *revIter;    // container of all necessary tables

    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > TIndex;
    typedef Iter<TIndex, VSTree<TopDown<TSpec> > > TBase;
    typedef typename HistoryStack_<Iter>::Type     TStack;
    typedef Iter                                   iterator;

    TStack            history;    // contains all previously visited intervals (allows to go up)

//____________________________________________________________________________

    SEQAN_HOST_DEVICE
    Iter() :
        TBase()
    {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index):
        revIter(0),
        TBase(_index) {}

    SEQAN_HOST_DEVICE
    Iter(TIndex &_index, MinimalCtor):
        TBase(_index, MinimalCtor()) {}

    SEQAN_HOST_DEVICE
    Iter(Iter const &_origin):
        TBase((TBase const &)_origin),
        history(_origin.history) {}

//____________________________________________________________________________

    SEQAN_HOST_DEVICE inline
    Iter const &
    operator = (Iter const &_origin)
    {
        revIter = 0;
        *(TBase*)(this) = _origin;
        history = _origin.history;
        return *this;
    }

    void setRevIter(TRevIndexIter &_revIter)
    {
        revIter = &_revIter;
    }
};

// ============================================================================
// Functions
// ============================================================================

// specialized method for updating both search intervals of the fm iterators
template <typename TText, typename TOccSpec, typename TLengthSum, typename TSpec, typename TChar>
inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<TSpec> > > const & it,
                           typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type const & vDesc,
                           Pair<typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type> & _range,
                           TChar c)
{
    typedef typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type TSize;
    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type TLF;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    _range = range(index, vDesc);
    TSize smaller1 = 0, smaller2 = 0;
    _range.i1 = lf(_range.i1, c, smaller1);
    _range.i2 = lf(_range.i2, c, smaller2);

    if (_range.i1 < _range.i2)
    {
        if (_isRoot(vDesc))
        {
            value(*it.revIter).range.i1 = _range.i1;
            value(*it.revIter).range.i2 = _range.i2;
        }
        else
        {
            value(*it.revIter).range.i1 += smaller2 - smaller1;
            value(*it.revIter).range.i2 = value(*it.revIter).range.i1 + (_range.i2 - _range.i1);
        }

        return true;
    }

    return false;
}
template <typename TText, typename TOccSpec, typename TLengthSum, typename TSpec, typename TChar>
inline bool _getNodeByChar(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > >, VSTree<TopDown<TSpec> > > const & it,
                           typename VertexDescriptor<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type const & vDesc,
                           Pair<typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type> & _range1,
                           Pair<typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type> & _range2,
                           TChar c)
{
    typedef typename Size<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > >::Type TSize;
    typedef Index<TText, FMIndex<TOccSpec, FMIndexConfig<TOccSpec, TLengthSum, FMBidirectional> > > TIndex;
    typedef typename Fibre<TIndex, FibreLF>::Type TLF;

    TIndex const & index = container(it);
    TLF const & lf = indexLF(index);

    _range1 = range(index, vDesc);
    TSize smaller1 = 0, smaller2 = 0;
    _range1.i1 = lf(_range1.i1, c, smaller1);
    _range1.i2 = lf(_range1.i2, c, smaller2);

    if (_range1.i1 < _range1.i2)
    {
        if (_isRoot(vDesc))
        {
            _range2.i1 = _range1.i1;
            _range2.i2 = _range1.i2;
        }
        else
        {
            _range2.i1 = value(*it.revIter).range.i1 + smaller2 - smaller1;
            _range2.i2 = _range2.i1 + _range1.i2 - _range1.i1;
        }

        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function goUp()                                                   [Iterator]
// ----------------------------------------------------------------------------

// go up one edge
template <typename TText, typename TOccSpec, typename TSpec2, typename TLengthSum, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_goUp(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TSpec2, TLengthSum, FMBidirectional> > >, VSTree<TopDown<TSpec> > > & it)
{
    if (isRoot(it)) return false;

    _historyPop(it);
    _historyPop(*it.revIter);
    return true;
}

template <typename TText, typename TOccSpec, typename TSpec2, typename TLengthSum, typename TSpec>
SEQAN_HOST_DEVICE inline bool
_goUp(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TSpec2, TLengthSum, FMBidirectional> > >, VSTree<TopDown<ParentLinks<TSpec> > > > & it)
{
    if (isRoot(it)) return false;

    _historyPop(it);
    _historyPop(*it.revIter);
    return true;
}

// ----------------------------------------------------------------------------
// Function goRoot()                                                 [Iterator]
// ----------------------------------------------------------------------------

template < typename TText, typename TOccSpec, typename TSpec2, typename TLengthSum, typename TSpec>
SEQAN_HOST_DEVICE inline void goRoot(Iter<Index<TText, FMIndex<TOccSpec, FMIndexConfig<TSpec2, TLengthSum, FMBidirectional> > >, VSTree<TSpec> > &it)
{
    // avoid trying to access it.revIter during construction of it, when it.revIter is not set
    if (it.revIter)
    {
        _historyClear(*it.revIter);
        clear(*it.revIter);
        if (!empty(indexSA(container(*it.revIter))))
            _setSizeInval(value(*it.revIter).range.i2);
    }

    _historyClear(it);
    clear(it);                            // start in root node with range (0,infty)
    if (!empty(indexSA(container(it))))
        _setSizeInval(value(it).range.i2);    // infty is equivalent to length(index) and faster to compare
}

}

#endif /* INDEX_BIFM_STREE_H_ */
