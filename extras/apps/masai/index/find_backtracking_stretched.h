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
// Approximate string matching via backtracking on VSTrees
// ==========================================================================

#ifndef SEQAN_EXTRAS_APPS_MASAI_FIND_BACKTRACKING_STRETCHED_H_
#define SEQAN_EXTRAS_APPS_MASAI_FIND_BACKTRACKING_STRETCHED_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance>
class Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, Stretched<> > >
{
protected:
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<Stretched<Preorder> > >::Type    TIndexIterator;
    typedef String<TIndexIterator, Block<> >                            TParentStack;

    typedef typename Fibre<TIndex, FibreSA>::Type                       TSA;
    typedef typename Iterator<TSA const, Standard>::Type                TIterator;
    typedef typename Size<TIndex>::Type                                 TSize;

    typedef typename Fibre<TIndex, EsaText>::Type                       TSAText;
    typedef typename Infix<TSAText const>::Type                         TPrefix;
    typedef PrefixAligner_<TPrefix, TDistance>                          TPrefixAligner;

    typedef typename BacktrackingState_<TPrefix, TDistance>::Type       TState;
    typedef String<TState, Block<> >                                    TStateStack;

    typedef String<bool, Block<> >                                      TEndStack;

public:
    bool                index_iterator_at_root;

    Holder<TIndex>      data_host;
    TIndexIterator      index_iterator;
    TParentStack        index_parents;
    Pair<TSize>         index_range;

    TIterator           data_iterator;
    TSize               data_length;
    Pair<TIterator>     range;

    TPrefixAligner      prefix_aligner;
    TStateStack         state;
    TEndStack           atEnd;

    unsigned            exact;
    bool                search;

    // TODO(esiragusa): Remove depth, isLeaf(it) should return true when repLength(it) >= depth
    TSize               depth;

    Pattern() {}

    Pattern(TIndex & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    Pattern(TIndex const & index, TSize depth) :
        data_host(index),
        index_iterator(index),
        depth(depth)
    {
        setHost(*this, index);
    }

    ~Pattern()
    {
        // Empty stacks
        clear(this->index_parents);
        clear(this->state);
        clear(this->atEnd);
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TNeedle, typename TSpec, typename TDistance>
void _moveIteratorAtRoot(Pattern<Index<TNeedle, TSpec>, Backtracking<TDistance, Stretched<> > > & me)
{
    typedef Index<TNeedle, TSpec>                                       TIndex;
    typedef typename Iterator<TIndex, TopDown<Stretched<Preorder> > >::Type   TIndexIterator;
    me.index_iterator = TIndexIterator(host(me));
}

}

#endif  // #ifndef SEQAN_EXTRAS_APPS_MASAI_FIND_BACKTRACKING_STRETCHED_H_
