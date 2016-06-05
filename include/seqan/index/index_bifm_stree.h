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

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TDirection>
inline void update(Iter<Index<TText, BidirectionalIndex<FMIndex<TOccSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > & it, TDirection)
{
    typedef typename IfC<IsSameType<TDirection, Tag<BidirectionalFwd_> >::VALUE, Rev, Fwd>::Type TOppositeDirection;

    typedef typename IfC<IsSameType<TDirection, Tag<BidirectionalFwd_> >::VALUE, TText, typename RevTextFibre<TText>::Type>::Type TDirText;
    typedef typename IfC<IsSameType<TDirection, Tag<BidirectionalFwd_> >::VALUE, typename RevTextFibre<TText>::Type, TText>::Type TOppDirText;

    typedef Iter<Index<TDirText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > TDirIter;
    typedef Iter<Index<TOppDirText, FMIndex<TOccSpec, TIndexSpec> >, VSTree<TopDown<TSpec> > > TOppDirIter;

    _historyPush(_iter(it, TOppositeDirection()));

    TDirIter & dirIter = _iter(it, TDirection());
    TOppDirIter & oppDirIter = _iter(it, TOppositeDirection());

    value(oppDirIter).range.i1 += value(dirIter).smaller;
    value(oppDirIter).range.i2 = value(oppDirIter).range.i1 + value(dirIter).range.i2 - value(dirIter).range.i1;

    value(oppDirIter).repLen = value(dirIter).repLen; // do not increment in case of goRight
}


template <typename TText, typename TOccSpec, typename TIndexSpec, typename TSpec, typename TString, typename TSize, typename TDirection>
inline bool
_goDownString(Iter<Index<TText, BidirectionalIndex<FMIndex<TOccSpec, TIndexSpec> > >, VSTree<TopDown<TSpec> > > &it,
              TString const & string,
              TSize & lcp,
              TDirection)
{
    typedef Index<TText, FMIndex<TOccSpec, TIndexSpec> >        TIndex;
    typedef typename Size<TIndex>::Type                         TSize2;
    typedef Pair<TSize2>                                        TRange;
    typedef typename Iterator<TString const, Standard>::Type    TStringIter;

    TStringIter stringIt = begin(string, Standard());
    TStringIter stringEnd = end(string, Standard());

    _historyPush(_iter(it, TDirection()));

    for (lcp = 0; stringIt != stringEnd; ++stringIt, ++lcp)
    {
        TRange _range;
        TSize2 _smaller = 0;

        // NOTE(esiragusa): this should be faster only for texts over small alphabets consisting of few/long sequences.
        if (isLeaf(_iter(it, TDirection())) || !_getNodeByChar(_iter(it, TDirection()), value(_iter(it, TDirection())), _range, _smaller, value(stringIt))) break;

        value(_iter(it, TDirection())).range = _range;
        value(_iter(it, TDirection())).smaller = _smaller;
        update(it, TDirection());
    }

    value(_iter(it, TDirection())).repLen += lcp;

    if (lcp) value(_iter(it, TDirection())).lastChar = value(stringIt - 1);

    return stringIt == stringEnd;
}

}

#endif /* INDEX_BIFM_STREE_H_ */
