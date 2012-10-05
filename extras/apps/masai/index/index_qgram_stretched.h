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

#ifndef SEQAN_EXTRAS_APPS_MASAI_INDEX_QGRAM_STRETCHED_H_
#define SEQAN_EXTRAS_APPS_MASAI_INDEX_QGRAM_STRETCHED_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(esiragusa): Rename this tag to NonRefined or Unsorted and spec QGramIndex
template <typename TSpec = void>
struct Stretched;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline bool isLeaf(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<Stretched<TSpec> > > > const & it)
{
    return (value(it).hash.i1 + 1 >= value(it).hash.i2) && (value(it).range.i1 + 1 >= value(it).range.i2);
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec>
inline typename Size<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> > >::Type
repLength(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<Stretched<TSpec> > > > const & it)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;

    if (isLeaf(it))
        return suffixLength(saAt(value(it).range.i1, container(it)), container(it));

    return value(it).repLen;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder>
inline bool _goDown(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<Stretched<TSpec> > > > & it,
                    VSTreeIteratorTraits<TDfsOrder, True> const)
{
    typedef Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >   TIndex;
    typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
    typedef typename Host<TShape>::Type                         TAlphabet;

    if (isLeaf(it))
        return false;

    // TODO(esiragusa): check nodeHullPredicate

    typename VertexDescriptor<TIndex>::Type nodeDesc;

    if (value(it).hash.i1 + 1 >= value(it).hash.i2)
    {
        nodeDesc = value(it);
        nodeDesc.range.i2 = nodeDesc.range.i1 + 1;
        nodeDesc.repLen++;
        _historyPush(it);
        value(it) = nodeDesc;
        return true;
    }

    // TODO(esiragusa): Remove workaround for alphabets with quality values
    for (unsigned char c = 0; c < +ValueSize<TAlphabet>::VALUE; ++c)
    {
        if (_getNodeByChar(it, c, nodeDesc))
        {
            _historyPush(it);
            value(it) = nodeDesc;
            return true;
        }
    }

    return false;
}

template <typename TText, typename TShapeSpec, typename TIndexSpec, typename TSpec, typename TDfsOrder, typename THideEmptyEdges>
inline bool _goRight(Iter<Index<TText, IndexQGram<TShapeSpec, TIndexSpec> >, VSTree<TopDown<Stretched<TSpec> > > > & it,
                     VSTreeIteratorTraits<TDfsOrder, THideEmptyEdges> const)
{
    if (isRoot(it))
        return false;

    if (isLeaf(it))
    {
        if (value(it).range.i2 >= nodeUp(it).range.i2)
            return false;

        value(it).range.i1 = value(it).range.i2;
        value(it).range.i2++;
        return true;
    }

    return _getNextNode(it);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_EXTRAS_APPS_MASAI_INDEX_QGRAM_STRETCHED_H_
