// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_
#define TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

template <typename TIter>
void waveletTreeStructureIteratorBegin(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    SEQAN_ASSERT_EQ(getCharacter(begin(waveletTreeStructure, typename Spec<TIter>::Type())), 'G');
}

// NOTE(esiragusa): IteratorContainer test should not use operator==()
//template <typename TIter>
//void waveletTreeStructureIteratorContainer(TIter & /*tag*/)
//{
//    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
//    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
//    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
//    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
//
//    String<TChar> text = "ACGTNACGTNACGTN";
//    TRightArrayBinaryTree waveletTreeStructure(text);
//
//    TIter it(waveletTreeStructure, 0);
//
//    SEQAN_ASSERT(container(it) == waveletTreeStructure);
//}

template <typename TIter>
void waveletTreeStructureIteratorEnd(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it = end(waveletTreeStructure, typename Spec<TIter>::Type());
    goToPosition(it, getPosition(it) - 1);

    SEQAN_ASSERT_EQ(getCharacter(it), 'N');
}

template <typename TIter>
void waveletTreeStructureIteratorGetCharacter(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    {
        TIter it(waveletTreeStructure, 0);
        SEQAN_ASSERT_EQ(getCharacter(it), 'G');
    }
    {
        TIter it(waveletTreeStructure, 1);
        SEQAN_ASSERT_EQ(getCharacter(it), 'C');
    }
    {
        TIter it(waveletTreeStructure, 2);
        SEQAN_ASSERT_EQ(getCharacter(it), 'T');
    }
    {
        TIter it(waveletTreeStructure, 3);
        SEQAN_ASSERT_EQ(getCharacter(it), 'N');
    }
}

template <typename TIter>
void waveletTreeStructureIteratorGetChildPos(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it = begin(waveletTreeStructure, typename Spec<TIter>::Type());

    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 2u);

    goLeftChild(it);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 0u);

    goRight(it);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 3u);
}

template <typename TIter>
void waveletTreeStructureIteratorGetNumChildVertices(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('G', 4u));
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('C', 0u));
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('T', 1u));
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('N', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(getSubTreeSize(it), 3u);

    goDown(it);
    SEQAN_ASSERT_EQ(getSubTreeSize(it), 0u);

    goRight(it);
    SEQAN_ASSERT_EQ(getSubTreeSize(it), 1u);

    goDown(it);
    SEQAN_ASSERT_EQ(getSubTreeSize(it), 0u);
}

template <typename TIter>
void waveletTreeStructureIteratorGetPosition(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it = begin(waveletTreeStructure, typename Spec<TIter>::Type());

    SEQAN_ASSERT_EQ(getPosition(it), 0u);

    goLeftChild(it);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    goRight(it);
    SEQAN_ASSERT_EQ(getPosition(it), 2u);
}

template <typename TIter>
void waveletTreeStructureIteratorGoChild(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(goLeftChild(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    SEQAN_ASSERT_EQ(goLeftChild(it), false);
    SEQAN_ASSERT_EQ(goRightChild(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    goRight(it);

    SEQAN_ASSERT_EQ(goLeftChild(it), false);
    SEQAN_ASSERT_EQ(goRightChild(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 3u);
}

template <typename TIter>
void waveletTreeStructureIteratorGoDown(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(goDown(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    SEQAN_ASSERT_EQ(goDown(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    goUp(it);
    goRightChild(it);
    SEQAN_ASSERT_EQ(goDown(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 3u);
}

template <typename TIter>
void waveletTreeStructureIteratorGoRight(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(goRight(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 0u);

    goDown(it);
    SEQAN_ASSERT_EQ(goRight(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 2u);

    SEQAN_ASSERT_EQ(goRight(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 2u);

    goDown(it);
    SEQAN_ASSERT_EQ(goRight(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 3u);
}

template <typename TIter>
void waveletTreeStructureIteratorGoToPosition(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(goToPosition(it, 1u), true);
    SEQAN_ASSERT_EQ(getPosition(it), 1u);

    SEQAN_ASSERT_EQ(goToPosition(it, 3u), true);
    SEQAN_ASSERT_EQ(getPosition(it), 3u);
}

template <typename TIter>
void waveletTreeStructureIteratorGoUp(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(goUp(it), false);
    SEQAN_ASSERT_EQ(getPosition(it), 0u);

    goDown(it);

    SEQAN_ASSERT_EQ(goUp(it), true);
    SEQAN_ASSERT_EQ(getPosition(it), 0u);
}

template <typename TIter>
void waveletTreeStructureIteratorIsLeaf(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(isLeaf(it), false);

    goDown(it);

    SEQAN_ASSERT_EQ(isLeaf(it), true);

    goRight(it);

    SEQAN_ASSERT_EQ(isLeaf(it), false);

    goDown(it);

    SEQAN_ASSERT_EQ(isLeaf(it), true);
}

template <typename TIter>
void waveletTreeStructureIteratorIsRoot(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(isRoot(it), true);
    goLeftChild(it);
    SEQAN_ASSERT_EQ(isRoot(it), false);
    goUp(it);
    goRightChild(it);
    SEQAN_ASSERT_EQ(isRoot(it), false);
    goRightChild(it);
    SEQAN_ASSERT_EQ(isRoot(it), false);
    goUp(it);
    goUp(it);
    SEQAN_ASSERT_EQ(isRoot(it), true);
}

template <typename TIter>
void _waveletTreeStructureIteratorSetAndGoRight(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    typedef String<TChar> TString;

    {
        TString text = "ACGTNACGTNACGTN";

        PrefixSumTable<TChar> pst(text);
        String<Pair<unsigned> > borderString;
        appendValue(borderString, Pair<unsigned>(0, 4));

        TRightArrayBinaryTree waveletTreeStructure;
        appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('G', 2u));
        appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('C', 0u));

        TIter it(waveletTreeStructure, 0);
        goDown(it);
        appendValue(borderString, Pair<unsigned>(0, 1));

        SEQAN_ASSERT_EQ(_setAndGoRight(it, borderString, pst), true);
        SEQAN_ASSERT_EQ(getPosition(it), 2u);
        SEQAN_ASSERT_EQ(borderString[1], Pair<unsigned>(2, 4));
    }
    {
        TString text = "TNACGTNACGTN";

        PrefixSumTable<TChar> pst(text);
        String<Pair<unsigned> > borderString;
        appendValue(borderString, Pair<unsigned>(0, 4));

        TRightArrayBinaryTree waveletTreeStructure;
        appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('T', 2u));
        appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('G', 2u));
        appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('C', 0u));

        TIter it(waveletTreeStructure, 0);
        goDown(it);
        appendValue(borderString, Pair<unsigned>(0, 2));

        SEQAN_ASSERT_EQ(_setAndGoRight(it, borderString, pst), true);
        SEQAN_ASSERT_EQ(getPosition(it), 3u);
        SEQAN_ASSERT_EQ(borderString[1], Pair<unsigned>(3, 4));
    }
}

template <typename TIter>
void waveletTreeStructureIteratorSetCharacter(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    String<TChar> text = "ACGTNACGTNACGTN";
    TRightArrayBinaryTree waveletTreeStructure(text);

    TIter it(waveletTreeStructure, 0);

    setCharacter(it, 'A');
    SEQAN_ASSERT_EQ(getCharacter(it), 'A');
}

template <typename TIter>
void _waveletTreeStructureSetChildVertices(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    typedef String<TChar> TString;

    TString text = "ACGTNACGTNACGTN";

    PrefixSumTable<TChar> pst(text);
    String<Pair<unsigned> > borderString;
    appendValue(borderString, Pair<unsigned>(0, 4));

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    _setChildVertices(it, borderString, pst);

    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);
    SEQAN_ASSERT_EQ(getCharacter(it), 'G');

    waveletTreeStructure.treeVertices[0].i2 = 4;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('C', 0u));
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('A', 0u));

    goRightChild(it);
    appendValue(borderString, Pair<unsigned>(2, 4));

    _setChildVertices(it, borderString, pst);

    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 3u);
    SEQAN_ASSERT_EQ(getCharacter(it), 'T');
}

template <typename TIter>
void _waveletTreeStructureSetLeftChildPos(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(_setLeftChildPos(it), true);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);

    waveletTreeStructure.treeVertices[0].i2 = 1;
    SEQAN_ASSERT_EQ(_setLeftChildPos(it), false);
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertices[0].i2, 1u);

    waveletTreeStructure.treeVertices[0].i2 = 2;
    SEQAN_ASSERT_EQ(_setLeftChildPos(it), true);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);

    waveletTreeStructure.treeVertices[0].i2 = 3;
    SEQAN_ASSERT_EQ(_setLeftChildPos(it), false);
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertices[0].i2, 3u);
}

template <typename TIter>
void _waveletTreeStructureSetPosition(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('A', 0u));
    resize(waveletTreeStructure.treeVertices, 10);

    TIter it(waveletTreeStructure, 0);

    goToPosition(it, 9u);

    SEQAN_ASSERT_EQ(getPosition(it), 9u);
}

template <typename TIter>
void _waveletTreeStructureSetRightChildPos(TIter & /*tag*/)
{
    typedef typename Container<TIter>::Type TRightArrayBinaryTree;
    typedef typename Fibre<TRightArrayBinaryTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    //typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
    //typedef String<TChar> TString;

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertices, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(_setRightChildPos(it, 0u), true);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 1u);

    waveletTreeStructure.treeVertices[0].i2 = 1;
    SEQAN_ASSERT_EQ(_setRightChildPos(it, 0u), true);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 1u);

    waveletTreeStructure.treeVertices[0].i2 = 2;
    SEQAN_ASSERT_EQ(_setRightChildPos(it, 11u), true);
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 11u);

    waveletTreeStructure.treeVertices[0].i2 = 3;
    SEQAN_ASSERT_EQ(_setRightChildPos(it, 11u), false);
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertices[0].i2, 3u);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_begin)
{
    using namespace seqan;

    typedef RightArrayBinaryTree<Dna5, void> TWaveletTreeStructure;
    typedef TopDown<> TIterSpec;
    typedef typename Iterator<TWaveletTreeStructure, TIterSpec>::Type TIter;

    TIter tag; // (waveletTreeStructure, 0);
    waveletTreeStructureIteratorBegin(tag);
}

// NOTE(esiragusa): IteratorContainer test should not use operator==()
//SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_container)
//{
//    using namespace seqan;
//
//    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
//    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
//    TIter tag(waveletTreeStructure, 0);
//    waveletTreeStructureIteratorContainer(tag);
//}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_end)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorEnd(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_get_character)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGetCharacter(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_get_child_pos)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGetChildPos(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_get_num_child_vertices)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGetNumChildVertices(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_get_position)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGetPosition(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_go_child)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGoChild(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_go_down)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGoDown(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_go_right)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGoRight(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_go_to_position)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGoToPosition(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_go_up)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGoUp(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_is_leaf)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorIsLeaf(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_is_root)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorIsRoot(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_and_go_right)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    _waveletTreeStructureIteratorSetAndGoRight(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_character)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorSetCharacter(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_child_vertices_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    _waveletTreeStructureSetChildVertices(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_left_child_pos_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    _waveletTreeStructureSetLeftChildPos(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_position_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    _waveletTreeStructureSetPosition(tag);
}


SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_right_child_pos_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    _waveletTreeStructureSetRightChildPos(tag);
}

#endif  // TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_

