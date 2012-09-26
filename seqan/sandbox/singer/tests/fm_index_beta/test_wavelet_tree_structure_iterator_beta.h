// ==========================================================================
//                               fm_index_beta
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_
#define TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/fm_sequence.h>

using namespace seqan;

template <typename TIter>
void waveletTreeStructureIteratorBegin(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

	String<TChar> text = "ACGTNACGTNACGTN";
	TRightArrayBinaryTree waveletTreeStructure(text);

	SEQAN_ASSERT_EQ(getCharacter(begin(waveletTreeStructure, typename Spec<TIter>::Type())), 'G');
}

template <typename TIter>
void waveletTreeStructureIteratorContainer(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;

	String<TChar> text = "ACGTNACGTNACGTN";
	TRightArrayBinaryTree waveletTreeStructure(text);

	TIter it(waveletTreeStructure, 0);	

	SEQAN_ASSERT(container(it) == waveletTreeStructure);
}

template <typename TIter>
void waveletTreeStructureIteratorEnd(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
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
void waveletTreeStructureIteratorGetNumChildVertieces(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('G', 4u));
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('C', 0u));
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('T', 1u));
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('N', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(getNumChildVertieces(it), 3u);

    goDown(it);
    SEQAN_ASSERT_EQ(getNumChildVertieces(it), 0u);
  
    goRight(it);
    SEQAN_ASSERT_EQ(getNumChildVertieces(it), 1u);
    
    goDown(it);
    SEQAN_ASSERT_EQ(getNumChildVertieces(it), 0u);
}

template <typename TIter>
void waveletTreeStructureIteratorGetPosition(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

	String<TChar> text = "ACGTNACGTNACGTN";
	TRightArrayBinaryTree waveletTreeStructure(text);

	TIter it(waveletTreeStructure, 0);

	SEQAN_ASSERT_EQ(goToPosition(it, 1u), true);
	SEQAN_ASSERT_EQ(getPosition(it), 1u);
	
	SEQAN_ASSERT_EQ(goToPosition(it, 3u), true);
	SEQAN_ASSERT_EQ(getPosition(it), 3u);

	SEQAN_ASSERT_EQ(goToPosition(it, 10u), false);
	SEQAN_ASSERT_EQ(getPosition(it), 3u);	
}

template <typename TIter>
void waveletTreeStructureIteratorGoUp(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

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
void waveletTreeStructureIteratorSetAndGoRight_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

    {
        TString text = "ACGTNACGTNACGTN";

        PrefixSumTable<TChar> pst(text);
        String<Pair<unsigned> > borderString;
        appendValue(borderString, Pair<unsigned>(0, 4));

        TRightArrayBinaryTree waveletTreeStructure;
        appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('G', 2u));
        appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('C', 0u));

        TIter it(waveletTreeStructure, 0);
        goDown(it);
        appendValue(borderString, Pair<unsigned>(0, 1));

        SEQAN_ASSERT_EQ(setAndGoRight_(it, borderString, pst), true);
        SEQAN_ASSERT_EQ(getPosition(it), 2u);
        SEQAN_ASSERT_EQ(borderString[1], Pair<unsigned>(2, 4));
    }
    {
        TString text = "TNACGTNACGTN";

        PrefixSumTable<TChar> pst(text);
        String<Pair<unsigned> > borderString;
        appendValue(borderString, Pair<unsigned>(0, 4));

        TRightArrayBinaryTree waveletTreeStructure;
        appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('T', 2u));
        appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('G', 2u));
        appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('C', 0u));

        TIter it(waveletTreeStructure, 0);
        goDown(it);
        appendValue(borderString, Pair<unsigned>(0, 2));

        SEQAN_ASSERT_EQ(setAndGoRight_(it, borderString, pst), true);
        SEQAN_ASSERT_EQ(getPosition(it), 3u);
        SEQAN_ASSERT_EQ(borderString[1], Pair<unsigned>(3, 4));
    }
}

template <typename TIter>
void waveletTreeStructureIteratorSetCharacter(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

	String<TChar> text = "ACGTNACGTNACGTN";
	TRightArrayBinaryTree waveletTreeStructure(text);

	TIter it(waveletTreeStructure, 0);

    setCharacter(it, 'A');
	SEQAN_ASSERT_EQ(getCharacter(it), 'A');
}

template <typename TIter>
void waveletTreeStructureSetChildVertices_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;

	TString text = "ACGTNACGTNACGTN";

    PrefixSumTable<TChar> pst(text);
    String<Pair<unsigned> > borderString;
    appendValue(borderString, Pair<unsigned>(0, 4));

    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    setChildVertieces_(it, borderString, pst);
    
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);
    SEQAN_ASSERT_EQ(getCharacter(it), 'G');
   
    waveletTreeStructure.treeVertieces[0].i2 = 4;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('C', 0u));
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    goRightChild(it);
    appendValue(borderString, Pair<unsigned>(2, 4));
    
    setChildVertieces_(it, borderString, pst);

    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 3u);
    SEQAN_ASSERT_EQ(getCharacter(it), 'T');
}

template <typename TIter>
void waveletTreeStructureSetLeftChildPos_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;
 
    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(setLeftChildPos_(it), true);  
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);

    waveletTreeStructure.treeVertieces[0].i2 = 1;
    SEQAN_ASSERT_EQ(setLeftChildPos_(it), false);  
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertieces[0].i2, 1u);

    waveletTreeStructure.treeVertieces[0].i2 = 2;
    SEQAN_ASSERT_EQ(setLeftChildPos_(it), true);  
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);

    waveletTreeStructure.treeVertieces[0].i2 = 3;
    SEQAN_ASSERT_EQ(setLeftChildPos_(it), false);  
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertieces[0].i2, 3u);
}

template <typename TIter>
void waveletTreeStructureSetPosition_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;
 
    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));
    resize(waveletTreeStructure.treeVertieces, 10);

    TIter it(waveletTreeStructure, 0);

    setPosition_(it, 9u);

    SEQAN_ASSERT_EQ(getPosition(it), 9u);  
}

template <typename TIter>
void waveletTreeStructureSetRightChildPos_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;
 
    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(setRightChildPos_(it,0), true);  
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 1u);

    waveletTreeStructure.treeVertieces[0].i2 = 1;
    SEQAN_ASSERT_EQ(setRightChildPos_(it,0), true);  
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 0u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 1u);

    waveletTreeStructure.treeVertieces[0].i2 = 2;
    SEQAN_ASSERT_EQ(setRightChildPos_(it,11), true);  
    SEQAN_ASSERT_EQ(getLeftChildPos(it), 1u);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 11u);

    waveletTreeStructure.treeVertieces[0].i2 = 3;
    SEQAN_ASSERT_EQ(setRightChildPos_(it,11), false);  
    SEQAN_ASSERT_EQ(waveletTreeStructure.treeVertieces[0].i2, 3u);
}

template <typename TIter>
void waveletTreeStructureSetRightChildPosOnly_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;
 
    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);
    setRightChildPosOnly_(it);
    SEQAN_ASSERT_EQ(getRightChildPos(it), 1u); 
}

template <typename TIter>
void waveletTreeStructureSetVertexToLeaf_(TIter & /*tag*/)
{
	typedef typename Container<TIter>::Type TRightArrayBinaryTree;
	typedef typename Fibre<TRightArrayBinaryTree, FibreTreeVertieces>::Type TWaveletTreeVertieces;
	typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
	typedef typename Value<TWaveletTreeVertex, 1>::Type TChar;
	typedef String<TChar> TString;
 
    TRightArrayBinaryTree waveletTreeStructure;
    appendValue(waveletTreeStructure.treeVertieces, TWaveletTreeVertex('A', 0u));

    TIter it(waveletTreeStructure, 0);

    SEQAN_ASSERT_EQ(setVertexToLeaf_(it), true);  
   
    waveletTreeStructure.treeVertieces[0].i2 = 1;
    SEQAN_ASSERT_EQ(setVertexToLeaf_(it), false);  
    
    waveletTreeStructure.treeVertieces[0].i2 = 2;
    SEQAN_ASSERT_EQ(setVertexToLeaf_(it), false);  
    
    waveletTreeStructure.treeVertieces[0].i2 = 3;
    SEQAN_ASSERT_EQ(setVertexToLeaf_(it), false);  
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

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_container)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorContainer(tag);
}

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

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_get_num_child_vertieces)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorGetNumChildVertieces(tag);
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
    waveletTreeStructureIteratorSetAndGoRight_(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_character)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureIteratorSetCharacter(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_child_vertieces_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetChildVertices_(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_left_child_pos_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetLeftChildPos_(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_position_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetPosition_(tag);
}


SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_right_child_pos_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetRightChildPos_(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_right_child_pos_only_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetRightChildPosOnly_(tag);
}

SEQAN_DEFINE_TEST(wavelet_tree_structure_iterator_set_vertex_to_leaf_)
{
    using namespace seqan;

    RightArrayBinaryTree<Dna5, void> waveletTreeStructure;
    typedef typename Iterator<RightArrayBinaryTree<Dna5, void>, TopDown<ParentLinks<> > >::Type TIter;
    TIter tag(waveletTreeStructure, 0);
    waveletTreeStructureSetVertexToLeaf_(tag);
}

#endif  // TESTS_WAVELT_TREE_STRUCTURE_ITERATOR_BETA_H_

