// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
//       from this software withoFIut specific prior written permission.
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

#ifndef INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H
#define INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H

//SEQAN_NO_DDDOC:do not generate documentation for this file

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
struct RightArrayBinaryTreeIterator;


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec> const, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TIterSpec>
struct Spec<Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TIterSpec> > >
{
    typedef TIterSpec Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Spec<Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TIterSpec> > >
{
    typedef TIterSpec Type;
};

// ============================================================================
// Classes
// ============================================================================
/*!
 * @class RightArrayBinaryTreeIterator RightArrayBinaryTree Iterator
 * 
 * @extends Iter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief An iterator for @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iter<RightArrayBinaryTree, TSpec >
 * 
 * @tparam TSpec Specialisation Tag. Types: TopDownIterator
 * @tparam RightArrayBinaryTree The @link RightArrayBinaryTree @endlink. Types:
 *                              @link RightArrayBinaryTree @endlink
 */
/**
.Spec.RightArrayBinaryTree Iterator:
..summary:An iterator for @Class.RightArrayBinaryTree@.
..cat:Iter
..general:Class.Iter
..signature:Iter<RightArrayBinaryTree, TSpec >
..param.RightArrayBinaryTree:The @Class.RightArrayBinaryTree@.
...type:Class.RightArrayBinaryTree
..param.TSpec:Specialisation Tag.
...type:Spec.TopDown Iterator
..include:seqan/index.h
*/

template <typename TTree, typename TIterSpec>
class Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > >
{
    typedef typename Fibre<TTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    TPos position;
    TTree * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}

    Iter(TTree & treeStructure, TPos pos = 0) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TTree, typename TSpec>
class Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TSpec> > > > :
    public Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > >
{
    typedef Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > > TBase;
    typedef typename Fibre<TTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    String<TPos, Block<> > history;

    Iter() :
        TBase(),
        history()
    {}
    
    Iter(TTree & treeStructure, TPos pos = 0) :
        TBase(treeStructure,  pos),
        history()
    {
        appendValue(history, pos);
    }

};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#begin
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The begin (root) of a @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iterator begin(rightArrayBinaryTree, iterSpec)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree.
 *
 * @param iterSpec A specialisation tag. Types: TopDown<>, TopDown<ParentLinks<> >
 * 
 * @return TReturn An iterator to the first item in <tt>object</tt>.
 *                 Metafunctions: Metafunction.Iterator
 */
///.Function.begin.param.object.type:Class.RightArrayBinaryTree
template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> const & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type(waveletTreeStructure);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type(waveletTreeStructure);
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#container
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Container of an iterator.
 * 
 * @signature Container container(iterator)
 * 
 * @param iterator An iterator.
 * 
 * @return TReturn The container that <tt>iterator</tt> traverses.
 */
///.Function.container.param.iterator.type:Class.RightArrayBinaryTree
template <typename TTree, typename TIterSpec>
inline TTree &
container(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    return *it.waveletTreeStructure;
}

template <typename TTree, typename TIterSpec>
inline TTree &
container(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return *it.waveletTreeStructure;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#end
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The end (rigthmost laef) of a @link RightArrayBinaryTree @endlink.
 * 
 * @signature Iterator end(rightArrayBinaryTree, iterSpec)
 * 
 * @param rightArrayBinaryTree The right-array-binary tree.
 *
 * @param iterSpec A specialisation tag. Types: TopDown<>, TopDown<ParentLinks<> >
 * 
 * @return TReturn An iterator to the first item in <tt>object</tt>.
 *                 Metafunctions: Metafunction.Iterator
 */
///.Function.end.param.object.type:Class.RightArrayBinaryTree
template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type
end(RightArrayBinaryTree<TChar, TSpec> const & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type(waveletTreeStructure, length(waveletTreeStructure));
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type
end(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type(waveletTreeStructure, length(waveletTreeStructure));
}

// ----------------------------------------------------------------------------
// Function getCharacter()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getCharacter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief This function returns the pivot character of the node the iterator
 *        currently points to.
 * 
 * @signature getCharacter(iterator)
 * 
 * @param iterator The iterator.
 */
/**
.Function.getCharacter
..summary:This function returns the pivot character of the node the iterator currently points to.
..signature:getCharacter(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator
..include:seqan/index.h
*/
template <typename TTree, typename TIterSpec>
inline typename Value<typename Value<TTree>::Type, 1>::Type
getCharacter(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & iter)
{
    return iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1;
}

// ----------------------------------------------------------------------------
// Function getLeftChildPos()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getLeftChildPos
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the
 *        left child node.
 * 
 * @signature getLeftChildPos(iterator)
 * 
 * @param iterator The iterator.
 */
/**
.Function.getLeftChildPos
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the left child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TTree, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & iter)
{
    if (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function getSubTreeSize()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getSubTreeSize
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the number of vertices in the subtree starting at the position
 *        an iterator points to.
 * 
 * @signature getSubTreeSize(iterator)
 * 
 * @param iterator The iterator.
 */
/**
.Function.getSubTreeSize
..summary:Returns the number of vertices in the subtree starting at the position an iterator points to.
..signature:getSubTreeSize(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TTree, typename TIterSpec>
inline unsigned getSubTreeSize(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > > _it(container(it));
    unsigned originalPos = getPosition(it);
    goToPosition(_it, originalPos);
    while (goRightChild(_it) || goLeftChild(_it))
        continue;

    return getPosition(_it) - originalPos;
}

// ----------------------------------------------------------------------------
// Function getPosition()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getPosition
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position of the iterator in the host.
 * 
 * @signature getPosition(iterator)
 * 
 * @param iterator The iterator.
 */
/**
.Function.getPosition
..summary:Returns the position of the iterator in the host.
..signature:getPosition(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TTree, typename TIterSpec>
inline unsigned int getPosition(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return it.position;
}

// ----------------------------------------------------------------------------
// Function getRightChildPos()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getRightChildPos
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the
 *        right child node.
 * 
 * @signature getLeftChildPos(iterator)
 * 
 * @param iterator The iterator.
 */
/**
.Function.getRightChildPos
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the right child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TTree, typename TIterSpec>
inline unsigned int getRightChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    if (it.waveletTreeStructure->treeVertices[getPosition(it)].i2 > 2)
    {
        return it.waveletTreeStructure->treeVertices[getPosition(it)].i2 - 2;
    }
    if (it.waveletTreeStructure->treeVertices[getPosition(it)].i2 == 1)
    {
        return getPosition(it) + 1;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _historyPush()
// ----------------------------------------------------------------------------

template <typename TTree, typename TIterSpec, typename TPos>
inline void _historyPush(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & /*tag*/ , TPos /*tag*/)
{}

template <typename TTree, typename TIterSpec, typename TPos>
inline void _historyPush(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TPos pos)
{
    appendValue(it.history, pos);
}

// ----------------------------------------------------------------------------
// Function goDown()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goDown
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates down the leftmost edge in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature bool goDown(iterator)
 *
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if an edge to go down exists,
 *                 otherwise <tt>false</tt>. Types: bool
 */
///.Function.goDown.param.iterator.type:Spec.RightArrayBinaryTree Iterator
template <typename TTree, typename TIterSpec>
inline bool goDown(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    if (goLeftChild(iter)) return true;
    if (goRightChild(iter)) return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _goDownConstruction()
// ----------------------------------------------------------------------------

template <typename TTree, typename TIterSpec>
inline bool _goDownConstruction(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    if (goDown(it))
    {
        resize(container(it).treeVertices, length(container(it).treeVertices) + 1);
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function goLeftChild()
// ----------------------------------------------------------------------------

/*!
 * @fn RightArrayBinaryTreeIterator#goLeftChild
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the iterator to the left child of the current node if it exists
 *        and returns true, otherwise the iterator does not change position and
 *        the function returns false.
 * 
 * @signature bool goLeftChild(iterator)
 * 
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if the edge to go down exists,
 *                 otherwise <tt>false</tt>.
 */
/**
.Function.goLeftChild
..summary:Sets the iterator to the left child of the current node if it exists and returns true, otherwise the iterator does not change position and the function returns false.
..signature:bool goLeftChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..remarks:$goLeftChild(iterator)$ goes down the left edge if it exist.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goLeftChild(it); // go to left child of root node
*/
template <typename TTree, typename TIterSpec>
inline bool goLeftChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    unsigned leftChildPos = getLeftChildPos(it);
    if (leftChildPos == 0)
        return false;
    
    if (!goToPosition(it, leftChildPos))
        return false;

    _historyPush(it, leftChildPos);
    return true;
}

// ----------------------------------------------------------------------------
// Function goRight()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goRight
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates to the next sibling in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature goRight(iterator)
 * 
 * @param iterator The iterator
 * 
 * @return TReturn <tt>true</tt> if the iterator could be moved, otherwise
 *                 <tt>false</tt>. Types: nolink:bool
 */
/**
.Function.goRight
..param.iterator:
...type:Spec.RightArrayBinaryTree Iterator
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goDown(it); // go to left child of root node
goRight(it); // go to right child of root node
*/
template <typename TTree, typename TIterSpec>
inline bool goRight(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    unsigned pos = getPosition(it);
    if (goUp(it))
    {
        if (goRightChild(it))
        {
            if (pos != getPosition(it))
                return true;
        }
        else
            goLeftChild(it);
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function goRightChild()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goRightChild
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Sets the iterator to the right child of the current node if it exists
 *        and returns true, otherwise the iterator does not change position and
 *        the function returns false.
 * 
 * @signature bool goRightChild(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if the edge to go down exists,
 *                 otherwise <tt>false</tt>.
 */
/**
.Function.goRightChild
..summary:Sets the iterator to the right child of the current node if it exists and returns true, otherwise the iterator does not change position and the function returns false.
..signature:bool goRightChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..remarks:$goRightChild(iterator)$ goes down the right edge if it exist.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
*/
template <typename TTree, typename TIterSpec>
inline bool goRightChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    unsigned rightChildPos = getRightChildPos(it);
    if (rightChildPos == 0)
        return false;

    if (!goToPosition(it, rightChildPos))
        return false;

    _historyPush(it, rightChildPos);
    return true;
}

// template <typename TTree, typename TIterSpec>
// inline bool goRightChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
// {
//     unsigned rightChildPos = getRightChildPos(iter);
//     if (rightChildPos == 0)
//         return false;
// 
//     appendValue(iter.position, rightChildPos);
//     return true;
// }

// ----------------------------------------------------------------------------
// Function goToPosition()
// ----------------------------------------------------------------------------

// TODO(singer): Make this work!
/*
.Function.goToPosition
..summary:Move the iterator to a specified position.
..signature:bool goToPosition(iterator, pos)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..param.pos:A position.
..remarks:$goToPosition(iterator)$ goes to position pos regardless of pos being a valid position.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goToPosition(it, 2); // go to right child of root node
*/
template <typename TTree, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it, TPos pos)
{
    it.position = pos;
    return true;
}

// ----------------------------------------------------------------------------
// Function goUp()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goUp
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Iterates up one edge to the parent in a @link RightArrayBinaryTree @endlink.
 * 
 * @signature goUp(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if the iterator could be moved, otherwise
 *                 <tt>false</tt>. Types: nolink:bool
 */
/**
.Function.goUp.param.iterator.type:Spec.TopDownHistory Iterator
*/
/*
..returns:$true$ if the current node is not the root node.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
goUp(it); // go to root node
*/
template <typename TTree, typename TIterSpec>
inline bool goUp(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it)
{
    unsigned treeLevel = length(it.history);

    if (isRoot(it))
        return false;

//     for (unsigned i = 0; i < length(it.history); ++i)
//         std::cerr << (int)it.history[i] << " ";
//     std::cerr << std::endl;
    resize(it.history, treeLevel - 1);
//     for (unsigned i = 0; i < length(it.history); ++i)
//         std::cerr << (int)it.history[i] << " ";
//     std::cerr << std::endl;
    goToPosition(it, back(it.history));
//     std::cerr << "done" << std::endl;

    return true;
}

// ----------------------------------------------------------------------------
// Function _goUpStructureConstruction()
// ----------------------------------------------------------------------------

// This function implements the functionality of go up and
// resizes the borderString of the structure construction.
template <typename TTree, typename TIterSpec, typename TBorderString>
inline bool _goUpStructureConstruction(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TBorderString & borderString)
{
    if (goUp(it))
    {
        resize(borderString, length(it.history));
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function isLeaf()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#isLeaf
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Tests whether a given node is a leaf or not.
 * 
 * @signature isLeaf(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn True if the node is a leaf.
 */
/**
.Function.isLeaf
..class:Spec.RightArrayBinaryTree Iterator
..param.iterator.type:Spec.RightArrayBinaryTree Iterator
..example:
...text:Code example for the @Spec.RightArrayBinaryTree Iterator@:
...code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
goUp(it); // go to root node
*/
template <typename TTree, typename TIterSpec>
inline bool isLeaf(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    return container(iter).treeVertices[getPosition(iter)].i2 == 0;
}

// ----------------------------------------------------------------------------
// Function _setAndGoRight()
// ----------------------------------------------------------------------------

// This function creates the right sibling of the current node
// and goes to that one.
// Note: It acn only be called, if the right sibling really exists!
template <typename TTree, typename TIterSpec, typename TBorderString>
inline bool _setAndGoRight(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it, TBorderString & borderString)
{
    typedef typename Value<typename Value<TTree>::Type, 1>::Type TChar;

    if (isRoot(it) || (back(borderString).i2 == borderString[length(borderString) - 2].i2))
        return false;

    goUp(it);

    if (borderString[length(borderString) - 2].i2 == ordValue(getCharacter(it)))
    {
        goLeftChild(it);
        return false;
    }

    resize(container(it).treeVertices, length(container(it).treeVertices) + 1);
    TChar pivot = getCharacter(it);
    _setRightChildPos(it, length(container(it).treeVertices) - 1);
    goRightChild(it);
    back(borderString).i1 = ordValue(pivot);
    back(borderString).i2 = borderString[length(borderString) - 2].i2;

    return true;
}

// ----------------------------------------------------------------------------
// Function setCharacter()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#setCharacter
 * 
 * @headerfile seqan/index.h
 * 
 * @brief The function sets the character of the
 *        node the iterator points to to character.
 * 
 * @signature void setCharacter(iterator, character)
 * 
 * @param character The character to be assigned to a node.
 * @param iterator The iterator.
 */
/**
.Function.setCharacter
..signature:bool setCharacter(iterator, character)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..param.character:The character to be assigned to a node.
..summary:$setCharacter(iterator, character)$ sets the character of the node the iterator points to to character.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goRightChild(it); // go to right child of root node
setCharacter(it,'T'); // sets the character of the root's
                      // right child to 'T'
*/
template <typename TTree, typename TIterSpec, typename TChar2>
inline void setCharacter(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter,
                         TChar2 character)
{
    container(iter).treeVertices[getPosition(iter)].i1 = character;
}

// ----------------------------------------------------------------------------
// Function _getPivotPosition()
// ----------------------------------------------------------------------------

// This function returns the position of the character which ensures that the sum of occurrences of the characters from
// beginPos to the computed pos and the sum of occurrences from the computed pos to endPos are about the same.
template <typename TPrefixSums, typename TBeginPos, typename TEndPos>
unsigned _getPivotPosition(TPrefixSums const & sums, TBeginPos beginPos, TEndPos endPos)
{
    TBeginPos realBeginPos = beginPos + 1;
    TEndPos realEndPos = endPos + 1;
    unsigned lengthRange = realEndPos - realBeginPos + 1;
    unsigned pivotPos = realBeginPos + lengthRange / 2 - 1;

    unsigned tooSmallValues = sums[beginPos];
    long currentMin = sums[realEndPos] + 1;

    if (sums[pivotPos] - tooSmallValues >= sums[realEndPos] - sums[pivotPos])
    {
        while ((pivotPos >= realBeginPos) && std::abs((long)(sums[pivotPos] - tooSmallValues) - (long)((sums[realEndPos] - sums[pivotPos]))) <= currentMin)
        {
            currentMin = std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos])));
            --pivotPos;
        }
        ++pivotPos;
    }
    else
    {
        while (std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos]))) < currentMin && (pivotPos < realEndPos))
        {
            currentMin = std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos])));
            ++pivotPos;
        }
        --pivotPos;
    }

    return pivotPos;
}

// ----------------------------------------------------------------------------
// Function _setChildVertices()
// ----------------------------------------------------------------------------

// This function sets the left child of the current node, or the right if there is no left child.
template <typename TTree, typename TIterSpec, typename TBorderString, typename TPrefixSums>
void _setChildVertices(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it,
                       TBorderString & borderString,
                       TPrefixSums & sums)
{
    typedef typename Value<TBorderString>::Type TBorderStringValue;
    unsigned leftBorder = back(borderString).i1;
    unsigned rightBorder = back(borderString).i2;
    unsigned pivotPosition = _getPivotPosition(sums, leftBorder, rightBorder);

    setCharacter(it, pivotPosition);

    if (leftBorder == pivotPosition - 1)
    {
        // set the right child to be the only one
        container(it).treeVertices[getPosition(it)].i2 = 1;
        appendValue(borderString, TBorderStringValue(pivotPosition, back(borderString).i2));
        return;
    }

    _setLeftChildPos(it);

    appendValue(borderString, TBorderStringValue(back(borderString).i1, pivotPosition - 1));
}

// ----------------------------------------------------------------------------
// Function _setLeftChildPos()
// ----------------------------------------------------------------------------

// This functions sets the pointer to the left child.
template <typename TTree, typename TIterSpec>
inline bool _setLeftChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    switch (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2)
    {
    case 0:
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = 2;
        return true;

    case 2:
        return true;

    default:
        return false;
    }
}

// ----------------------------------------------------------------------------
// Function _setRightChildPos()
// ----------------------------------------------------------------------------

// This functions sets the pointer to the left child.
template <typename TTree, typename TPos, typename TIterSpec>
inline bool _setRightChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter, TPos rightChildPosition)
{
    switch (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2)
    {
    case 0:
        SEQAN_ASSERT_EQ_MSG(rightChildPosition, 0u, "Wrong right child position!");
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = 1;
        return true;

    case 2:
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = rightChildPosition + 2;
        return true;

    case 1:
        SEQAN_ASSERT_MSG(rightChildPosition == 0u, "Wrong right child position!");
        return true;

    default:
        return false;
    }
}

// ----------------------------------------------------------------------------
// Function isRoot()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#isRoot
 * 
 * @headerfile seqan/index.h
 * 
 * @brief Test whether a iterator points to the root node.
 * 
 * @signature bool isRoot(iterator)
 * 
 * @param iterator The iterator.
 * 
 * @return TReturn <tt>true</tt> if <tt>iterator</tt> points to the root of the
 *                 tree, otherwise <tt>false</tt>. Types: nolink:bool
 */
/**
.Function.isRoot
..class:Spec.RightArrayBinaryTree Iterator
..param.iterator.type:Spec.RightArrayBinaryTree Iterator
..example
...text:Code example for the @Spec.RightArrayBinaryTree Iterator@:
...code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

isRoot(it) // returns true
*/
template <typename TTree, typename TIterSpec>
inline bool isRoot(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return getPosition(it) == 0;
}

/*
template <typename TTree, typename TIterSpec, typename TString>
inline void _writeGraphImpl(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter, TString name)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > iter2 = iter;
    std::ofstream stream(toCString(name), std::ios::app);
    unsigned pos = getLeftChildPos(iter);
    if (pos)
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1) << " -> " << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[pos].i1) << ";" << std::endl;
        goLeftChild(iter);
        writeGraphImpl(iter, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1) << " -> " << "leave1" << (unsigned)ordValue(getPosition(iter)) << ";" << std::endl;
    }

    pos = getRightChildPos(iter2);
    if (pos)
    {
        stream << (unsigned) ordValue(iter2.waveletTreeStructure->treeVertices[getPosition(iter2)].i1) << " -> " << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertices[pos].i1) << ";" << std::endl;

        goRightChild(iter2);
        writeGraphImpl(iter2, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertices[getPosition(iter2)].i1) << " -> " << "leave2" << (unsigned)ordValue(getPosition(iter2)) << ";" << std::endl;
    }
    stream.close();
}

template <typename TTree>
inline void _writeGraph(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{

    typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<> > >::Type iter(treeStructure, 0);

    String<char> name = "testfile.dot";
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}
*/

}
#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H
