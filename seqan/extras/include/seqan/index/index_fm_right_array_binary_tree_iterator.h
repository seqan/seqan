// ==========================================================================
//                                  wavelet tree
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

namespace seqan {

template <typename TSpec>
struct RightArrayBinaryTreeIterator;

// ==========================================================================
// Metafunctions
// ==========================================================================
template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > >
{
    typedef RightArrayBinaryTree<TChar, TSpec> Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > >
{
    typedef RightArrayBinaryTree<TChar, TSpec> const Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > >:
    Container<Iter<RightArrayBinaryTree<TChar, TSpec>, TopDown<> > >
{};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Container<Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > >:
    Container<Iter<RightArrayBinaryTree<TChar, TSpec> const, TopDown<> > >
{};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<> > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec> const, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<> > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<TIterSpec> > >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<> > > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec> const, TopDown<ParentLinks<TIterSpec> > >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<> > > > Type;
};

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

// ==========================================================================
// Classes
// ==========================================================================

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

///.Metafunction.Iterator.param.T.type:Class.RightArrayBinaryTree

template <typename TChar, typename TSpec, typename TIterSpec>
class Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    TPos position;
    RightArrayBinaryTree<TChar, TSpec> const * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}

    template <typename TPos>
    Iter(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, const TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TChar, typename TSpec, typename TIterSpec>
class Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    TPos position;
    RightArrayBinaryTree<TChar, TSpec> * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}


    template <typename TPos>
    Iter(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TChar, typename TSpec>
class Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<> > > >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    String<TPos> position;
    RightArrayBinaryTree<TChar, TSpec> * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}

    template <typename TPos>
    Iter(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TPos pos) :
        position(),
        waveletTreeStructure(&treeStructure)
    {
        appendValue(position, pos);
    }

};

template <typename TChar, typename TSpec>
class Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<> > > >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    String<TPos> position;
    RightArrayBinaryTree<TChar, TSpec> const * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}

    template <typename TPos>
    Iter(const RightArrayBinaryTree<TChar, TSpec> & treeStructure, const TPos pos) :
        position(),
        waveletTreeStructure(&treeStructure)
    {
        appendValue(position, pos);
    }

};

// ==========================================================================
// Functions
// ==========================================================================

///.Function.begin.param.object.type:Class.RightArrayBinaryTree
template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> const & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type(waveletTreeStructure, 0);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type(waveletTreeStructure, 0);
}

///.Function.container.param.object.type:Class.RightArrayBinaryTree
template <typename TChar, typename TSpec, typename TIterSpec>
inline RightArrayBinaryTree<TChar, TSpec> &
container(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    return *it.waveletTreeStructure;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline RightArrayBinaryTree<TChar, TSpec> const &
container(Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    return *it.waveletTreeStructure;
}

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

/**
.Function.getCharacter
..signature:getCharacter(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator
..include:seqan/index.h
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline TChar getCharacter(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> const & iter)
{
    return iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i1;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline TChar getCharacter(Iter<const RightArrayBinaryTree<TChar, TSpec>, TIterSpec> const & iter)
{
    return iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i1;
}

/**
.Function.getLeftChildPos
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the left child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<const RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

/**
.Function.getNumChildVertieces
..summary:Returns the number of vertices in the subtree starting at the position an iterator points to.
..signature:getNumChildVertieces(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned getNumChildVertieces(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> it)
{
    unsigned originalPos = getPosition(it);
    while (goRightChild(it) || goLeftChild(it))
        continue;

    unsigned newPos = getPosition(it);
    goToPosition(it, originalPos);

    return newPos - originalPos;
}

/**
.Function.getPosition
..summary:Returns the position of the iterator in the host.
..signature:getPosition(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> const & iter)
{
    return iter.position;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<const RightArrayBinaryTree<TChar, TSpec>, TIterSpec> const & iter)
{
    return iter.position;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > const & iter)
{
    return iter.position[length(iter.position) - 1];
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getPosition(Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > const & iter)
{
    return iter.position[length(iter.position) - 1];
}

/**
.Function.getRightChildPos
..summary:Returns the position in @Class.RightArrayBinaryTree@ of the right child vertex.
..signature:getLeftChildPos(it)
..param.it:The iterator.
...type:Spec.RightArrayBinaryTree Iterator.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getRightChildPos(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 > 2)
    {
        return iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 == 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline unsigned int getRightChildPos(Iter<const RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 > 2)
    {
        return iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 == 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

/*
.Function.goDown
..param.iterator
...type:Spec.RightArrayBinaryTree Iterator
*/
/**
.Function.goDown.param.iterator.type:Spec.RightArrayBinaryTree Iterator
*/
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goDown(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (goLeftChild(iter))
        return true;

    if (goRightChild(iter))
        return true;

    return false;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goDown(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & iter)
{
    if (goLeftChild(iter))
        return true;

    if (goRightChild(iter))
        return true;

    return false;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goDownConstruction(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it)
{
    if (goLeftChild(it))
    {
        resize(container(it).treeVertieces, length(container(it).treeVertieces) + 1);
        return true;
    }

    if (goRightChild(it))
    {
        resize(container(it).treeVertieces, length(container(it).treeVertieces) + 1);
        return true;
    }

    return false;
}

/**
.Function.goLeftChild
..signature:bool goLeftChild(iterator)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.TopDown Iterator
...type:Spwc.RightArrayBinaryTree Iterator
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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    setPosition_(iter, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    setPosition_(iter, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    appendValue(iter.position, leftChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goLeftChild(Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
{
    unsigned leftChildPos = getLeftChildPos(iter);
    if (leftChildPos == 0)
        return false;

    appendValue(iter.position, leftChildPos);
    return true;
}

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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRight(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{

    unsigned pos = getPosition(iter);
    if (goUp(iter))
    {
        if (goRightChild(iter))
        {
            if (pos != getPosition(iter))
                return true;
        }
        else
        {
            goToPosition(iter, pos);
        }
    }

    return false;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRight(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & iter)
{
    unsigned pos = getPosition(iter);
    if (goUp(iter))
    {
        if (goRightChild(iter))
        {
            if (pos != getPosition(iter))
                return true;
        }
        else
        {
            goToPosition(iter, pos);
        }

    }

    return false;
}

/**
.Function.goRightChild
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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    setPosition_(iter, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<const RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    setPosition_(iter, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    appendValue(iter.position, rightChildPos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goRightChild(Iter<const RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
{
    unsigned rightChildPos = getRightChildPos(iter);
    if (rightChildPos == 0)
        return false;

    appendValue(iter.position, rightChildPos);
    return true;
}

/**
.Function.goToPosition
..signature:bool goToPosition(iterator, pos)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..param.pos:A position.
..remarks:$goToPosition(iterator)$ goes to position pos if it exist.
..returns:$true$ if the edge or path to go down exists, otherwise $false$.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
RightArrayBinaryTree<Dna5> waveletTreeStructure(genome);

Iterator<RightArrayBinaryTree<Dna5>, TopDown<> >::Type it;
it = begin(waveletTreeStructure); // go to root node

goToPosition(it, 2); // go to right child of root node
*/

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter, TPos pos)
{
    if (pos >= length(container(iter).treeVertieces))
        return false;

    setPosition(iter, pos);
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter, TPos pos)
{
    if (pos >= length(container(iter).treeVertieces))
        return false;

    appendValue(iter.position, pos);
    return true;
}

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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goUp(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(it.position, treeLevel - 1);
    setPosition_(it, it.position[treeLevel - 2]);

    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool goUp(Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(it.position, treeLevel - 1);
    setPosition_(it, it.position[treeLevel - 2]);

    return true;
}

// This function implements the functionality of go up and
// resizes the borderString of the structure construction.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString>
inline bool goUpStructureConstruction_(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TBorderString & borderString)
{
    unsigned treeLevel = length(it.position);

    if (isRoot(it))
        return false;

    resize(borderString, treeLevel - 1);
    resize(it.position, treeLevel - 1);
    setPosition_(it, it.position[treeLevel - 2]);

    return true;
}

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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isLeaf(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & iter)
{
    return !((*iter.waveletTreeStructure).treeVertieces[getPosition(iter)].i2);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isLeaf(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    return !((*iter.waveletTreeStructure).treeVertieces[getPosition(iter)].i2);
}

// This function creates the right sibling of the current node
// and goes to that one.
// Note: It acn only be called, if the right sibling really exists!
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString, typename TPrefixSumTable>
inline bool setAndGoRight_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it, TBorderString & borderString, TPrefixSumTable & pst)
{
    if (isRoot(it) || (borderString[length(borderString) - 1].i2 == borderString[length(borderString) - 2].i2))
        return false;

    goUp(it);

    if (borderString[length(borderString) - 2].i2 == ordValue(getCharacter(it)))
    {
        goLeftChild(it);
        return false;
    }

    resize(container(it).treeVertieces, length(container(it).treeVertieces) + 1);
    TChar pivot = getCharacter(it);
    // SEQAN_ASSERT_MSG(setRightChildPos_(it, length(container(it).treeVertieces) - 1), "You just deleted inserted vertieves!");
    setRightChildPos_(it, length(container(it).treeVertieces) - 1);
    goRightChild(it);
    borderString[length(borderString) - 1].i1 = getCharacterPosition(pst, pivot);
    borderString[length(borderString) - 1].i2 = borderString[length(borderString) - 2].i2;

    return true;
}

/**
.Function.setCharacter
..signature:bool setCharacter(iterator, character)
..param.iterator:An iterator of a wavelet tree structure.
...type:Spec.RightArrayBinaryTree Iterator
..param.character:The character to be assigned to a node.
..remarks:$setCharacter(iterator, character)$ sets the character of the node the iterator points to to character.
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
template <typename TChar, typename TSpec, typename TIterSpec, typename TChar2>
inline void setCharacter(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter,
                         TChar2 character)
{
    iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i1 = character;
}

// This function sets the left child of the current node, or the right if there is no left child.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString, typename TCharPST, typename TSpecPST>
void setChildVertieces_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                        TBorderString & borderString,
                        PrefixSumTable<TCharPST, TSpecPST> & pst)
{
    typedef typename Value<TBorderString>::Type TBorderStringValue;
    unsigned leftBorder = borderString[length(borderString) - 1].i1;
    unsigned rightBorder = borderString[length(borderString) - 1].i2;
    unsigned pivotPosition = getPivotPosition(pst, leftBorder, rightBorder);

    setCharacter(it, getCharacter(pst, pivotPosition));

    if (leftBorder == pivotPosition - 1)
    {
        setRightChildPosOnly_(it);
        appendValue(borderString, TBorderStringValue(pivotPosition, borderString[length(borderString) - 1].i2));
        return;
    }

//     SEQAN_ASSERT_MSG(setLeftChildPos_(it), "The right child has just been deleted!");
    setLeftChildPos_(it);

    appendValue(borderString, TBorderStringValue(borderString[length(borderString) - 1].i1, pivotPosition - 1));
}

// This functions sets the pointer to the left child.
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool setLeftChildPos_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    switch (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2)
    {
    case (0):
        iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 = 2;
        return true;

    case (2):
        return true;

    default:
        return false;
    }

//     if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 == 0)
//     {
//         iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 = 2;
//         return true;
//     }
//     if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 == 2)
//     {
//         return true;
//     }
//     //std::cerr << "ERROR: The right child has just been deleted!" << std::endl;
//     return false;
}

// This function sets the position of iter to pos.
template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it, TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, length(getFibre(container(it), FibreTreeVertieces())), "The position does not exist");
    it.position = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition_(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & it, TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, length(getFibre(container(it), FibreTreeVertieces())), "The position does not exist");
    it.position = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition_(Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, length(getFibre(container(it), FibreTreeVertieces())), "The position does not exist");
    it.position[length(it.position) - 1] = pos;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TPos>
inline void setPosition_(Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TPos pos)
{
    SEQAN_ASSERT_LT_MSG(pos, length(getFibre(container(it), FibreTreeVertieces())), "The position does not exist");
    it.position[length(it.position) - 1] = pos;
}

// This functions sets the pointer to the left child.
template <typename TChar, typename TSpec, typename TPos, typename TIterSpec>
inline bool setRightChildPos_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter, TPos rightChildPosition)
{
    switch (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2)
    {
    case (0):
        SEQAN_ASSERT_MSG(rightChildPosition == 0u, "Wrong right child position!");
        iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 = 1;
        return true;

    case (2):
        iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 = rightChildPosition + 2;
        return true;

    case (1):
        SEQAN_ASSERT_MSG(rightChildPosition == 0u, "Wrong right child position!");
        return true;

    default:
        return false;
    }
}

// This function sets the pointer to the right child such that it is cleat that there is no left child.
template <typename TChar, typename TSpec, typename TIterSpec>
inline void setRightChildPosOnly_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 = 1;
}

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
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isRoot(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it)
{
    return !getPosition(it);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline bool isRoot(Iter<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec> & it)
{
    return !getPosition(it);
}

// This function sets the current node to be a node.
template <typename TChar, typename TSpec, typename TIterSpec>
inline bool setVertexToLeaf_(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter)
{
    if (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2 != 0)
    {
        std::cerr << "You just deleted ";
        switch (iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i2)
        {
        case (1):
            std::cerr << "the right sub tree!" << std::endl;
            return false;

        case (2):
            std::cerr << "the left sub tree!" << std::endl;
            return false;

        default:
            std::cerr << "both sub trees!" << std::endl;
            return false;
        }
    }
    return true;
}

template <typename TChar, typename TSpec, typename TIterSpec, typename TString>
inline void writeGraphImpl(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & iter, TString name)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> iter2 = iter;
    std::ofstream stream(toCString(name), std::ios::app);
    unsigned pos = getLeftChildPos(iter);
    if (pos)
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i1) << " -> " << (unsigned)ordValue(iter.waveletTreeStructure->treeVertieces[pos].i1) << ";" << std::endl;
        goLeftChild(iter);
        writeGraphImpl(iter, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertieces[getPosition(iter)].i1) << " -> " << "leave1" << (unsigned)ordValue(getPosition(iter)) << ";" << std::endl;
    }

    pos = getRightChildPos(iter2);
    if (pos)
    {
        stream << (unsigned) ordValue(iter2.waveletTreeStructure->treeVertieces[getPosition(iter2)].i1) << " -> " << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertieces[pos].i1) << ";" << std::endl;

        goRightChild(iter2);
        writeGraphImpl(iter2, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertieces[getPosition(iter2)].i1) << " -> " << "leave2" << (unsigned)ordValue(getPosition(iter2)) << ";" << std::endl;
    }
    stream.close();
}

template <typename TChar, typename TSpec>
inline void writeGraph(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
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

}
#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H
