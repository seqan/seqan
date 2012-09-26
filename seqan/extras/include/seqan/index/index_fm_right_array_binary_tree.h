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

#ifndef INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
#define INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H

namespace seqan {

#include <seqan/index_fm.h>

// ==========================================================================
//Forwards
// ==========================================================================
template <typename TText, typename TSpec>
class WaveletTree;

template <typename TChar, typename TSpec = void>
class RightArrayBinaryTree;

//////////////////////////////////////////////////////////////////////////////
// RightArrayBinaryTree fibres

/**
.Tag.RightArrayBinaryTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RightArrayBinaryTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a RightArrayBinaryTree.
..cat:RightArrayBinaryTree

..tag.FibreTreeVertieces:The string encoding the wavelet tree structure.


..see:Metafunction.Fibre
..see:Function.getFibre
..see:Spec.IndexEsa
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.RightArrayBinaryTree Fibres

struct FibreTreeVertieces_;
struct FibreWaveletTreeStructure_;

typedef Tag<FibreTreeVertieces_> const FibreTreeVertieces;
typedef Tag<FibreWaveletTreeStructure_> const FibreWaveletTreeStructure;

// ==========================================================================
// Metafunctions
// ==========================================================================

// TODO (singer): Dna, Dna5 and AS require onle unsigned char
template <typename TChar, typename TSpec>
struct Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>
{
    typedef typename IfC < BitsPerValue<TChar>::VALUE<17,
                                                      unsigned short,
                                                      unsigned int>::Type TPos;
    typedef String<Pair<TChar, TPos> > Type;
};

template <typename TChar, typename TSpec>
struct Reference<RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Value<RightArrayBinaryTree<TChar, TSpec> >::Type & Type;
};

template <typename TChar, typename TSpec>
struct Reference<const RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Value<RightArrayBinaryTree<TChar, TSpec> >::Type const Type;
};

template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    typedef Pair<TChar, TPos> Type;
};

template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> const>
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type TWaveletTreeVertieces;
    typedef typename Value<TWaveletTreeVertieces>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

    typedef Pair<TChar, TPos> Type;
};

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

/**
.Class.RightArrayBinaryTree:
..cat:WaveletTree
..summary:A special format to encode the structure of a wavelet tree. The structure is very space efficient because only one position is stored which encodes where the left and right subtree of a given node exist.
..signature:RightArrayBinaryTree<TValue, TSpec>
..param.TSpec:The value type, that is the type of the stored characters.
..param.TSpec:The wavelet tree structure specialisation.
...default:void.
..include:seqan/index.h
*/
template <typename TChar, typename TSpec>
class RightArrayBinaryTree
{
public:
    typename Fibre<RightArrayBinaryTree, FibreTreeVertieces>::Type treeVertieces;
    TChar minCharValue;

    RightArrayBinaryTree() :
        treeVertieces(),
        minCharValue()
    {}

    template <typename TText>
    RightArrayBinaryTree(TText const & text) :
        treeVertieces(),
        minCharValue()
    {
        computeRightArrayBinaryTree(*this,
                                    text);
    }

    inline RightArrayBinaryTree & operator=(RightArrayBinaryTree const & other)
    {
        treeVertieces = other.treeVertieces;
        minCharValue = other.minCharValue;
        return *this;
    }

    inline bool operator==(const RightArrayBinaryTree & b) const
    {
        return treeVertieces == b.treeVertieces;
    }

};

// ==========================================================================
// Functions
// ==========================================================================

/**
.Function.clear
..param.object:
...type:Class.RightArrayBinaryTree
*/
template <typename TChar, typename TSpec>
inline void clear(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{
    clear(treeStructure.treeVertieces);
}

/**
.Function.computeRightArrayBinaryTree
..summary:Computes the wavelet tree structure of a text.
..signature:computeRightArrayBinaryTree(waveletTreeStructure, text)
..param.waveletTreeStructure:A wavelet tree structure.
...type:Class.RightArrayBinaryTree
..param.text:A text.
...type:Class.String
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RightArrayBinaryTree<Dna5> waveletTreeStructure;
computeRightArrayBinaryTree(genome);
*/
template <typename TChar, typename TSpec, typename TText>
inline void computeRightArrayBinaryTree(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TText const & text)
{
    PrefixSumTable<TChar, void> pst(text);

    typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<> > >::Type it(waveletTreeStructure, 0u);

    computeRightArrayBinaryTree(it, pst);
}

// This function computes the wavelet tree structure contained in the lfTable.
template <typename TSpec, typename TPrefixSumTable, typename TText>
inline void computeRightArrayBinaryTree(LfTable<WaveletTree<TText, TSpec>, TPrefixSumTable> & lfTable)
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type TRightArrayBinaryTree;

    TRightArrayBinaryTree & rightArrayBinaryTree = lfTable.occTable.waveletTreeStructure;

    typename Iterator<TRightArrayBinaryTree, TopDown<ParentLinks<void> > >::Type it(rightArrayBinaryTree, 0u);

    computeRightArrayBinaryTree(it, lfTable.prefixSumTable);
}

// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TIterSpec, typename TPrefixSumTable>
inline void computeRightArrayBinaryTree(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                        TPrefixSumTable & pst)
{
    typedef RightArrayBinaryTree<TChar, TSpec> TRightArrayBinaryTree;

    TRightArrayBinaryTree & waveletTreeStructure = container(it);

    unsigned alpSize = getAlphabetSize(pst);
    String<Pair<unsigned> > borderString;
    appendValue(borderString, Pair<unsigned>(0, alpSize - 1));

    resize(waveletTreeStructure, 1);

    computeRightArrayBinaryTree(it, borderString, pst);

}

// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TBorderString, typename TPrefixSumTable, typename TIterSpec>
inline void computeRightArrayBinaryTree(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                        TBorderString & borderString,
                                        TPrefixSumTable & pst)
{
    do
    {
        if (borderString[length(borderString) - 1].i2 - borderString[length(borderString) - 1].i1 + 1 < 3 ||
            getPrefixSum(pst, borderString[length(borderString) - 1].i1) == getPrefixSum(pst, borderString[length(borderString) - 1].i2 + 1))
        {
            setCharacter(it, getCharacter(pst, (borderString[length(borderString) - 1].i1 + 1)));
            setVertexToLeaf_(it);
        }
        else
            setChildVertieces_(it, borderString, pst);

        if (!goDownConstruction(it) && !setAndGoRight_(it, borderString, pst))
            while (goUpStructureConstruction_(it, borderString) && !setAndGoRight_(it, borderString, pst))
                ;
    }
    while (!isRoot(it));
}

/**
.Function.empty
..param.object:
...type:Class.RightArrayBinaryTree
*/
template <typename TChar, typename TSpec>
inline bool empty(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{
    return empty(getFibre(treeStructure, FibreTreeVertieces()));
}

template <typename TChar, typename TSpec>
inline bool empty(RightArrayBinaryTree<TChar, TSpec> const & treeStructure)
{
    return empty(getFibre(treeStructure, FibreTreeVertieces()));
}

/**
.Function.getFibre
..param.container:
...type:Class.RightArrayBinaryTree
..param.fibreTag:
...type:Tag.WaveletTree Fibres
*/
template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type &
getFibre(RightArrayBinaryTree<TChar, TSpec>&treeStructure, FibreTreeVertieces)
{
    return treeStructure.treeVertieces;
}

template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type const &
getFibre(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, FibreTreeVertieces)
{
    return treeStructure.treeVertieces;
}

// This function returns the number of different entries in the wavelet tree structure.
template <typename TChar, typename TSpec>
inline unsigned length(RightArrayBinaryTree<TChar, TSpec> & tree)
{
    return length(tree.treeVertieces);
}

// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize>
inline void resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size)
{
    resize(treeStructure.treeVertieces, size);
}

// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize>
inline void resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size,
                   typename Value<typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeVertieces>::Type>::Type value)
{
    resize(treeStructure.treeVertieces, size, value);
}

template <typename TChar, typename TSpec>
inline bool open(
    RightArrayBinaryTree<TChar, TSpec> & treeStructure,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<TChar> minString;
    name = fileName;    append(name, ".rtv");   open(getFibre(treeStructure, FibreTreeVertieces()), toCString(name), openMode);
    name = fileName;    append(name, ".rtm");   open(minString, toCString(name), openMode);
    treeStructure.minCharValue = minString[0];
    return true;
}


template <typename TChar, typename TSpec>
inline bool open(
    RightArrayBinaryTree<TChar, TSpec> & treeStructure,
    const char * fileName)
{
    return open(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & treeStructure,
    const char * fileName,
    int openMode)
{
    String<char> name;

    String<TChar> minString;
    appendValue(minString, treeStructure.minCharValue);
    name = fileName;    append(name, ".rtv");   save(getFibre(treeStructure, FibreTreeVertieces()), toCString(name), openMode);
    name = fileName;    append(name, ".rtm");   save(minString, toCString(name), openMode);
    return true;
}

template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & treeStructure,
    const char * fileName)
{
    return save(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}


/**
.Function.getAlphabet
..summary:Determines the characters
..signature:computeRightArrayBinaryTree(waveletTreeStructure)
..param.waveletTreeStructure:A wavelet tree structure.
...type:Class.RightArrayBinaryTree
..remarks:This function determines all
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RightArrayBinaryTree<Dna5> waveletTreeStructure;
computeRightArrayBinaryTree(genome);
*/

/*template <typename TChar, typename TSpec>
inline String<TChar> getAlphabet(RightArrayBinaryTree<TChar, TSpec> & structure)
{
    typename Iterator<RightArrayBinaryTree<TChar, TSpec> >::Type it(structure, 0);
    String<TChar> alphabet;
    appendValue(alphabet, structure.minCharValue);
    getNexVertex(alphabet, it);
    return alphabet;
}*/


/*template <typename TPos>
unsigned getNumRequieredVertieces(Pair<TPos> const & range)
{
    return range.i2 - range.i1;
}*/

/*template <typename TBWT, typename TWaveletTreeSpec, typename TPrefixSumTable>
inline void computeTreeEntries(LfTable<WaveletTree<TBWT, TWaveletTreeSpec>, TPrefixSumTable> & lfTable)
{
    resize(lfTable.prefixSumTable, getAlphabetSize(lfTable.prefixSumTable), 0);
    computeTreeEntries(infix(lfTable.prefixSumTable, 0, getAlphabetSize(lfTable.prefixSumTable)));
}


template <typename TFreq, typename TChar, typename TSpec, typename TNumChildVertieces>
inline void computeTreeEntries(
    const String<Pair<TChar, TFreq> > & freq,
    RightArrayBinaryTree<TChar, TSpec> & structure,
    TNumChildVertieces & numChildVertieces)
{
    structure.minCharValue = freq[0].i1;
    typename Iterator<RightArrayBinaryTree<TChar, TSpec> >::Type iter(structure, 0);
    resize(structure, length(freq) - 1);
    computeTreeEntries(freq, iter, 0, (length(freq) - 1), numChildVertieces);
}*/

/*template <typename TChar, typename TSpec, typename TPos>
inline void addDollarVertex(RightArrayBinaryTree<TChar, TSpec> & structure, TPos minPos)
{
    //adjusting the array
    resize(structure, length(structure.treeVertieces) + 1);

    for (unsigned i = length(structure.treeVertieces) - 1; i > minPos; --i)
    {
        structure.treeVertieces[i] = structure.treeVertieces[i - 1];
    }

    for (unsigned i = 0; i < length(structure.treeVertieces); ++i)
    {
        if (structure.treeVertieces[i].i2 > minPos + 2)
            ++structure.treeVertieces[i].i2;
    }

    //the next line can be found in wavelet_tree.h in addDollarVertex
    //structure.treeVertieces[minPos].i2 = 1;
    structure.treeVertieces[minPos + 1].i1 = structure.treeVertieces[minPos].i1;
    structure.treeVertieces[minPos + 1].i2 = 0;


}*/

//////////////////////////////////////////////////////////////////////////////////////////////


/*template <typename TFreq, typename TChar, typename TSpec, typename TNumOfChildVertieces>
inline void computeTreeEntries(
    String<Pair<TChar, TFreq> > const & freq,
    Iter<RightArrayBinaryTree<TChar, TSpec>, Standard> & iter,
    unsigned lowerPosInFreq,
    unsigned upperPosInFreq,
    TNumOfChildVertieces & numChildVertieces)
{
    typedef TFreq TSize;

    unsigned oldLowerPosInFreq = lowerPosInFreq;
    unsigned oldUpperPosInFreq = upperPosInFreq;
    TSize leftSize = freq[lowerPosInFreq].i2;
    TSize rightSize = freq[upperPosInFreq].i2;

    if ((lowerPosInFreq == upperPosInFreq - 1) || (lowerPosInFreq == upperPosInFreq))
    {
        setCharacter(iter, freq[upperPosInFreq].i1);
        setVertexToLeaf(iter);
        ++numChildVertieces;
        return;
    }

    //determine the pivot element such that the number of 0 and 1 are as equal as possible
    while (upperPosInFreq - 1 > lowerPosInFreq)
    {
        if (leftSize < rightSize)
        {
            ++lowerPosInFreq;
            leftSize += freq[lowerPosInFreq].i2;
        }
        else
        {
            --upperPosInFreq;
            rightSize += freq[upperPosInFreq].i2;
        }
    }

    setCharacter(iter, freq[upperPosInFreq].i1);
    if (!leftSize && !rightSize)
    {
        setVertexToLeaf(iter);
        ++numChildVertieces;
        return;
    }

    //there must be either a left or right subtree or both
    typename Iterator<RightArrayBinaryTree<TChar, TSpec> >::Type iter2 = iter;

    //is there only one element on the left side?
    if (oldLowerPosInFreq == lowerPosInFreq)
    {
        setRightChildPos(iter2, iter2.position + numChildVertieces + 1);
        TNumOfChildVertieces numChildVertieces2 = 0;
        goRightChild(iter2);
        //++counter;
        computeTreeEntries(freq, iter2, upperPosInFreq, oldUpperPosInFreq, numChildVertieces2);
        numChildVertieces += numChildVertieces2 + 1;
        return;
    }
    else
    {
        setLeftChildPos(iter);
        goLeftChild(iter);
        computeTreeEntries(freq, iter, oldLowerPosInFreq, lowerPosInFreq, numChildVertieces);
    }

    //is there only one element on the right side?
    if (oldUpperPosInFreq == upperPosInFreq)
    {
        setLeftChildPos(iter2);
        ++numChildVertieces;
    }
    else
    {
        setRightChildPos(iter2, iter2.position + numChildVertieces + 1);
        TNumOfChildVertieces numChildVertieces2 = 0;
        goRightChild(iter2);
        computeTreeEntries(freq, iter2, upperPosInFreq, oldUpperPosInFreq, numChildVertieces2);
        numChildVertieces += numChildVertieces2 + 1;
        return;
    }
}*/



//template <typename TFreq, typename TChar, typename TSpec, typename TNumChildVertieces>
//inline void computeTreeEntries(
//		const String<TFreq> & freq,
//		RightArrayBinaryTree<TChar, TSpec> & structure,
//		TNumChildVertieces & numChildVertieces)
//{
//	String<Pair<TChar,TFreq> > newFreq;
//	for(TFreq i = 0; i < length(freq); ++i)
//		if(freq[i])
//			appendValue(newFreq, Pair<TChar, TFreq>((TChar)i, freq[i]));
//	structure.minCharValue = newFreq[0].i1;
//	typename Iterator<RightArrayBinaryTree<TChar, TSpec> >::Type iter(structure, 0);
//
//	resize(structure, length(newFreq) - 1);
//	computeTreeEntries(newFreq, iter, 0, (length(newFreq)-1), numChildVertieces);
//}

//template <typename TFreq, typename TChar, typename TSpec>
//inline void computeTreeEntries(RightArrayBinaryTree<TChar, TSpec> & structure,
//							   const String<TFreq> & freq)
//{
//	unsigned numChildVertieces = 0;
//	computeTreeEntries(freq, structure, numChildVertieces);
//}

/*template <typename TChar, typename TSpec>
inline bool open(
    RightArrayBinaryTree<TChar, TSpec> & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    open(getFibre(structure, FibreTreeVertieces()), toCString(name), openMode);
    return true;
}

template <typename TChar, typename TSpec>
inline bool open(RightArrayBinaryTree<TChar, TSpec> & structure,
                 const char * fileName)
{
    return open(structure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    save(getFibre(structure, FibreTreeVertieces()), toCString(name), openMode);
    return true;
}

template <typename TChar, typename TSpec>
inline bool save(
    RightArrayBinaryTree<TChar, TSpec> const & structure,
    const char * fileName)
{
    return save(structure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}*/

}

#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
