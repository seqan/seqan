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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREESTRUCTURE_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

template <typename TValue, typename TPointer, typename TSpec = void>
struct WaveletTreeStructure;
struct FibreTreeNodes_;

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================
typedef Tag<FibreTreeNodes_> const FibreTreeNodes;

template <typename TChar, typename TPointer, typename TSpec>
struct Fibre<WaveletTreeStructure<TChar, TPointer, TSpec>, FibreTreeNodes>
{
	//TODO
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
	//typedef unsigned TValue;
    typedef String<Pair<TChar, TPointer> > Type;
    typedef TPointer Value;
};

template <typename TChar, typename TPointer, typename TSpec>
struct WaveletTreeStructure
{
    typename Fibre<WaveletTreeStructure, FibreTreeNodes>::Type treeNodes;
    TChar minCharValue;

    WaveletTreeStructure() :
    	treeNodes(),
    	minCharValue()
    {};

    template <typename TFreqString>
    WaveletTreeStructure(TFreqString & freqString) :
    	treeNodes(),
    	minCharValue()
    {
        resize(treeNodes, length(freqString) - 1);
        computeTreeEntries(freqString,
                           *this);

    }

    inline WaveletTreeStructure & operator=(WaveletTreeStructure const & other)
    {
    	treeNodes = other.treeNodes;
    	minCharValue = other.minCharValue;
    	return *this;
    }

    inline bool operator==(const WaveletTreeStructure & b) const
    {
        return treeNodes == b.treeNodes;
    }
};

template <typename TValue>
inline void initStartTreeStructureValue(TValue & numChildNodes,
								 TValue & lowest,
								 TValue & highest,
								 TValue & posInTree,
								 TValue const sigmaSize)
{
	numChildNodes = 0;
	lowest = 0;
	highest = sigmaSize - 1;
	posInTree = 0;
}


// ==========================================================================
// Functions
// ==========================================================================

template <typename TChar, typename TPointer, typename TSpec>
inline typename Fibre<WaveletTreeStructure<TChar, TPointer, TSpec>, FibreTreeNodes>::Type &
getFibre(WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure, FibreTreeNodes)
{
    return treeStructure.treeNodes;
}

template <typename TChar, typename TPointer, typename TSpec>
inline typename Fibre<WaveletTreeStructure<TChar, TPointer, TSpec>, FibreTreeNodes>::Type const &
getFibre(WaveletTreeStructure<TChar, TPointer, TSpec> const & treeStructure, FibreTreeNodes)
{
    return treeStructure.treeNodes;
}

/*template< typename TBitString, typename TText, typename TSpec >
typename Fibre< WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec >, FibreSplitValues >::Type const &
getFibre(WaveletTree< TBitString, TSplitValue, TPosInSubTree, TSpec > const &tree, const FibreSplitValues)
{
    return tree.splitValues;
}*/

template <typename TChar, typename TPointer, typename TSpec>
inline unsigned length(WaveletTreeStructure<TChar, TPointer, TSpec> & tree)
{
    return length(tree.treeNodes);
}

template <typename TChar, typename TPointer, typename TSpec, typename TSize>
inline void resize(WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure, TSize size)
{
    resize(treeStructure.treeNodes, size);
}

template <typename TChar, typename TPointer, typename TSpec, typename TSize>
inline void resize(WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure, TSize size,
			typename Value<typename Fibre<WaveletTreeStructure<TChar, TPointer, TSpec>,FibreTreeNodes>::Type>::Type value)
{
    resize(treeStructure.treeNodes, size, value);
}

template <typename TChar, typename TPointer, typename TSpec>
inline void clear(WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure)
{
    clear(treeStructure.treeNodes);
}

template <typename TChar, typename TPointer, typename TSpec>
struct Iter<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard>
{
    TPointer position;
    const WaveletTreeStructure<TChar, TPointer, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(const WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure, const TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};



template <typename TChar, typename TPointer, typename TSpec>
struct Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard>
{
    TPointer position;
    WaveletTreeStructure<TChar, TPointer, TSpec> * waveletTreeStructure;

    template <typename TPos>
    Iter(WaveletTreeStructure<TChar, TPointer, TSpec> & treeStructure, TPos pos) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TChar, typename TPointer, typename TSpec>
struct Iterator<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard>
{
    typedef Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> Type;
};

template <typename TChar, typename TPointer, typename TSpec>
struct Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard>
{
    typedef Iter<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard> Type;
};

template <typename TChar, typename TPointer, typename TSpec>
struct Iterator<WaveletTreeStructure<TChar, TPointer, TSpec>, Rooted>:
    Iterator<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard>{};

template <typename TChar, typename TPointer, typename TSpec>
struct Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> const, Rooted>:
    Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard>{};


template <typename TChar, typename TPointer, typename TSpec>
struct Value<WaveletTreeStructure<TChar, TPointer, TSpec> >
{
    typedef Pair<TChar, TPointer> Type;
};

template <typename TChar, typename TPointer, typename TSpec>
struct Value<WaveletTreeStructure<TChar, TPointer, TSpec> const>
{
    typedef Pair<TChar, TPointer> const Type;
};

template <typename TChar, typename TPointer, typename TSpec>
struct Reference<WaveletTreeStructure<TChar, TPointer, TSpec> >
{
    typedef typename Value<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type Type;
};

template <typename TChar, typename TPointer, typename TSpec>
struct Reference<const WaveletTreeStructure<TChar, TPointer, TSpec> >
{
    typedef typename Value<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type const Type;
};


template <typename TChar, typename TPointer, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> const>::Type
begin(WaveletTreeStructure<TChar, TPointer, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type(waveletTreeStructure, 0);
}

template <typename TChar, typename TPointer, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> const>::Type
end(WaveletTreeStructure<TChar, TPointer, TSpec> const & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type(waveletTreeStructure, length(waveletTreeStructure.treeNodes));
}

template <typename TChar, typename TPointer, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type
begin(WaveletTreeStructure<TChar, TPointer, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type(waveletTreeStructure.treeNodes, 0);
}

template <typename TChar, typename TPointer, typename TSpec>
inline typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type
end(WaveletTreeStructure<TChar, TPointer, TSpec> & waveletTreeStructure)
{
    return typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type(waveletTreeStructure.treeNodes, length(waveletTreeStructure.treeNodes));
}

template <typename TChar, typename TPointer, typename TSpec, typename TPos>
inline void setPosition(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter, TPos pos)
{
    iter.position = pos;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getPosition(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    return iter.position;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getPosition(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    return iter.position;
}

template <typename TChar, typename TPointer, typename TSpec>
inline void setCharacter(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter,
                  TChar character)
{
    iter.waveletTreeStructure->treeNodes[iter.position].i1 = character;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TChar getCharacter(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    return iter.waveletTreeStructure->treeNodes[iter.position].i1;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TChar getCharacter(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    return iter.waveletTreeStructure->treeNodes[iter.position].i1;
}

template <typename TChar, typename TPointer, typename TSpec>
inline void setLeftChildPos(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = 2;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
    {
        return;
    }
    std::cerr << "ERROR: The right child has just been deleted!" << std::endl;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getLeftChildPos(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getLeftChildPos(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TChar, typename TPointer, typename TSpec>
inline void setRightChildPosOnly(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
}

template <typename TChar, typename TPointer, typename TSpec, typename TPos>
inline void setRightChildPos(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter, TPos rightChildPosition)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 0)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = 1;
        return;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
    {
        iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
        return ;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return;
    }
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = rightChildPosition + 2;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getRightChildPos(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return iter.position + 1;
    }
    return 0;
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getRightChildPos(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 > 2)
    {
        return iter.waveletTreeStructure->treeNodes[iter.position].i2 - 2;
    }
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
    {
        return iter.position + 1;
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TChar, typename TPointer, typename TSpec>
inline void setNodeToLeaf(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    if (iter.waveletTreeStructure->treeNodes[iter.position].i2 != 0)
    {
        std::cerr << "You just deleted ";
        if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 1)
        {
            std::cerr << "the right sub tree!" << std::endl;
        }
        if (iter.waveletTreeStructure->treeNodes[iter.position].i2 == 2)
        {
            std::cerr << "the left sub tree!" << std::endl;
        }
        else
        {
            std::cerr << "both sub trees!" << std::endl;
        }
    }
    iter.waveletTreeStructure->treeNodes[iter.position].i2 = 0;
}

template <typename TChar, typename TPointer, typename TSpec, typename TPos>
inline void addDollarNode(WaveletTreeStructure<TChar, TPointer, TSpec> & structure, TPos minPos)
{
	//adjusting the array
	resize(structure, length(structure.treeNodes) + 1);

	for(TPointer i = length(structure.treeNodes) - 1; i > minPos; --i)
	{
		structure.treeNodes[i] = structure.treeNodes[i - 1];
	}

	for(TPointer i = 0; i < length(structure.treeNodes); ++i)
	{
		if(structure.treeNodes[i].i2 > minPos + 2)
			++structure.treeNodes[i].i2;
	}

	//the next line can be found in wavelet_tree.h in addDollarNode
	//structure.treeNodes[minPos].i2 = 1;
	structure.treeNodes[minPos + 1].i1 = structure.treeNodes[minPos].i1;
	structure.treeNodes[minPos + 1].i2 = 0;


}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TChar, typename TPointer, typename TSpec>
inline void goLeft(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    iter.position = getLeftChildPos(iter);
}

template <typename TChar, typename TPointer, typename TSpec>
inline void goLeft(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    iter.position = getLeftChildPos(iter);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TChar, typename TPointer, typename TSpec>
inline void goRight(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    iter.position = getRightChildPos(iter);
}

template <typename TChar, typename TPointer, typename TSpec>
inline void goRight(Iter<const WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter)
{
    iter.position = getRightChildPos(iter);
}

//////////////////////////////////////////////////////////////////////////////////////////////
template <typename TText, typename TIter>
inline void getNexNode(TText & alphabet, TIter & it)
{
	TIter copyIt = it;
	if(getLeftChildPos(it))
	{
		goLeft(it);
		getNexNode(alphabet, it);
	}
	appendValue(alphabet, getCharacter(copyIt));
	if(getRightChildPos(it))
	{
		goRight(copyIt);
		getNexNode(alphabet, copyIt);
	}
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool isLeaf(Iter<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard> & iter){
	return !((*iter.waveletTreeStructure).treeNodes[iter.position].i2);
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool isLeaf(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter){
	//std::cerr << (int)iter.position << " " << (int)(*iter.waveletTreeStructure).treeNodes[iter.position].i2 << std::endl;
	return !((*iter.waveletTreeStructure).treeNodes[iter.position].i2);
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getNodePosition(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter, TChar character){

	//TPointer pos = 0;
	while(!isLeaf(iter))
	{
		if(character < getCharacter(iter))
			goLeft(iter);
		else
			goRight(iter);
	}
	return getPosition(iter);
}

template <typename TChar, typename TPointer, typename TSpec>
inline TPointer getNodePosition(Iter<WaveletTreeStructure<TChar, TPointer, TSpec> const, Standard> & iter, TChar character){

	while(!isLeaf(iter))
		{
			if(character < getCharacter(iter))
				goLeft(iter);
			else
				goRight(iter);
		}
		return getPosition(iter);
}

template <typename TChar, typename TPointer, typename TSpec>
inline String<TChar> getAlphabet(WaveletTreeStructure<TChar, TPointer, TSpec> & structure)
{
	typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type it(structure, 0);
	String<TChar> alphabet;
	appendValue(alphabet, structure.minCharValue);
	getNexNode(alphabet, it);
	//std::cerr << alphabet << " " << length(alphabet) << std::endl;
	return alphabet;
}


template <typename TChar, typename TPointer, typename TSpec, typename TPosInSubTree>
inline void printTreeLevel(WaveletTreeStructure<TChar, TPointer, TSpec> & tree, TPosInSubTree posInTree)
{
    TPosInSubTree leftChild = tree.treeNodes[posInTree].i2;
    TPosInSubTree rightChild = tree.treeNodes[posInTree].i3;
    std::cout << "(" << posInTree << ": " << tree.treeNodes[posInTree].i1 << ", " << leftChild << ", " << rightChild << ") | ";
    if (leftChild)
    {
        printTreeLevel(tree, leftChild);
    }
    if (rightChild)
    {
        printTreeLevel(tree, rightChild);
    }
}

template <typename TChar, typename TPointer, typename TSpec>
inline void printTree(WaveletTreeStructure<TChar, TPointer, TSpec> & tree)
{
    printTreeLevel(tree, 0);
}

template <typename TChar, typename TPointer, typename TSpec, typename TString>
inline void writeGraphImpl(Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter, TString name)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter2 = iter;
    std::ofstream stream(toCString(name), std::ios::app);
    TPointer pos = getLeftChildPos(iter);
    if (pos)
    {
        stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << ordValue(iter.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
        goLeft(iter);
        writeGraphImpl(iter, name);
    }
    else
    {
        stream << ordValue(iter.waveletTreeStructure->treeNodes[iter.position].i1) << " -> " << "leave1" << ordValue(iter.position) << ";" << std::endl;
    }

    pos = getRightChildPos(iter2);
    if (pos)
    {
        stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << ordValue(iter2.waveletTreeStructure->treeNodes[pos].i1) << ";" << std::endl;
        goRight(iter2);
        writeGraphImpl(iter2, name);
    }
    else
    {
        stream << ordValue(iter2.waveletTreeStructure->treeNodes[iter2.position].i1) << " -> " << "leave2" << ordValue(iter2.position) << ";" << std::endl;
    }
    stream.close();
}

template <typename TChar, typename TPointer, typename TSpec>
inline void writeGraph(WaveletTreeStructure<TChar, TPointer, TSpec> & tree)
{

    typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter(tree, 0);

    String<char> name = "testfile.dot";
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}

template <typename TChar, typename TPointer, typename TSpec, typename TString>
inline void writeGraph(WaveletTreeStructure<TChar, TPointer, TSpec> & tree, TString name)
{

    typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter(tree, 0);
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}

template <typename TStringLengthString, typename TPrefixSumTable, typename TAlphabetString, typename TIter>//, typename TTreeSplitValue>
inline void computeSingleStringLengthFromTree(TStringLengthString & lengthString,
                                       TPrefixSumTable & prefixSumTable,
                                       TAlphabetString & alphabet,
                                       TIter & iter,
                                       unsigned lowerBound,
                                       unsigned upperBound)
//										TreeSplitValue lowerBound,
//                                       TTreeSplitValue upperBound)
{
	typedef typename Value<TAlphabetString>::Type TChar;

    lengthString[iter.position] = prefixSumTable[upperBound + 1] - prefixSumTable[lowerBound];
    TIter iter2 = iter;

    TChar pivot = getCharacter(iter);

    unsigned newSplit, newSplit2;
    for(unsigned i = 0; i < length(alphabet); ++i)
    {
    	if(alphabet[i] == pivot)
    	{
    		newSplit = i - 1;
    		newSplit2 = i;
    	}
    }

    goLeft(iter);
    if (getPosition(iter))
    {
        computeSingleStringLengthFromTree(lengthString,
                                          prefixSumTable,
                                          alphabet,
                                          iter,
                                          lowerBound,
                                          newSplit);
    }

    //newSplit = pivot;
    goRight(iter2);
    if (getPosition(iter2))
    {
        computeSingleStringLengthFromTree(lengthString,
                                          prefixSumTable,
                                          alphabet,
                                          iter2,
                                          newSplit2,
                                          upperBound);
    }
}

template <typename TStringLengthString, typename TPrefixSumTable, typename TChar, typename TPointer, typename TSpec>
inline void computeStringLengthFromTree(TStringLengthString & lengthString,
                                 TPrefixSumTable & prefixSumTable,
                                 WaveletTreeStructure<TChar, TPointer, TSpec> & structure)
{
    resize(lengthString, length(structure.treeNodes), 0);
    typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter(structure, 0);

    String<TChar> alphabet = getAlphabet(structure);

    computeSingleStringLengthFromTree(lengthString,
                                      prefixSumTable,
                                      alphabet,
                                      iter,
                                      (TPointer)0,
                                      (TPointer)(length(prefixSumTable) - 2));
}


template <typename TFreq, typename TChar, typename TPointer, typename TSpec, typename TNumOfChildNodes>
inline void computeTreeEntries(
    String<Pair<TChar, TFreq> > const & freq,
    Iter<WaveletTreeStructure<TChar, TPointer, TSpec>, Standard> & iter,
    TPointer lowerPosInFreq,
    TPointer upperPosInFreq,
    TNumOfChildNodes & numChildNodes)
{
    typedef TFreq TSize;

	TPointer oldLowerPosInFreq = lowerPosInFreq;
	TPointer oldUpperPosInFreq = upperPosInFreq;
    TSize leftSize = freq[lowerPosInFreq].i2;
    TSize rightSize = freq[upperPosInFreq].i2;

    if (lowerPosInFreq == upperPosInFreq - 1)
    {
        setCharacter(iter, freq[upperPosInFreq].i1);
        setNodeToLeaf(iter);
        ++numChildNodes;
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
    	setNodeToLeaf(iter);
    	++numChildNodes;
    	return;
    }

    //there must be either a left or right subtree or both
    typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter2 = iter;

    //is there only one element on the left side?
    if (oldLowerPosInFreq == lowerPosInFreq)
    {
   		setRightChildPos(iter2, iter2.position + numChildNodes + 1);
   		TNumOfChildNodes numChildNodes2 = 0;
   		goRight(iter2);
   		//++counter;
    	computeTreeEntries(freq, iter2, upperPosInFreq, oldUpperPosInFreq, numChildNodes2);
    	numChildNodes += numChildNodes2 + 1;
    	return;
    }
    else
    {
        setLeftChildPos(iter);
        goLeft(iter);
        computeTreeEntries(freq, iter, oldLowerPosInFreq, lowerPosInFreq, numChildNodes);
    }

    //is there only one element on the right side?
    if (oldUpperPosInFreq == upperPosInFreq)
    {
        setLeftChildPos(iter2);
        ++numChildNodes;
    }
    else
    {
    	setRightChildPos(iter2, iter2.position + numChildNodes + 1);
        TNumOfChildNodes numChildNodes2 = 0;
        goRight(iter2);
        computeTreeEntries(freq, iter2, upperPosInFreq, oldUpperPosInFreq, numChildNodes2);
        numChildNodes += numChildNodes2 + 1;
        return;
    }
}



template <typename TFreq, typename TChar, typename TPointer, typename TSpec, typename TNumChildNodes>
inline void computeTreeEntries(
    const String<Pair<TChar, TFreq> > & freq,
    WaveletTreeStructure<TChar, TPointer, TSpec> & structure,
    TNumChildNodes & numChildNodes)
{
	structure.minCharValue = freq[0].i1;
	typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter(structure, 0);
	resize(structure, length(freq) - 1);
	computeTreeEntries(freq, iter, (TPointer)0, (TPointer)(length(freq)-1), numChildNodes);
}

template <typename TFreq, typename TChar, typename TPointer, typename TSpec, typename TNumChildNodes>
inline void computeTreeEntries(
    const String<TFreq> & freq,
    WaveletTreeStructure<TChar, TPointer, TSpec> & structure,
    TNumChildNodes & numChildNodes)
{
	String<Pair<TChar,TFreq> > newFreq;
	for(TFreq i = 0; i < length(freq); ++i)
		if(freq[i])
			appendValue(newFreq, Pair<TChar, TFreq>((TChar)i, freq[i]));
	structure.minCharValue = newFreq[0].i1;
	typename Iterator<WaveletTreeStructure<TChar, TPointer, TSpec> >::Type iter(structure, 0);

	resize(structure, length(newFreq) - 1);
	computeTreeEntries(newFreq, iter, (TPointer)0, (TPointer)(length(newFreq)-1), numChildNodes);
}

template <typename TFreq, typename TChar, typename TPointer, typename TSpec>
inline void computeTreeEntries(
    const String<TFreq> & freq,
    WaveletTreeStructure<TChar, TPointer, TSpec> & structure)
{

	TPointer numChildNodes = 0;
	computeTreeEntries(freq, structure, numChildNodes);
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool open(
    WaveletTreeStructure<TChar, TPointer, TSpec> & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    open(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
    return true;
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool open(WaveletTreeStructure<TChar, TPointer, TSpec> & structure,
				 const char * fileName)
{
    return open(structure, fileName, DefaultOpenMode<WaveletTreeStructure<TChar, TPointer, TSpec> >::VALUE);
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool save(
    WaveletTreeStructure<TChar, TPointer, TSpec> const & structure,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".treestruct");    save(getFibre(structure, FibreTreeNodes()), toCString(name), openMode);
    return true;
}

template <typename TChar, typename TPointer, typename TSpec>
inline bool save(
    WaveletTreeStructure<TChar, TPointer, TSpec> const & structure,
    const char * fileName)
{
    return save(structure, fileName, DefaultOpenMode<WaveletTreeStructure<TChar, TPointer, TSpec> >::VALUE);
}

}

#endif
