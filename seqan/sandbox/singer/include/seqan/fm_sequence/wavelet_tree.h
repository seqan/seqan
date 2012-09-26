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

#ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_
#define SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================

struct DollarSubstituted_;
struct MultiDollarSubstituted_;
struct FibreBitStrings_;
struct FibreSplitValues_;
struct FibreDollarPositions_;
struct AlphabetExternalChar_;
struct FibreTreeNodes_;

template <typename TText, typename TSpec>
struct WaveletTree;

// ==========================================================================
//Tags, Classes, Enums
// ==========================================================================

typedef Tag<DollarSubstituted_> DollarSubstituted;
typedef Tag<MultiDollarSubstituted_> MultiDollarSubstituted;
typedef Tag<AlphabetExternalChar_> const AlphabetExternalChar;
typedef Tag<FibreBitStrings_> const FibreBitStrings;
typedef Tag<FibreSplitValues_> const FibreSplitValues;
typedef Tag<FibreDollarPositions_> const FibreDollarPositions;
typedef Tag<FibreTreeNodes_> const FibreTreeNodes;

template <typename TText, typename TSpec = void>
struct WaveletTree
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type TBitStrings;
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreSplitValues>::Type TWaveletTreeStructure;

    TBitStrings bitStrings;
    TWaveletTreeStructure splitValues;

    WaveletTree(){}

    WaveletTree(TText const & text)
    {
        createWaveletTree(* this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable)
    {
        createWaveletTree(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        createWaveletTree(*this, text, freqTable, prefixSumTable);
    }

    inline WaveletTree & operator=(WaveletTree const & other)
    {
    	bitStrings = other.bitStrings;
    	splitValues = other.splitValues;
    	return *this;
    }

    inline bool operator==(const WaveletTree & b) const
    {
        typedef typename Size<TText>::Type                            TSize;
        bool test = true;
        if (length(bitStrings) == length(b.bitStrings))
        {
            for (TSize i = 0; i < length(bitStrings); ++i)
            {
                if (!(bitStrings[i] == b.bitStrings[i]))
                {
                    test = false;
                }
            }
        }
        else
        {
            test = false;
        }

        return test &&
               splitValues == b.splitValues;
    }


};

template <typename TText>
struct WaveletTree<TText, DollarSubstituted>
{
    typedef typename Fibre<WaveletTree, FibreBitStrings>::Type    TBitStrings;
    typedef typename Size<TText>::Type                            TSize;
    typedef typename Fibre<WaveletTree<TText, DollarSubstituted>, FibreSplitValues>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TBitStrings                         bitStrings;
    TWaveletTreeStructure       	splitValues;
    TSize                               dollarPosition;
    TChar                              dollarSub;

    WaveletTree(){}

    WaveletTree(TText const & text)
    {
        dollarPosition = (TSize) - 1;
        dollarSub = 0;
        createWaveletTree(*this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable)
    {
        dollarPosition = (TSize) - 1;
        dollarSub = 0;
        createWaveletTree(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        dollarPosition = (TSize) - 1;
        dollarSub = 0;
        createWaveletTree(*this, text, freqTable, prefixSumTable);
    }

    inline bool operator==(const WaveletTree & b) const
    {

        typedef typename Size<TText>::Type                            TSize;
        bool test = true;
        if (length(bitStrings) == length(b.bitStrings))
        {
            for (TSize i = 0; i < length(bitStrings); ++i)
            {
                if (!(bitStrings[i] == b.bitStrings[i]))
                {
                    test = false;
                }
            }
        }
        else
        {
            test = false;
        }

        return test &&
               splitValues == b.splitValues &&
               dollarPosition == b.dollarPosition &&
               dollarSub == b.dollarSub;
    }

};

template <typename TText>
struct WaveletTree<TText, MultiDollarSubstituted>
{
    typedef typename Fibre<WaveletTree, FibreBitStrings>::Type                                TBitStrings;
    typedef typename Size<TText>::Type                                                        TSize;
    typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreSplitValues>::Type TWaveletTreeStructure;
    typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>::Type TDollarString;
    //typedef typename Fibre<TWaveletTreeStructure, FibreTreeNodes>::Value TValue;


    typedef typename Fibre<TWaveletTreeStructure, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TBitStrings                             bitStrings;
    TWaveletTreeStructure		          splitValues;

    TDollarString                           dollarPosition;
    TChar                                  dollarSub;

    WaveletTree(){}

    WaveletTree(TText const & text)
    {
        dollarSub = 0;
        resize(dollarPosition, length(text));
        createWaveletTree(*this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable)
    {
        dollarSub = 0;
        resize(dollarPosition, length(text));
        createWaveletTree(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        dollarSub = 0;
        resize(dollarPosition, length(text));
        createWaveletTree(*this, text, freqTable, prefixSumTable);
    }

    inline bool operator==(const WaveletTree & b) const
    {
        typedef typename Size<TText>::Type TSize;
        bool test = true;
        if (length(bitStrings) == length(b.bitStrings))
        {
            for (TSize i = 0; i < length(bitStrings); ++i)
            {
                if (!(bitStrings[i] == b.bitStrings[i]))
                {
                    test = false;
                }
            }
        }
        else
        {
            test = false;
        }

        return test &&
               splitValues == b.splitValues &&
               dollarPosition == b.dollarPosition &&
               dollarSub == b.dollarSub;
    }

};

// ==========================================================================
//Metafunctions
// ==========================================================================
template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

//template <typename TText, typename TSpec>
//struct Fibre<WaveletTree<TText, TSpec>, FibreSplitValues>
//{
//    typedef WaveletTreeStructure<TText> Type;
//};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreSplitValues>
{
    typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TPointer;
    typedef WaveletTreeStructure<typename Value<TText>::Type, TPointer, void> Type;
};

template <typename TText>
struct Fibre<WaveletTree<TText, AlphabetExternalChar>, FibreSplitValues>
{
    typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE + 1>::Type TPointer;
    typedef WaveletTreeStructure<typename Value<TText>::Type, TPointer, void> Type;
};

template <typename TText>
struct Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>
{
    typedef typename Value<typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreBitStrings>::Type>::Type Type;
};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> >
{
	typedef typename Value<TText>::Type Type;
};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> const>
{
	typedef typename Value<TText>::Type const Type;
};

// ==========================================================================
//Functions
// ==========================================================================
template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type &
getFibre(WaveletTree<TText, TSpec> & tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type const &
getFibre(const WaveletTree<TText, TSpec> & tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreSplitValues>::Type &
getFibre(WaveletTree<TText, TSpec> & tree, FibreSplitValues)
{
    return tree.splitValues;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreSplitValues>::Type const &
getFibre(WaveletTree<TText, TSpec> const & tree, const FibreSplitValues)
{
    return tree.splitValues;
}

template <typename TText>
inline typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>::Type &
getFibre(WaveletTree<TText, MultiDollarSubstituted> & tree, const FibreDollarPositions)
{
    return tree.dollarPosition;
}

template <typename TText>
inline typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>::Type const &
getFibre(WaveletTree<TText, MultiDollarSubstituted> const & tree, const FibreDollarPositions)
{
    return tree.dollarPosition;
}

template <typename TText, typename TSpec>
inline void clear(WaveletTree<TText, TSpec> & tree)
{
    for (unsigned i = 0; i < length(tree.bitStrings); ++i)
    {
        clear(tree.bitStrings[i]);
    }
    resize(tree.bitStrings, 0);
    clear(tree.splitValues);
}

template <typename TText, typename TSpec>
inline unsigned getNumNodes(WaveletTree<TText, TSpec> & tree)
{
	return length(tree.splitValues.treeNodes);
}


template <
    typename TText,
    typename TWaveletTreeSpec,
    typename TCharIn,
    typename TPos>
inline unsigned getOccImpl(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                           const TCharIn character,
                           const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TPos sum = pos + 1;
    TPointer treePos = 0;


    typename Iterator<const TSplitValues>::Type it(tree.splitValues, treePos);
    FibreSplitValues tag = FibreSplitValues();
    TChar charInTree = tree.splitValues.minCharValue;
    do
    {
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        if (character < getFibre(tree, tag).treeNodes[treePos].i1)
        {
            sum -= addValue;
            goLeft(it);
        }
        else
        {
        	charInTree = getCharacter(it);
            sum = addValue;
            goRight(it);
        }
        treePos = getPosition(it);
    }
    while ((bool)treePos && (bool)sum);

    if(character == charInTree)
    	return sum;
    return 0;
}

template <typename TText, typename TChar, typename TPos, typename TWaveletTreeSpec>
inline unsigned getOcc(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                       const TChar character,
                       const TPos pos)
{
    //typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
   // typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
    //TValue ordChar = ordValue(character);
    return getOccImpl(tree, character, pos);
}

template <typename TText, typename TChar, typename TPos>
inline unsigned getOcc(const WaveletTree<TText, DollarSubstituted> & tree,
                       const TChar character,
                       const TPos pos)
{
//	typedef typename Fibre<WaveletTree<TText, DollarSubstituted>, FibreSplitValues>::Type TSplitValues;
//    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
//
//	TValue ordChar = ordValue(character);
    unsigned occ = getOccImpl(tree, character, pos);
    //static_cast<Nothing>(character);
    if (character == tree.dollarSub && pos >= tree.dollarPosition)
    {
        return occ - 1;
    }
    return occ;
}

template <typename TText, typename TChar, typename TPos>
inline unsigned getOcc(const WaveletTree<TText, MultiDollarSubstituted> & tree,
                       const TChar character,
                       const TPos pos)
{
	typedef typename Fibre<WaveletTree<TText, DollarSubstituted>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;

	TValue ordChar = ordValue(character);
    unsigned occ = getOccImpl(tree, ordChar, pos);
    unsigned numDollar = getRank(tree.dollarPosition, pos);
    if (ordChar == tree.dollarSub)
    {
        return occ - numDollar;
    }
    return occ;
}

//template <typename TText, typename TWaveletTreeSpec, typename TPos>
//inline typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type
//getCharacter(const WaveletTree<TText, TWaveletTreeSpec> & tree,
//             const TPos pos)
//{
//    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
//    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
//
//    TPos sum = pos + 1;
//    TValue treePos = 0;
//
//    typename Iterator<const TSplitValues>::Type iter(tree.splitValues, treePos);
//    bool direction;
//    TValue character;
//    do
//    {
//        direction = getBit(tree.bitStrings[treePos], sum - 1);
//        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
//        if (direction)
//        {
//            character = getCharacter(iter) + 1;
//            sum = addValue;
//            goRight(iter);
//        }
//        else
//        {
//            character = getCharacter(iter);
//            sum -= addValue;
//            goLeft(iter);
//        }
//        treePos = getPosition(iter);
//    }
//    while (treePos);
//
//    return character;
//}

template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline typename Value<TText>::Type
getCharacter(const WaveletTree<TText, TWaveletTreeSpec> & tree,
             const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;

    TPos sum = pos + 1;
    TPointer treePos = 0;

    //std::cerr << "getChar0" << std::endl;

    typename Iterator<const TSplitValues>::Type iter(tree.splitValues, treePos);
    bool direction;
    TChar character = tree.splitValues.minCharValue;
    //std::cerr << "getChar1" << std::endl;
    do
    {
       // std::cerr << "getChar2 " << (int)treePos << " " << length(tree.bitStrings) << " " << pos << " "  << sum -1 << std::endl;
        direction = getBit(tree.bitStrings[treePos], sum - 1);
        //std::cerr << "getChar2.1" << std::endl;
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        //std::cerr << "getChar2.2" << std::endl;
        if (direction)
        {
            character = getCharacter(iter);// + 1;
            sum = addValue;
            goRight(iter);
        }
        else
        {
            //character = getCharacter(iter);
            sum -= addValue;
            goLeft(iter);
        }
        //std::cerr << "getChar3" << std::endl;
        treePos = getPosition(iter);
        //std::cerr << "getChar4" << std::endl;
    }
    while (treePos);

    return character;
}

template <typename TText, typename TWaveletTreeSpec>
inline typename Value<TText>::Type
getDollarSub(const WaveletTree<TText, TWaveletTreeSpec> & tree)
{
    return tree.dollarSub;
}

template <typename TText, typename TWaveletTreeSpec>
inline TText getAlphabet(WaveletTree<TText, TWaveletTreeSpec> & tree)
{
	getAlphabet(tree.splitValues);
}

inline unsigned getTreeLevel(unsigned treePosition)
{
    return floor(log(treePosition + 1) / log(2));
}

template <typename TText, typename TPos>
inline void setDollarPosition(WaveletTree<TText, DollarSubstituted> & tree,
                              TPos position)
{
    tree.dollarPosition = position;
}

template <typename TText, typename TPos>
inline void setDollarPosition(WaveletTree<TText, MultiDollarSubstituted> & tree,
                              TPos position)
{
    setBit(tree.dollarPosition, position, 1);
}

template <typename TText, typename TWaveletTreeSpec>
inline void fillWaveletTree(
    WaveletTree<TText, TWaveletTreeSpec> & tree,                //the set of bit strings to be filled
    const TText & text)                                         //the original string
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreRankSupportBitString>::Type TFibreBitString;
    typedef typename Size<TFibreBitString>::Type TSize;
	typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;


    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TSplitValues >::Type iter(tree.splitValues, 0);
        bool bit;

        do
        {
            if (value(text, i) < getCharacter(iter))
            {
                bit = 0;
                //appendBitOnly(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
                appendBit(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
                goLeft(iter);
            }
            else
            {
                bit = 1;
                //appendBitOnly(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
                appendBit(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
                goRight(iter);
            }
        }
        while (getPosition(iter) != 0);
    }

    TFibreRankSupportBitStrings & bitStrings = getFibre(tree, FibreBitStrings());
    for (TSize i = 0; i < length(bitStrings); ++i)
    {
        TFibreRankSupportBitString & temp = bitStrings[i];
        completeRankSupportBitString(temp);
    }
}
//this function is used to fill all bit strings in a wavelet tree AND counting
template <typename TText, typename TWaveletTreeSpec, typename TDollarPosition>
inline void fillWaveletTree(
    WaveletTree<TText, TWaveletTreeSpec> & tree,                //the set of bit strings to be filled
    const TText & text,
    TDollarPosition dollarPosition)                                         //the original string
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreRankSupportBitString>::Type TFibreBitString;
    typedef typename Size<TFibreBitString>::Type TSize;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;

    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TSplitValues>::Type iter(tree.splitValues, 0);
        bool bit;

        do
        {
        	//std::cerr << value(text, i) << std::endl;
        	//std::cerr << getCharacter(iter) << std::endl;
        	if(i == dollarPosition)
        	{
        		bit = 1;
        		//appendBitOnly(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		appendBit(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		goRight(iter);
        	}
        	else if (value(text, i) < getCharacter(iter))
        	{
        		bit = 0;
        		//appendBitOnly(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		appendBit(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		goLeft(iter);


        	}
        	else
        	{
        		bit = 1;
        		//appendBitOnly(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		appendBit(getFibre(tree, FibreBitStrings())[getPosition(iter)], bit);
        		goRight(iter);
        	}
        }
        while (getPosition(iter) != 0);
    }

    TFibreRankSupportBitStrings & bitStrings = getFibre(tree, FibreBitStrings());
    for (TSize i = 0; i < length(bitStrings); ++i)
    {
        TFibreRankSupportBitString & temp = bitStrings[i];
        completeRankSupportBitString(temp);
    }
}

//this function is used to fill all bit strings in a wavelet tree
template <typename TText, typename TWaveletTreeSpec>
inline void fillWaveletTree(
    WaveletTree<TText, TWaveletTreeSpec> & tree,                //the set of bit strings to be filled
    typename Iterator<WaveletTree<TText, TWaveletTreeSpec> >::Type iter,
    const TText & text,                                 //the original string
    const typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Value lowerBound,
    const typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Value upperBound)
{

    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Value TValue;
    typename Iterator<WaveletTree<TText, TWaveletTreeSpec> >::Type iterRight = iter;

    fillBitString(
        tree,
        iter,
        text,
        lowerBound,
        upperBound);

    TValue leftUpperBound = getCharacter(iter);
    goLeft(iter);
    if (getPosition(iter))
    {
        fillWaveletTree(
            tree,
            iter,
            text,
            lowerBound,
            leftUpperBound);
    }

    TValue rightLowerBound = ordValue(getCharacter(iterRight)) + 1;
    goRight(iterRight);
    if (getPosition(iterRight))
    {
        fillWaveletTree(
            tree,
            iterRight,
            text,
            rightLowerBound,
            upperBound);
    }
}

//template <typename TText, typename TWaveletTreeSpec, typename TDollarPos>
//inline void completeDollarNode(WaveletTree<TText, TWaveletTreeSpec> & tree, unsigned posInTree, TDollarPos dollarPosInLeaf, bool bit)
//{
//	setBit(tree.bitStrings[posInTree + 1], dollarPosInLeaf, 1);
//}
//
//template <typename TText, typename TWaveletTreeSpec, typename TDollarPos>
//inline void completeDollarNode(WaveletTree<TText, TWaveletTreeSpec> & tree, unsigned posInTree, RankSupportBitString<void> & dollarPosInLeaf, bool bit)
//{
//	for(unsigned i = 0; i < length(dollarPosInLeaf)
//	setBit(tree.bitStrings[posInTree + 1], dollarPosInLeaf, 1);
//}

template <typename TText, typename TWaveletTreeSpec, typename TDollarChar, typename TDollarPos>
inline void addDollarNode(WaveletTree<TText, TWaveletTreeSpec> & tree, TDollarChar dollarSub, TDollarPos dollarPos)
{
	typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
	typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
	typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
	typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
	typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;
	typedef typename Size<typename Value<TText>::Type>::Type TSize;

	TPointer minPos = getNodePosition(tree, dollarSub);

	addDollarNode(tree.splitValues, minPos);

	resize(tree.bitStrings, length(tree.bitStrings) + 1);
	for(TPointer i = length(tree.bitStrings) - 1; i > minPos; --i)
	{
		tree.bitStrings[i] = tree.bitStrings[i - 1];
	}

	//TChar freqChar = (TChar)minPos;
	TSize newNodeLegth;
	TSize dollarPosInLeaf = getOcc(tree, dollarSub, dollarPos) - 1;

	//do we need to add 0 or 1
	if(dollarSub < tree.splitValues.treeNodes[minPos].i1)
	{
		tree.splitValues.treeNodes[minPos].i2 = 2;
		newNodeLegth = length(tree.bitStrings[minPos-1]) - getRank(tree.bitStrings[minPos-1], length(tree.bitStrings[minPos-1]) - 1);
		clear(tree.bitStrings[minPos + 1]);
		resize(tree.bitStrings[minPos + 1], newNodeLegth, 0);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos + 1], i, 0);
		setBit(tree.bitStrings[minPos + 1], dollarPosInLeaf, 1);


	}
	else
	{
		tree.splitValues.treeNodes[minPos].i2 = 1;
		newNodeLegth = getRank(tree.bitStrings[minPos], length(tree.bitStrings[minPos]) - 1);
		clear(tree.bitStrings[minPos + 1]);
		resize(tree.bitStrings[minPos + 1], newNodeLegth, 0);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos + 1], i, 1);
		setBit(tree.bitStrings[minPos + 1], dollarPosInLeaf, 0);

	}
	completeRankSupportBitString(tree.bitStrings[minPos + 1]);

}

template <typename TText, typename TWaveletTreeSpec, typename TDollarChar>
inline void addDollarNode(WaveletTree<TText, TWaveletTreeSpec> & tree, TDollarChar dollarSub, RankSupportBitString<void> & dollarPos)
{
	typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
	typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
	typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
	typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
	typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;
	typedef typename Size<typename Value<TText>::Type>::Type TSize;

	TPointer minPos = getNodePosition(tree, dollarSub);

	addDollarNode(tree.splitValues, minPos);

	resize(tree.bitStrings, length(tree.bitStrings) + 1);
	for(TPointer i = length(tree.bitStrings) - 1; i > minPos; --i)
	{
		tree.bitStrings[i] = tree.bitStrings[i - 1];
	}

	//TChar freqChar = (TChar)minPos;
	TSize newNodeLegth;
	String<unsigned> dollarPosInLeaf;

	for(unsigned i = 0; i < length(dollarPos); ++i)
	{
		if(getBit(dollarPos, i))
		{
			appendValue(dollarPosInLeaf, getOcc(tree, dollarSub, i)-1);
//			std::cerr << getOcc(tree, dollarSub, i) << std::endl;
//			char c;
//			std::cin>>c;
		}

	}

	std::cerr << "length(dollarPosInLeaf): " << length(dollarPosInLeaf) << std::endl;

	//do we need to add 0 or 1
	if(dollarSub < tree.splitValues.treeNodes[minPos].i1)
	{
		tree.splitValues.treeNodes[minPos].i2 = 2;
		newNodeLegth = length(tree.bitStrings[minPos]) - getRank(tree.bitStrings[minPos], length(tree.bitStrings[minPos]) - 1);
		clear(tree.bitStrings[minPos + 1]);
		resize(tree.bitStrings[minPos + 1], newNodeLegth, 0);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos + 1], i, 0);

		for(unsigned i = 0; i < length(dollarPosInLeaf); ++i)
		{
			setBit(tree.bitStrings[minPos + 1], dollarPosInLeaf[i], 1);
		}

	}
	else
	{
		tree.splitValues.treeNodes[minPos].i2 = 1;
		newNodeLegth = getRank(tree.bitStrings[minPos], length(tree.bitStrings[minPos]) - 1);
		clear(tree.bitStrings[minPos + 1]);
		resize(tree.bitStrings[minPos + 1], newNodeLegth, 0);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos + 1], i, 1);

		for(unsigned i = 0; i < length(dollarPosInLeaf); ++i)
		{
			setBit(tree.bitStrings[minPos + 1], dollarPosInLeaf[i], 0);
		}
	}
	completeRankSupportBitString(tree.bitStrings[minPos + 1]);

}

template <typename TText, typename TWaveletTreeSpec, typename TFreq>
inline void addDollarNode(WaveletTree<TText, TWaveletTreeSpec> & tree, String<TFreq> & freq)
{
	typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
	typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
	typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
	typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
	typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;
	typedef typename Size<typename Value<TText>::Type>::Type TSize;

	TPointer minPos = 0;
	TPointer i = 0;
	for(; i < length(freq); ++i)
	{
		if(freq[i] > 0)
		{
			minPos = i;
			break;
		}
	}

	for(; i < length(freq) - 1; ++i)
	{
		if(freq[i] < freq[minPos] && freq[i] > 0)
			minPos = i;
	}

	minPos = getNodePosition(tree, (TChar)minPos);

	addDollarNode(tree.splitValues, minPos);

	resize(tree.bitStrings, length(tree.bitStrings) + 1);
	for(TPointer i = length(tree.bitStrings) - 1; i > minPos; --i)
	{
		tree.bitStrings[i] = tree.bitStrings[i - 1];
	}

	TChar freqChar = (TChar)minPos;
	TPointer newNodeLegth;
	//do we need to add 0 or 1
	if(freqChar < tree.splitValues.treeNodes[minPos].i1)
	{
		newNodeLegth = length(tree.bitStrings[minPos - 1]) - getRank(tree.bitStrings[minPos - 1], length(tree.bitStrings[minPos - 1]) - 1);
		resize(tree.bitStrings[minPos], newNodeLegth);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos], i, 0);
	}
	else
	{
		newNodeLegth = getRank(tree.bitStrings[minPos], length(tree.bitStrings[minPos]) - 1);
		resize(tree.bitStrings[minPos], newNodeLegth);
		for(unsigned i = 0; i < newNodeLegth; ++i)
			setBit(tree.bitStrings[minPos], i, 1);
	}

}

template <typename TText, typename TWaveletTreeSpec, typename TChar>
inline unsigned getNodePosition(WaveletTree<TText, TWaveletTreeSpec> & tree, TChar character)
{
	typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
	typename Iterator<TSplitValues>::Type iter(tree.splitValues, 0);
	return getNodePosition(iter, character);
}

//template <typename TText, typename TWaveletTreeSpec, typename TFreqTable>
template <typename TText, typename TWaveletTreeSpec, typename TFreqTable>//, typename TPrefixSumTable>
inline void createWaveletTree(WaveletTree<TText, TWaveletTreeSpec> & tree,
							  TText const & bwt,
                              TFreqTable const & freq)//,
                              //TPrefixSumTable & prefixSumTable)
{

    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;
    typedef typename Value<TWaveletTreeStructureEntry, 2>::Type TPointer;
    typedef typename Size<typename Value<TText>::Type>::Type TSize;


    //generate the tree structure
    TSize sigmaSize = length(freq);
    TPointer numberOfTreeNodes = sigmaSize - 1;
    resize(tree.splitValues, numberOfTreeNodes);
    computeTreeEntries(freq,
    		tree.splitValues);
    numberOfTreeNodes = getNumNodes(tree);
    resize(tree.splitValues,numberOfTreeNodes);
    resize(tree.bitStrings, numberOfTreeNodes);

    fillWaveletTree(tree, bwt);
}

//template <typename TBWT, typename TWaveletTreeSpec, typename TFreqTable, typename TPrefixSumTable, typename TDollarPosition>
template <typename TBWT, typename TWaveletTreeSpec, typename TFreqTable, typename TDollarPosition>
inline void createWaveletTree(WaveletTree<TBWT, TWaveletTreeSpec> & tree,
                              TBWT const & bwt,
                              TFreqTable const & freq,
                             // TPrefixSumTable & prefixSumTable,
                              TDollarPosition dollarPosition)
{
	typedef typename Fibre<WaveletTree<TBWT, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
	typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
    typedef typename Size<typename Value<TBWT>::Type>::Type TSize;

    //generate the tree structure
    TSize sigmaSize = length(freq);
    unsigned numberOfTreeNodes = sigmaSize - 1;

    resize(tree.splitValues, numberOfTreeNodes);//, Pair<TValue, TValue>(0, 0));
    TValue numChildNodes = 0;

    typedef typename Fibre<WaveletTree<TBWT, TWaveletTreeSpec>, FibreSplitValues>::Type TSplitValues;
    typename Iterator<TSplitValues>::Type iter(tree.splitValues, 0);
    computeTreeEntries(freq,
    				   tree.splitValues,
                       numChildNodes);

    numberOfTreeNodes = numChildNodes;

    resize(tree.bitStrings, numberOfTreeNodes);

    fillWaveletTree(tree, bwt, dollarPosition);
}



template <typename TBWT, typename TWaveletTreeSpec>
inline void createWaveletTree(WaveletTree<TBWT, TWaveletTreeSpec> & tree,
                              TBWT const & bwt)
{

    typedef String<typename Size<TBWT>::Type> TFreqTable;
    typedef typename Value<TFreqTable>::Type TFreqValue;
    typedef typename Value<TBWT>::Type TChar;

    //TFreqTable freqTable;
    String<Pair<TChar,TFreqValue> > newFreq;
    getNumChars(bwt, newFreq);

//    String<Pair<TChar,TFreqValue> > newFreq;
//    	for(TFreqValue i = 0; i < length(freqTable); ++i)
//    		if(freqTable[i])
//    			appendValue(newFreq, Pair<TChar, TFreqValue>((TChar)i, freqTable[i]));

    createWaveletTree(tree, bwt, newFreq);
}

template <typename TText>
inline bool openDollarInformation(
    WaveletTree<TText, DollarSubstituted> & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Size<TText>::Type TSize;

    typedef typename Fibre<WaveletTree<TText, DollarSubstituted>, FibreSplitValues>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeNodes>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;

    String<Pair<TChar, TSize> > dollarValues;

    name = fileName;    append(name, ".dollar");
    if (!open(dollarValues, toCString(name), openMode))
    {
        return false;
    }
    tree.dollarSub = dollarValues[0].i1;
    tree.dollarPosition = dollarValues[0].i2;
    return true;
}

template <typename TText>
inline bool openDollarInformation(
    WaveletTree<TText, MultiDollarSubstituted> & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
	typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreSplitValues>::Value TValue;
    typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>::Type TDollarString;

    String<Pair<TValue, TDollarString> > dollarValues;

    name = fileName;    append(name, ".dollar");
    if (!open(dollarValues, toCString(name), openMode))
    {
        return false;
    }
    tree.dollarSub = dollarValues[0].i1;
    tree.dollarPosition = dollarValues[0].i2;
    return true;
}

template <typename TText>
inline bool open(
    WaveletTree<TText, DollarSubstituted> & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".tree");  open(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
    name = fileName;    append(name, ".split"); open(getFibre(tree, FibreSplitValues()), toCString(name), openMode);
    openDollarInformation(tree, fileName, openMode);
    return true;
}

template <typename TText, typename TSpec>
inline bool open(
    WaveletTree<TText, TSpec> & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".tree");  open(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
    name = fileName;    append(name, ".split"); open(getFibre(tree, FibreSplitValues()), toCString(name), openMode);
    return true;
}

template <typename TText>
inline bool open(
    WaveletTree<TText, DollarSubstituted> & tree,
    const char * fileName)
{
    return open(tree, fileName, DefaultOpenMode<WaveletTree<TText, DollarSubstituted> >::VALUE);
}

template <typename TText, typename TSpec>
inline bool saveDollarInformation(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName,
    int openMode)
{
	WaveletTree<TText, TSpec> dummyTree = tree;
	return true;
}

template <typename TText>
inline bool saveDollarInformation(
    WaveletTree<TText, MultiDollarSubstituted> const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
    typedef typename Fibre<WaveletTree<TText, MultiDollarSubstituted>, FibreDollarPositions>::Type TDollarString;

    String<Pair<TValue, TDollarString> > dollarValues;
    append(dollarValues, Pair<TValue, TDollarString>(tree.dollarSub, tree.dollarPosition));

    name = fileName;    append(name, ".dollar"); save(dollarValues, toCString(name), openMode);
    return true;
}

template <typename TText>
inline bool saveDollarInformation(
    WaveletTree<TText, DollarSubstituted> const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Fibre<WaveletTree<TText, DollarSubstituted>, FibreSplitValues>::Type TSplitValues;
    typedef typename Fibre<TSplitValues, FibreTreeNodes>::Value TValue;
    typedef typename Size<TText>::Type                            TSize;

    String<Pair<TValue, TSize> > dollarValues;
    append(dollarValues, Pair<TValue, TSize>(tree.dollarSub, tree.dollarPosition));

    name = fileName;    append(name, ".dollar"); save(dollarValues, toCString(name), openMode);
    return true;
}

template <typename TText, typename TSpec>
inline bool save(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".tree");  save(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
    name = fileName;    append(name, ".split"); save(getFibre(tree, FibreSplitValues()), toCString(name), openMode);
    saveDollarInformation(tree, fileName, openMode);
    return true;
}

template <typename TText>
inline bool save(
    WaveletTree<TText, DollarSubstituted> const & tree,
    const char * fileName)
{
    return save(tree, fileName, DefaultOpenMode<WaveletTree<TText, DollarSubstituted> >::VALUE);
}

}
#endif  // #ifndef SANDBOX_MY_SANDBOX_APPS_FMINDEX_WAVELETTREE_H_
