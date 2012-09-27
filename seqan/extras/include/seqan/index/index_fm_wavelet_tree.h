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

#ifndef INDEX_FM_WAVELET_TREE_H_
#define INDEX_FM_WAVELET_TREE_H_

namespace seqan {

// ==========================================================================
//Forwards
// ==========================================================================
template <typename TSpec>
struct SingleDollar;

template <typename TSpec>
struct MultiDollar;

template <typename TSpec = void>
struct FmiDollarSubstituted;

/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a WaveletTree.
..cat:WaveletTree

..tag.FibreBitStrings:The string set containing a bit string for each node.

..tag.FibreWaveletTreeStructure:The wavelet tree structure of the wavelet tree.

..tag.FibreDollarPositions:The bit string encoding the position of the dollar sign.
...remarks:This fibre is only available if the wavelet tree is used as the
occurrence table data structure of a FM index.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres
struct FibreBitStrings_;
struct FibreWaveletTreeStructure_;
struct FibreDollarPosition_;
struct FibreTreeStructureEncoding_;

template <typename TText, typename TSpec>
class WaveletTree;

// ==========================================================================
// Tags
// ==========================================================================

typedef Tag<FibreWaveletTreeStructure_> const FibreWaveletTreeStructure;
typedef Tag<FibreBitStrings_> const FibreBitStrings;
typedef Tag<FibreDollarPosition_> const FibreDollarPosition;
typedef Tag<FibreTreeStructureEncoding_> const FibreTreeStructureEncoding;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>
{
    typedef typename Value<TText>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;
    typedef RightArrayBinaryTree<TUChar, void> Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>
{
    typedef Nothing Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec> const, FibreDollarPosition>
{
    typedef Nothing const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>
{
    typedef typename Size<TText>::Type Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > const, FibreDollarPosition>
{
    typedef typename Size<TText>::Type const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > >, FibreDollarPosition>
{
    typedef RankSupportBitString<void> Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const, FibreDollarPosition>
{
    typedef RankSupportBitString<void> const Type;
};

// ==========================================================================
template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> >
{
    typedef typename Value<TText>::Type Type;
};

template <typename TText, typename TSetSpec, typename TSpec>
struct Value<WaveletTree<StringSet<TText, TSetSpec>, TSpec> >
{
    typedef typename Value<TText>::Type Type;
};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> const>
{
    typedef typename Value<WaveletTree<TText, TSpec> >::Type const Type;
};

template <typename TText, typename TSetSpec, typename TSpec>
struct Value<WaveletTree<StringSet<TText, TSetSpec> const, TSpec> >
{
    typedef typename Value<WaveletTree<StringSet<TText, TSetSpec>, TSpec> >::Type const Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.WaveletTree:
..cat:Graph
..summary:A wavelet tree is a tree like binary encoding of a text.
..signature:WaveletTree<TText, TSpec>
..param.TText:The value type of the text.
..param.TSpec:The wavelet tree specialisation.
...tag:FmiDollarSubstituted
...default:void.
..include:seqan/index.h
*/
template <typename TText, typename TSpec = void>
class WaveletTree
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type TBitStrings;
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;

public:
    TBitStrings bitStrings;
    TWaveletTreeStructure waveletTreeStructure;

    WaveletTree() :
        bitStrings(),
        waveletTreeStructure()
    {}

    WaveletTree(TText const & text) :
        bitStrings(),
        waveletTreeStructure()
    {
        waveletTreeCreate(*this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable) :
        bitStrings(),
        waveletTreeStructure()
    {
        waveletTreeCreate(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable) :
        bitStrings(),
        waveletTreeStructure()
    {
        waveletTreeCreate(*this, text, freqTable, prefixSumTable);
    }

    inline WaveletTree & operator=(WaveletTree const & other)
    {
        bitStrings = other.bitStrings;
        waveletTreeStructure = other.waveletTreeStructure;
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
               waveletTreeStructure == b.waveletTreeStructure;
    }

};

template <typename TText, typename TSpec>
class WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >:
    public WaveletTree<TText, void>
{
    typedef WaveletTree<TText, void> TBase;
    typedef typename Value<TText>::Type TChar;

public:
    typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type dollarPosition;
    TChar dollarSubstitute;

    WaveletTree() :
        TBase(),
        dollarPosition(),
        dollarSubstitute()
    {}

    WaveletTree(TText const & text) :
        TBase(text),
        dollarPosition(),
        dollarSubstitute()
    {}

    inline bool operator==(const WaveletTree & b) const
    {
        return static_cast<typename WaveletTree::TBase>(*this) == static_cast<typename WaveletTree::TBase>(b) &&
               dollarPosition == b.dollarPosition &&
               dollarSubstitute == b.dollarSubstitute;
    }

};

template <typename TText, typename TSpec>
class WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > >:
    public WaveletTree<TText, void>
{
    typedef WaveletTree<TText, void> TBase;
    typedef typename Value<TText>::Type TChar;

public:
    typename Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > >, FibreDollarPosition>::Type dollarPosition;
    TChar dollarSubstitute;

    WaveletTree() :
        TBase(),
        dollarPosition(),
        dollarSubstitute()
    {}

    WaveletTree(TText const & text) :
        TBase(text),
        dollarPosition(),
        dollarSubstitute()
    {}

    inline bool operator==(const WaveletTree & b) const
    {
        return static_cast<typename WaveletTree::TBase>(*this) == static_cast<typename WaveletTree::TBase>(b) &&
               dollarPosition == b.dollarPosition &&
               dollarSubstitute == b.dollarSubstitute;
    }

};

// ==========================================================================
//Functions
// ==========================================================================

/**
.Function.clear
..param.object:
...type:Class.WaveletTree
*/
template <typename TText, typename TSpec>
inline void clear(WaveletTree<TText, TSpec> & tree)
{
    for (unsigned i = 0; i < length(tree.bitStrings); ++i)
    {
        clear(tree.bitStrings[i]);
    }
    resize(tree.bitStrings, 0);
    clear(tree.waveletTreeStructure);
}

template <typename TText, typename TSpec, typename TPos>
inline bool dollarPosition(WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > & tree, TPos pos)
{
    return tree.dollarPosition == pos;
}

template <typename TText, typename TSpec, typename TPos>
inline bool dollarPosition(WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > const & tree, TPos pos)
{
    return tree.dollarPosition == pos;
}

template <typename TText, typename TSpec, typename TPos>
inline bool dollarPosition(WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const & tree, TPos pos)
{
    return getBit(tree.dollarPosition, pos);
}

/**
.Function.empty
..param.object:
...type:Class.WaveletTree
*/
template <typename TText, typename TSpec>
inline bool empty(WaveletTree<TText, TSpec> & tree)
{
    return empty(getFibre(tree, FibreWaveletTreeStructure()));
}

template <typename TText, typename TSpec>
inline bool empty(WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > & tree)
{
    return empty(getFibre(tree, FibreWaveletTreeStructure())) &&
           empty(getFibre(tree, FibreDollarPosition()));
}

/**
.Function.getCharacter
..summary:Returns the character of a specified position.
..signature:getCharacter(waveletTree, pos)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.pos:The position
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTree<String<Dna5> > waveletTree(genome);

std::cerr << getCharacter(waveletTree, 3) << std::endl; // T
std::cerr << getCharacter(waveletTree, 'a', 4) << std::endl; // A
*/
template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline typename Value<TText>::Type
getCharacterImpl(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                 const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type const TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type                             TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                                           TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                                         TChar;

    unsigned long sum = pos + 1;
    unsigned treePos = 0;

    typename Iterator<TWaveletTreeStructure, TopDown<> >::Type iter(tree.waveletTreeStructure, treePos);
    bool direction;
    TChar character = tree.waveletTreeStructure.minCharValue;
    do
    {
        direction = getBit(tree.bitStrings[treePos], sum - 1);
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        if (direction)
        {
            character = getCharacter(iter); // + 1;
            sum = addValue;
            if (!goRightChild(iter))
                break;
        }
        else
        {
            sum -= addValue;
            if (!goLeftChild(iter))
                break;
        }
        treePos = getPosition(iter);
    }
    while (true);

    return character;
}

template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline typename Value<TText>::Type
getCharacter(const WaveletTree<TText, TWaveletTreeSpec> & tree,
             const TPos pos)
{
    return getCharacterImpl(tree, pos);
}

template <typename TText, typename TSpec, typename TPos>
inline typename Value<TText>::Type
getCharacter(WaveletTree<TText, FmiDollarSubstituted<TSpec> > const & tree,
             TPos const pos)
{
    #if SEQAN_ENABLE_DEBUG
    if (dollarPosition(tree, pos))
        std::cout << "Note: the character of the requested position (the dollar position) is the dollar substitute." << std::endl;
    #endif

    return getCharacterImpl(tree, pos);
}

/**
.Function.getDollarPosition
..summary:Sets the dollar position..
..signature:getDollarPosition(waveletTree)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..include:seqan/index.h
*/
template <typename TText, typename TSpec>
inline unsigned long getDollarPosition(WaveletTree<TText, FmiDollarSubstituted<TSpec> > const & tree)
{
    return tree.dollarPosition;
}

template <typename TText, typename TSpec>
inline String<unsigned long> getDollarPosition(WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const & tree)
{
    String<unsigned long> dollarPositions;

    for (unsigned i = 0; i < length(tree.dollarPosition); ++i)
        if (getBit(tree.dollarPosition))
            appendValue(dollarPositions, i);
    return dollarPositions;
}

// TODO (singer): Decide whether we need this function.
/*
.Function.getAlphabet
..summary:Returns a string of characters present in the wavelet tree.
..signature:getAlphabet(waveletTree)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ATTAGGTGA";
WaveletTree<String<Dna5> > waveletTree(genome);

std::cerr << getAlphabet(waveletTree) << std::endl; // AGT
*/
// template <typename TText, typename TWaveletTreeSpec>
// inline TText getAlphabet(WaveletTree<TText, TWaveletTreeSpec> & tree)
// {
//     getAlphabet(tree.splitValues);
// }


/**
.Function.getFibre:
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.WaveletTree Fibres
*/
template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type &
getFibre(WaveletTree<TText, TSpec>&tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type const &
getFibre(const WaveletTree<TText, TSpec>&tree, const FibreBitStrings)
{
    return tree.bitStrings;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type &
getFibre(WaveletTree<TText, TSpec>&tree, FibreWaveletTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreWaveletTreeStructure>::Type const &
getFibre(WaveletTree<TText, TSpec> const & tree, const FibreWaveletTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>::Type &
getFibre(WaveletTree<TText, TSpec>&/*tag*/, FibreDollarPosition)
{
    return Nothing();
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>::Type const &
getFibre(WaveletTree<TText, TSpec> const & /*tag*/, const FibreDollarPosition)
{
    return Nothing();
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type &
getFibre(WaveletTree<TText, FmiDollarSubstituted<TSpec> >&tree, FibreDollarPosition)
{
    return tree.dollarPosition;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type const &
getFibre(WaveletTree<TText, FmiDollarSubstituted<TSpec> > const & tree, const FibreDollarPosition)
{
    return tree.dollarPosition;
}




// This functions computes the number of occurrences of a specified character
// up to a specified position.
template <
    typename TText,
    typename TWaveletTreeSpec,
    typename TCharIn,
    typename TPos>
inline unsigned getOccurrencesImpl(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                                   const TCharIn character,
                                   const TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;

    TPos sum = pos + 1;
    TPos treePos = 0;

    typename Iterator<TWaveletTreeStructure const, TopDown<> >::Type it(tree.waveletTreeStructure, treePos);
    TChar charInTree = tree.waveletTreeStructure.minCharValue;
    do
    {
        TPos addValue = getRank(tree.bitStrings[treePos], sum - 1);
        if (character < getCharacter(it)) //getFibre(tree, FibreWaveletTreeStructure()).treeVertieces[treePos].i1)
        {
            sum -= addValue;
            if (!goLeftChild(it))
                break;
        }
        else
        {
            charInTree = getCharacter(it);
            sum = addValue;
            if (!goRightChild(it))
                break;
        }
        treePos = getPosition(it);
    }
    while (sum);

    //static_cast<Nothing>(character);

    //std::cerr << "character: " << character << " " << ordValue(character) << " " << charInTree << " " << (int)charInTree << std::endl;

    if (character == charInTree)
        return sum;

    return 0;
}

/**
.Function.getOccurrences
..summary:Returns the number of occurrences of a specified character from the start
to a specified position.
..signature:getOccurrences(waveletTree, character, pos)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.character:The character.
..param.pos:The position
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTree<String<Dna5> > waveletTree(genome);

std::cerr << getOccurrences(waveletTree, 'a', 3) << std::endl; // 1
std::cerr << getOccurrences(waveletTree, 'a', 4) << std::endl; // 2
*/
template <typename TText, typename TChar, typename TPos, typename TWaveletTreeSpec>
inline unsigned getOccurrences(const WaveletTree<TText, TWaveletTreeSpec> & tree,
                               const TChar character,
                               const TPos pos)
{
    typedef typename MakeUnsigned<TChar>::Type TUChar;

    return getOccurrencesImpl(tree, static_cast<TUChar>(character), pos);
}

template <typename TText, typename TSpec, typename TChar, typename TPos>
inline unsigned getOccurrences(const WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > & tree,
                               const TChar character,
                               const TPos pos)
{
    typedef typename MakeUnsigned<TChar>::Type TUChar;

    unsigned occ = getOccurrencesImpl(tree, static_cast<TUChar>(character), pos);
    if (character == getDollarSubstitute(tree) && pos >= tree.dollarPosition)
    {
        return occ - 1;
    }
    return occ;
}

template <typename TText, typename TSpec, typename TChar, typename TPos>
inline unsigned getOccurrences(const WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > & tree,
                               const TChar character,
                               const TPos pos)
{
    typedef typename MakeUnsigned<TChar>::Type TUChar;



    unsigned occ = getOccurrencesImpl(tree, static_cast<TUChar>(character), pos);
    if (character == getDollarSubstitute(tree))
        return occ - getRank(getFibre(tree, FibreDollarPosition()), pos);

    return occ;
}

/**
.Function.getDollarSubstitute
..summary:Returns the character used to substitute the dollar sign.
..signature:getDollarSubstitute(waveletTree)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGGTT";
WaveletTree<String<Dna5>, FmiDollarSubstituted> waveletTree(genome);

std::cerr << getDollarSubstitute(waveletTree) << std::endl; // A
*/
template <typename TText, typename TWaveletTreeSpec>
inline typename Value<TText>::Type
getDollarSubstitute(WaveletTree<TText, TWaveletTreeSpec> const & /*tag*/)
{
    return Nothing();
}

template <typename TText, typename TWaveletTreeSpec>
inline typename Value<TText>::Type
getDollarSubstitute(const WaveletTree<TText, FmiDollarSubstituted<TWaveletTreeSpec> > & tree)
{
    return tree.dollarSubstitute;
}

// This function returns the distance to the root
// inline unsigned getTreeLevel_(unsigned treePosition)
// {
//     return floor(log(treePosition + 1) / log(2));
// }


/**
.Function.numVertices
..summary:Returns the number of vertices in a wavelet tree.
..signature:numVertices(waveletTree)
..param.waveletTree:The wavelet tree
...type:Class.WaveletTree
..returns: unsigned
*/
template <typename TText, typename TSpec>
inline unsigned numVertieces(WaveletTree<TText, TSpec> & tree)
{
    return length(tree.waveletTreeStructure.treeVertieces);
}

/**
.Function.setDollarSubstitute
..summary:Sets the character used to substitute the dollar sign.
..signature:setDollarSubstitute(waveletTree, character)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.character:The dollar substitute.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGGTT";
WaveletTree<String<Dna5>, FmiDollarSubstituted> waveletTree(genome);

std::cerr << getDollarSubstitute(waveletTree) << std::endl; // A

setDollarSubstitute(waveletTree, 'G');

std::cerr << getDollarSubstitute(waveletTree) << std::endl; // G
*/
template <typename TText, typename TWaveletTreeSpec, typename TChar>
inline void setDollarSubstitute(WaveletTree<TText, FmiDollarSubstituted<TWaveletTreeSpec> > & tree,
                                TChar dollarSubstitute)
{
    tree.dollarSubstitute = dollarSubstitute;
}

/**
.Function.setDollarPosition
..summary:Sets the dollar position..
..signature:setDollarPosition(waveletTree, pos)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.pos:The dollar position.
..include:seqan/index.h
*/
template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline void setDollarPosition(WaveletTree<TText, FmiDollarSubstituted<TWaveletTreeSpec> > & tree,
                              TPos const & position)
{
    //static_cast<Nothing>(position);
    //static_cast<Nothing>(tree.dollarPosition);

    tree.dollarPosition = position;
}

// template <typename TText, typename TPos>
// inline void setDollarPosition(WaveletTree<StringSet<TText>, FmiDollarSubstituted<> > & tree,
//                               TPos position)
// {
//     setBit(tree.dollarPosition, position, true);
// }

// This function is used to fill the bit strings of the wavelet tree.
template <typename TText, typename TWaveletTreeSpec, typename TText2>
inline void fillWaveletTree_(WaveletTree<TText, TWaveletTreeSpec> & tree,
                             TText2 const & text)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type             TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type                                       TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreBits>::Type                        TFibreBitString;
    typedef typename Size<TFibreBitString>::Type                                                    TSize;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type   TWaveletTreeStructure;

    typedef typename Value<TText>::Type TChar;
    typedef typename MakeUnsigned<TChar>::Type TUChar;

    resize(tree.bitStrings, length(tree.waveletTreeStructure));

    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TWaveletTreeStructure, TopDown<> >::Type it(tree.waveletTreeStructure, 0);
        bool bit;

        do
        {
            if ((ordValue(getValue(text, i))) < getCharacter(it))
            {
                bit = 0;
                appendValue(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goLeftChild(it))
                    break;
            }
            else
            {
                bit = 1;
                appendValue(getFibre(tree, FibreBitStrings())[getPosition(it)], bit);
                if (!goRightChild(it))
                    break;
            }
        }
        while (true);
    }

    TFibreRankSupportBitStrings & bitStrings = getFibre(tree, FibreBitStrings());
    for (TSize i = 0; i < length(bitStrings); ++i)
    {
        TFibreRankSupportBitString & temp = bitStrings[i];
        updateRanks_(temp);
    }
}

template <typename TWaveletTreeStructure, typename TText>
inline void waveletTreeStructureCreate(TWaveletTreeStructure & waveletTreeStructure, TText const & text)
{
    computeRightArrayBinaryTree(waveletTreeStructure, text);
}

template <typename TLfTable>
inline void waveletTreeStructureCreate(TLfTable & lfTable)
{
    computeRightArrayBinaryTree(lfTable);
}

/*
.Function.getVertexPosition
..summary:Returns the position of the vertex in the wavelet tree structure fibre given a specific character.
..signature:getVertexPosition(waveletTree, character)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.character.The character
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGGTT";
WaveletTree<String<Dna5>, FmiDollarSubstituted> waveletTree(genome);

std::cerr << getVertexPosition(waveletTree, 'G') << std::endl; // 0
*/
// template <typename TText, typename TWaveletTreeSpec, typename TChar>
// inline unsigned getVertexPosition(WaveletTree<TText, TWaveletTreeSpec> & tree, TChar character)
// {
//     typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreWaveletTreeStructure>::Type TSplitValues;
//     typename Iterator<TSplitValues>::Type iter(tree.splitValues, 0);
//     return getVertexPosition(iter, character);
// }

// This function is called to create the wavelet tree.
template <typename TText, typename TWaveletTreeSpec, typename TText2>
inline void waveletTreeCreate(WaveletTree<TText, TWaveletTreeSpec> & waveletTree,
                              TText2 const & text)
{
    waveletTreeStructureCreate(getFibre(waveletTree, FibreWaveletTreeStructure()), text);
    fillWaveletTree_(waveletTree, text);
}

template <typename TText, typename TWaveletTreeSpec, typename TPrefixSumTable, typename TText2, typename TDollarSub, typename TDollarPos>
inline void waveletTreeCreate(LfTable<WaveletTree<TText, TWaveletTreeSpec>, TPrefixSumTable> & lfTable,
                              TText2 const & text,
                              TDollarSub const & dollarSub,
                              TDollarPos const & dollarPos)
{
    setDollarSubstitute(lfTable.occTable, dollarSub);
    setDollarPosition(lfTable.occTable, dollarPos);

    waveletTreeStructureCreate(lfTable);
    fillWaveletTree_(getFibre(lfTable, FibreOccTable()), text);
}

template <typename TText, typename TSpec>
inline bool openDollarInformation(
    WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type TDollarString;
    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > >::Type TChar;

    String<Pair<TChar, TDollarString> > dollarValues;

    name = fileName;    append(name, ".dr");
    if (!open(dollarValues, toCString(name), openMode))
    {
        return false;
    }
    tree.dollarSubstitute = dollarValues[0].i1;
    tree.dollarPosition = dollarValues[0].i2;
    return true;
}

template <typename TText, typename TSpec>
inline bool openDollarInformation(
    WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<TSpec> > >::Type TChar;
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type TDollarString;

    String<TChar> dollarSub;
    name = fileName;    append(name, ".drs"); open(dollarSub, toCString(name), openMode);
    name = fileName;    append(name, ".drp"); open(tree.dollarPosition, toCString(name), openMode);
    tree.dollarSubstitute = dollarSub[0];
    return true;
}


template <typename TText, typename TSpec>
inline bool openDollarInformation(
    WaveletTree<TText, TSpec> & /*tag*/,
    char const * /*tag*/,
    int const/*tag*/)
{ 
    return true;
}

template <typename TText, typename TSpec>
inline bool open(
    WaveletTree<TText, TSpec> & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".wtc"); open(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
    name = fileName;    append(name, ".wts"); open(getFibre(tree, FibreWaveletTreeStructure()), toCString(name), openMode);
    openDollarInformation(tree, fileName, openMode);
    return true;
}


template <typename TText, typename TSpec>
inline bool open(
    WaveletTree<TText, TSpec> & tree,
    const char * fileName)
{
    return open(tree, fileName, DefaultOpenMode<WaveletTree<TText, TSpec> >::VALUE);
}

template <typename TText, typename TSpec>
inline bool saveDollarInformation(
    WaveletTree<TText, TSpec> const & /*tag*/,
    const char * /*tag*/,
    int /*tag*/)
{
	return true;
}
template <typename TText, typename TSpec>
inline bool saveDollarInformation(
    WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > >::Type TChar;
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type TDollarString;

    String<Pair<TChar, TDollarString> > dollarValues;
    append(dollarValues, Pair<TChar, TDollarString>(tree.dollarSubstitute, tree.dollarPosition));

    name = fileName;    append(name, ".dr"); save(dollarValues, toCString(name), openMode);
    return true;
}

template <typename TText, typename TSpec>
inline bool saveDollarInformation(
    WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<TSpec> > >::Type TChar;
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type TDollarString;

    String<TChar> dollarSub;
    append(dollarSub, tree.dollarSubstitute);

    name = fileName;    append(name, ".drs"); save(dollarSub, toCString(name), openMode);
    name = fileName;    append(name, ".drp"); save(tree.dollarPosition, toCString(name), openMode);
    return true;
}

template <typename TText, typename TSpec>
inline bool save(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".wtc");   save(getFibre(tree, FibreBitStrings()), toCString(name), openMode);
    
    name = fileName;    append(name, ".wts");   save(getFibre(tree, FibreWaveletTreeStructure()), toCString(name), openMode);
    saveDollarInformation(tree, fileName, openMode);
    
    return true;
}

template <typename TText, typename TSpec>
inline bool save(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName)
{
    return save(tree, fileName, DefaultOpenMode<WaveletTree<TText, TSpec> >::VALUE);
}
}
#endif  // INDEX_FM_WAVELET_TREE_H_
