// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
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

struct FibreBitStrings_;
struct FibreTreeStructure_;
struct FibreDollarPosition_;

template <typename TText, typename TSpec>
class WaveletTree;

// ==========================================================================
// Tags
// ==========================================================================

typedef Tag<FibreTreeStructure_> const FibreTreeStructure;
typedef Tag<FibreBitStrings_> const FibreBitStrings;
typedef Tag<FibreDollarPosition_> const FibreDollarPosition;

// ==========================================================================
// Metafunctions
// ==========================================================================
/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a WaveletTree.
..cat:WaveletTree

..tag.FibreBitStrings:The string set containing a bit string for each node.

..tag.FibreTreeStructure:The wavelet tree structure of the wavelet tree.

..tag.FibreDollarPositions:The bit string encoding the position of the dollar sign.
...remarks:This fibre is only available if the wavelet tree is used as the
occurrence table data structure of a FM index.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>
{
    typedef StringSet<RankSupportBitString<void> > Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec> const, FibreBitStrings>
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreTreeStructure>
{
    typedef typename Value<WaveletTree<TText, TSpec> >::Type    TChar;
    typedef typename MakeUnsigned<TChar>::Type                  TUChar;
    typedef RightArrayBinaryTree<TUChar, void>                  Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec> const, FibreTreeStructure>
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreTreeStructure>::Type const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>
{
    typedef Nothing Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, TSpec> const, FibreDollarPosition>
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>::Type const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>
{
    typedef typename Size<TText>::Type Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > const, FibreDollarPosition>
{
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type const Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > >, FibreDollarPosition>
{
    typedef RankSupportBitString<void> Type;
};

template <typename TText, typename TSpec>
struct Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const, FibreDollarPosition>
{
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > >, FibreDollarPosition>::Type const Type;
};

// ==========================================================================
template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> > :
    public Value<TText> {};

template <typename TText, typename TSpec>
struct Value<WaveletTree<TText, TSpec> const> :
    public Value<WaveletTree<TText, TSpec> > {};

template <typename TText, typename TSetSpec, typename TSpec>
struct Value<WaveletTree<StringSet<TText, TSetSpec>, TSpec> > :
    public Value<TText> {};

template <typename TText, typename TSetSpec, typename TSpec>
struct Value<WaveletTree<StringSet<TText, TSetSpec> const, TSpec> > :
    public Value<WaveletTree<StringSet<TText, TSetSpec>, TSpec> > {};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.WaveletTree:
..cat:Index
..summary:A wavelet tree is a tree like binary encoding of a text.
..signature:WaveletTree<TText, TSpec>
..param.TText:The value type of the text.
..param.TSpec:The wavelet tree specialisation.
...value:FmiDollarSubstituted
...value:void
...default:void.
..include:seqan/index.h
*/
template <typename TText, typename TSpec = void>
class WaveletTree
{
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreBitStrings>::Type    TBitStrings;
    typedef typename Fibre<WaveletTree<TText, TSpec>, FibreTreeStructure>::Type TWaveletTreeStructure;

public:
    TBitStrings bitStrings;
    TWaveletTreeStructure waveletTreeStructure;

    WaveletTree() {}

    explicit
    WaveletTree(TText const & text)
    {
        waveletTreeCreate(*this, text);
    }

    template <typename TFreqTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable)
    {
        waveletTreeCreate(*this, text, freqTable);
    }

    template <typename TFreqTable, typename TPrefixSumTable>
    WaveletTree(TText const & text, TFreqTable const & freqTable, TPrefixSumTable const & prefixSumTable)
    {
        waveletTreeCreate(*this, text, freqTable, prefixSumTable);
    }

    WaveletTree & operator=(WaveletTree const & other)
    {
        bitStrings = other.bitStrings;
        waveletTreeStructure = other.waveletTreeStructure;
        return *this;
    }

    bool operator==(WaveletTree const & b) const
    {
        typedef typename Size<TText>::Type TSize;
        
        if (length(bitStrings) == length(b.bitStrings))
        {
            for (TSize i = 0; i < length(bitStrings); ++i)
                if (!(bitStrings[i] == b.bitStrings[i]))
                    return false;
        }
        else
            return false;

        return waveletTreeStructure == b.waveletTreeStructure;
    }
};

template <typename TText, typename TSpec>
class WaveletTree<TText, FmiDollarSubstituted<TSpec> >:
    public WaveletTree<TText, void>
{
    typedef WaveletTree<TText, void> TBase;
    typedef typename Value<TText>::Type TChar;

public:
    typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type dollarPosition;
    TChar dollarSubstitute;

    WaveletTree() : dollarPosition(), dollarSubstitute() {}

    explicit
    WaveletTree(TText const & text) :
            TBase(text), dollarPosition(), dollarSubstitute()
    {}

    bool operator==(const WaveletTree & b) const
    {
        return static_cast<typename WaveletTree::TBase const &>(*this) == static_cast<typename WaveletTree::TBase const &>(b) &&
               dollarPosition == b.dollarPosition &&
               dollarSubstitute == b.dollarSubstitute;
    }

};

// ==========================================================================
//Functions
// ==========================================================================

///.Function.clear.param.type:Class.WaveletTree
template <typename TText, typename TSpec>
inline void clear(WaveletTree<TText, TSpec> & tree)
{
    clear(getFibre(tree, FibreBitStrings()));
    clear(getFibre(tree, FibreTreeStructure()));
}

// ==========================================================================
/**
.Function.dollarPosition
..summary:Returns whether a specified position is a dollar position.
..signature:dollarPosition(waveletTree, pos)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.pos:The position.
...type:Concept.UnsignedIntegerConcept
..include:seqan/index.h
*/

///.Function.dollarPosition.param.type:Class.WaveletTree
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

// ==========================================================================
///.Function.empty.param.type:Class.WaveletTree
template <typename TText, typename TSpec>
inline bool empty(WaveletTree<TText, TSpec> const & tree)
{
    return empty(getFibre(tree, FibreTreeStructure()));
}

// ==========================================================================
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

std::cout << getCharacter(waveletTree, 3) << std::endl; // T
std::cout << getCharacter(waveletTree, 'a', 4) << std::endl; // A
*/
template <typename TText, typename TWaveletTreeSpec, typename TPos>
inline typename Value<TText>::Type
_getCharacter(WaveletTree<TText, TWaveletTreeSpec> const & tree,
                 TPos pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreTreeStructure>::Type const    TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type                 TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                                       TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                                     TChar;

    unsigned treePos = 0;

    typename Iterator<TWaveletTreeStructure, TopDown<> >::Type iter(tree.waveletTreeStructure, treePos);
    TChar character = tree.waveletTreeStructure.minCharValue;
    do
    {
        TPos rank1 = getRank(tree.bitStrings[treePos], pos);
        if (getBit(tree.bitStrings[treePos], pos))
        {
            character = getCharacter(iter); 
            pos = rank1 - 1;  // -1 because strings start at 0
            if (!goRightChild(iter))
                break;
        }
        else
        {
            pos -= rank1;
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
getCharacter(WaveletTree<TText, TWaveletTreeSpec> const & tree,
             TPos pos)
{
    return _getCharacter(tree, pos);
}

// ==========================================================================
/**
.Function.getDollarPosition
..summary:Returns the dollar position.
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

// ==========================================================================
/**
.Function.WaveletTree#getFibre:
..summary:Returns a specific fibre of a container.
..signature:getFibre(container, fibreTag)
..class:Class.WaveletTree
..cat:Index
..param.container:The container holding the fibre.
...type:Class.WaveletTree
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.WaveletTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
..example.code:
Index< String<char> > index_esa("tobeornottobe");

String<char> & text = getFibre(indexEsa, EsaText());
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
inline typename Fibre<WaveletTree<TText, TSpec>, FibreTreeStructure>::Type &
getFibre(WaveletTree<TText, TSpec>&tree, FibreTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreTreeStructure>::Type const &
getFibre(WaveletTree<TText, TSpec> const & tree, const FibreTreeStructure)
{
    return tree.waveletTreeStructure;
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>::Type
getFibre(WaveletTree<TText, TSpec>&/*tag*/, FibreDollarPosition)
{
    return Nothing();
}

template <typename TText, typename TSpec>
inline typename Fibre<WaveletTree<TText, TSpec>, FibreDollarPosition>::Type
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

// ==========================================================================
// This functions computes the number of occurrences of a specified character
// up to a specified position.
template <
    typename TText,
    typename TWaveletTreeSpec,
    typename TCharIn,
    typename TPos>
inline unsigned _countOccurrences(WaveletTree<TText, TWaveletTreeSpec> const & tree,
                                   TCharIn const character,
                                   TPos const pos)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreTreeStructure>::Type TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type TChar;

    TPos sum = pos;
    TPos treePos = 0;

    typename Iterator<TWaveletTreeStructure const, TopDown<> >::Type it(tree.waveletTreeStructure, treePos);
    TChar charInTree = tree.waveletTreeStructure.minCharValue;
    do
    {
        TPos addValue = getRank(tree.bitStrings[treePos], sum);
        if (ordGreater(getCharacter(it), character))
        {
            if (addValue > sum) return 0;

            sum -= addValue;
            if (!goLeftChild(it))
                break;
        }
        else
        {
            if (addValue == 0) return 0;

            charInTree = getCharacter(it);
            sum = addValue - 1;
            if (!goRightChild(it))
                break;
        }
        treePos = getPosition(it);
    }
    while (true);

    if (ordEqual(charInTree, character))
        return sum + 1;

    return 0;
}

/**
.Function.WaveletTree#countOccurrences
..summary:Returns the number of occurrences of a specified character from the start
to a specified position.
..signature:countOccurrences(waveletTree, character, pos)
..param.waveletTree:The wavelet tree.
...type:Class.WaveletTree
..param.character:The character.
..param.pos:The position
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";
WaveletTree<String<Dna5> > waveletTree(genome);

std::cout << countOccurrences(waveletTree, 'a', 3) << std::endl; // 1
std::cout << countOccurrences(waveletTree, 'a', 4) << std::endl; // 2
*/
template <typename TText, typename TChar, typename TPos, typename TWaveletTreeSpec>
inline unsigned countOccurrences(WaveletTree<TText, TWaveletTreeSpec> const & tree,
                               TChar const character,
                               TPos const pos)
{
    return _countOccurrences(tree, character, pos);
}

template <typename TText, typename TSpec, typename TChar, typename TPos>
inline unsigned countOccurrences(WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > const & tree,
                               TChar const character,
                               TPos const pos)
{
    unsigned occ = _countOccurrences(tree, character, pos);
    if (ordEqual(getDollarSubstitute(tree), character) && pos >= tree.dollarPosition)
         --occ;
    return occ;
}

template <typename TText, typename TSpec, typename TChar, typename TPos>
inline unsigned countOccurrences(WaveletTree<TText, FmiDollarSubstituted<MultiDollar<TSpec> > > const & tree,
                               TChar const character,
                               TPos const pos)
{
    unsigned occ = _countOccurrences(tree, character, pos);
    if (ordEqual(getDollarSubstitute(tree), character))
        return occ - getRank(getFibre(tree, FibreDollarPosition()), pos);

    return occ;
}

// ==========================================================================
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
getDollarSubstitute(WaveletTree<TText, FmiDollarSubstituted<TWaveletTreeSpec> > const & tree)
{
    return tree.dollarSubstitute;
}

// ==========================================================================
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
    tree.dollarPosition = position;
}

// ==========================================================================
// This function is used to fill the bit strings of the wavelet tree.
template <typename TText, typename TWaveletTreeSpec, typename TText2>
inline void _fillWaveletTree(WaveletTree<TText, TWaveletTreeSpec> & tree,
                             TText2 const & text)
{
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreBitStrings>::Type     TFibreRankSupportBitStrings;
    typedef typename Value<TFibreRankSupportBitStrings>::Type                               TFibreRankSupportBitString;
    typedef typename Fibre<TFibreRankSupportBitString, FibreBits>::Type                     TFibreBitString;
    typedef typename Size<TFibreBitString>::Type                                            TSize;
    typedef typename Fibre<WaveletTree<TText, TWaveletTreeSpec>, FibreTreeStructure>::Type  TWaveletTreeStructure;
    //typedef typename Value<TText>::Type                                                     TChar;
    //typedef typename MakeUnsigned<TChar>::Type                                              TUChar;

    resize(tree.bitStrings, _length(tree.waveletTreeStructure));

    for (TSize i = 0; i < length(text); ++i)
    {
        typename Iterator<TWaveletTreeStructure, TopDown<> >::Type it(tree.waveletTreeStructure, 0);
        bool bit;

        do
        {
            if (ordGreater(getCharacter(it), getValue(text, i)))
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
        _updateRanks(bitStrings[i]);
    }
}

// ==========================================================================
/**
.Function.waveletTreeStructureCreate
..summary:This functions creates the wavelet tree structure.
..signature:void waveletTreeStructureCreate(waveletTreeStructure, text)
..param.waveletTreeStructure:The wavelet tree structure.
...type:Class.RightArrayBinaryTree
..param.text:A text to be transfered into a wavelet tree.
...type:Class.String
..include:seqan/index.h
*/
template <typename TWaveletTreeStructure, typename TText>
inline void waveletTreeStructureCreate(TWaveletTreeStructure & waveletTreeStructure, TText const & text)
{
    computeTree(waveletTreeStructure, text);
}

/**
.Function.waveletTreeStructureCreate
..signature:void waveletTreeStructureCreate(lfTable)
..param.lfTable:The LF table of the text to be transfered into a wavelet tree structure.
...type:Class.LfTable
*/
// create structure from prefix sum array
template <typename TLfTable>
inline void waveletTreeStructureCreate(TLfTable & lfTable)
{
    computeTree(lfTable);
}

// ==========================================================================
/**
.Function.waveletTreeCreate
..summary:This functions creates the wavelet tree structure.
..signature:void waveletTreeStructureCreate(waveletTree, text)
..param.waveletTree:The wavelet tree.
...type:ClassWaveletTree.
..param.text:A text to be transfered into a wavelet tree.
...type:Class.String
..include:seqan/index.h
*/

template <typename TText, typename TWaveletTreeSpec, typename TText2>
inline void waveletTreeCreate(WaveletTree<TText, TWaveletTreeSpec> & waveletTree,
                              TText2 const & text)
{
    waveletTreeStructureCreate(getFibre(waveletTree, FibreTreeStructure()), text);
    _fillWaveletTree(waveletTree, text);
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
    _fillWaveletTree(getFibre(lfTable, FibreOccTable()), text);
}

// ==========================================================================
template <typename TText, typename TSpec>
inline bool openDollarInformation(
    WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;

    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type TDollarString;
    typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > >::Type TChar;

    String<Pair<__int64, TDollarString, Pack> > dollarValues;

    name = fileName;    append(name, ".dr");
    if (!open(dollarValues, toCString(name), openMode) || empty(dollarValues))
    {
        return false;
    }
    tree.dollarSubstitute = static_cast<TChar>(dollarValues[0].i1);
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
    //typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type TDollarString;

    String<TChar> dollarSub;
    
    name = fileName;    append(name, ".drs"); if (!open(dollarSub, toCString(name), openMode)) return false;
    name = fileName;    append(name, ".drp"); if (!open(tree.dollarPosition, toCString(name), openMode)) return false;
    
    if (empty(dollarSub))
        return false;

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
    name = fileName;    append(name, ".wtc"); if (!open(getFibre(tree, FibreBitStrings()), toCString(name), openMode)) return false;
    name = fileName;    append(name, ".wts"); if (!open(getFibre(tree, FibreTreeStructure()), toCString(name), openMode)) return false;
    if (!openDollarInformation(tree, fileName, openMode)) return false;
    return true;
}

/**
.Function.open
..param.object:
...type:Class.WaveletTree
*/
template <typename TText, typename TSpec>
inline bool open(
    WaveletTree<TText, TSpec> & tree,
    const char * fileName)
{
    return open(tree, fileName, DefaultOpenMode<WaveletTree<TText, TSpec> >::VALUE);
}

// ==========================================================================
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

    //typedef typename Value<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > > >::Type TChar;
    typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<SingleDollar<TSpec> > >, FibreDollarPosition>::Type TDollarString;

    String<Pair<__int64, TDollarString, Pack> > dollarValues;
    appendValue(dollarValues, Pair<__int64, TDollarString, Pack>(ordValue(tree.dollarSubstitute), tree.dollarPosition));

    name = fileName;    append(name, ".dr"); if (!save(dollarValues, toCString(name), openMode)) return false;
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
    //typedef typename Fibre<WaveletTree<TText, FmiDollarSubstituted<TSpec> >, FibreDollarPosition>::Type TDollarString;

    String<TChar> dollarSub;
    appendValue(dollarSub, tree.dollarSubstitute);

    name = fileName;    append(name, ".drs"); if(!save(dollarSub, toCString(name), openMode)) return false;
    name = fileName;    append(name, ".drp"); if(!save(tree.dollarPosition, toCString(name), openMode)) return false;
    return true;
}

template <typename TText, typename TSpec>
inline bool save(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName,
    int openMode)
{
    String<char> name;
    name = fileName;    append(name, ".wtc");   if (!save(getFibre(tree, FibreBitStrings()), toCString(name), openMode)) return false;
    name = fileName;    append(name, ".wts");   if (!save(getFibre(tree, FibreTreeStructure()), toCString(name), openMode)) return false;
    if (!saveDollarInformation(tree, fileName, openMode)) return false;
    
    return true;
}
/**
.Function.save
..param.object:
...type:Class.WaveletTree
*/
template <typename TText, typename TSpec>
inline bool save(
    WaveletTree<TText, TSpec> const & tree,
    const char * fileName)
{
    return save(tree, fileName, DefaultOpenMode<WaveletTree<TText, TSpec> >::VALUE);
}
}
#endif  // INDEX_FM_WAVELET_TREE_H_
