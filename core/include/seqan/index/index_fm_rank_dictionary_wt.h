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

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_RANKDICTIONARY_WT
#define INDEX_FM_RANKDICTIONARY_WT

namespace seqan {

// ==========================================================================
// Tags
// ==========================================================================

template <typename TSpec = void>
struct WaveletTree {};

struct FibreTreeStructure_;
typedef Tag<FibreTreeStructure_>    const FibreTreeStructure;

// ==========================================================================
// Metafunctions
// ==========================================================================
/**
.Tag.WaveletTree Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Spec.WaveletTree@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a @Spec.WaveletTree@.
..cat:Spec.WaveletTree

..tag.FibreRanks:The string set containing a bit string for each node.

..tag.FibreTreeStructure:The wavelet tree structure of the wavelet tree.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/
/*!
 * @defgroup WaveletTreeFibres WaveletTree Fibres
 * 
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link WaveletTree @endlink.
 * 
 * @section Remarks
 * 
 * These tags can be used to get @link Fibre Fibres @endlink of a @link WaveletTree @endlink.
 * 
 * @tag WaveletTreeFibres#FibreTreeStructure
 * 
 * @brief The wavelet tree structure of the wavelet tree.
 * 
 * @tag WaveletTreeFibres#FibreRanks
 * 
 * @brief A string set containing a rank support bit string for each node in the tree.
 *
 * @see Fibre
 * @see Index#getFibre
 * 
 */

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

///.Metafunction.Fibre.param.TContainer.type:Class.RankDictionary
///.Metafunction.Fibre.param.TSpec.type:Tag.WaveletTree Fibres
/*
.Metafunction.Fibre:
..summary:Type of a specific FMIndex member (fibre).
..signature:Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, TFibreSpec>::Type
..class:Spec.FMIndex
..cat:Index
..param.TValue:The character type of the @Class.String@ the wavelet tree represents.
..param.TFibreSpec:Tag to specify the fibre.
...type:Tag.WaveletTree Fibres
..returns:Fibre type.
..remarks:Some containers, such as @Spec.FMIndex@, can be seen as a bundle consisting of various fibres. Because not 
every table is a fibre we did not call them tables, however, in many cases one can think of fibres as tables. The 
fibre interface was designed to unify the access to the members of the different fibres.
To get a reference or the type of a specific fibre use @Function.getFibre@ or @Metafunction.Fibre@.		
..include:seqan/index.h
*/
/*!
 * @class WaveletTree
 *
 * @extends RankDictionary
 * 
 * @headerfile seqan/index.h
 * 
 * @brief A wavelet tree is a tree like binary encoding of a text.
 * 
 * @signature template <typename TValue, typename TSpec>
 *            RankDictionary<TValue, WaveletTree<TSpec> >
 * 
 * @tparam TValue The alphabet type of the wavelet tree.
 * @tparam TSpec A tag for specialization purposes. Default: <tt>void</tt>
 * 
 * @section Remarks
 * 
 * The nodes of a wavelet tree consist of a bit string as well as a character c.
 * In each level of the tree, characters smaller than c are represented as a 0
 * while character greater or equal to c are represented with a 1. The
 * characters represented by a 0 form the string to be represented by the left
 * subtree while characters represented by a 1 form the string of the right
 * subtree. Therefore, only the bit string of the root node represents all
 * characters while all other nodes represent subsets.
 */



template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreRanks>
{
    typedef StringSet<RankDictionary<bool, TwoLevels<TSpec> > > Type;
};

template <typename TValue, typename TSpec>
struct Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreTreeStructure>
{
    typedef typename MakeUnsigned<TValue>::Type TUChar_;
    typedef RightArrayBinaryTree<TUChar_, void>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this - const version should be Value const (as by default)

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, WaveletTree<TSpec> > >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<RankDictionary<TValue, WaveletTree<TSpec> > const> :
    public Value<RankDictionary<TValue, WaveletTree<TSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBlock_<TValue, WaveletTree<TSpec> >
{
    typedef RankDictionary<TValue, WaveletTree<TSpec> >             TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   TSize_;
//    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>                 Type;
    typedef String<TSize_>                                          Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Spec WaveletTree
// ----------------------------------------------------------------------------

/**
.Spec.WaveletTree:
..cat:Index
..summary:A wavelet tree is a tree like binary encoding of a text.
..signature:WaveletTree<TValue>
..param.TValue:The value type of the wavelet tree.
..include:seqan/index.h
..remarks:The nodes of a wavelet tree consist of a bit string as well as a character c. In each level of the tree, 
characters smaller than c are represented as a 0 while character greater or equal to c are represented with a 1.
The characters represented by a 0 form the string to be represented by the left subtree while characters represented
by a 1 form the string of the right subtree. Therefore, only the bit string of the root node represents all characters while all other nodes represent subsets.
*/

template <typename TValue, typename TSpec>
struct RankDictionary<TValue, WaveletTree<TSpec> >
{
    typename Fibre<RankDictionary, FibreRanks>::Type            ranks;
    typename Fibre<RankDictionary, FibreTreeStructure>::Type    waveletTreeStructure;

    RankDictionary() {}

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#getFibre:
..summary:Returns a specific fibre of a dictionary.
..signature:getFibre(dictionary, fibreTag)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary holding the fibre.
...type:Class.RankDictionary
..param.fibreTag:A tag that identifies the @Metafunction.Fibre@.
...type:Tag.WaveletTree Fibres
..returns:A reference to the @Metafunction.Fibre@ object.
..include:seqan/index.h
*/

template <typename TValue, typename TSpec>
inline typename Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreTreeStructure>::Type &
getFibre(RankDictionary<TValue, WaveletTree<TSpec> > & dict, FibreTreeStructure)
{
    return dict.waveletTreeStructure;
}

template <typename TValue, typename TSpec>
inline typename Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreTreeStructure>::Type const &
getFibre(RankDictionary<TValue, WaveletTree<TSpec> > const & dict, FibreTreeStructure)
{
    return dict.waveletTreeStructure;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#clear
..class:Class.RankDictionary
..summary:Clears the rank dictionary.
..signature:clear(dictionary)
..param.dictionary:The rank dictionary to be cleared.
...type:Class.RankDictionary
..include:seqan/index.h
*/

template <typename TValue, typename TSpec>
inline void clear(RankDictionary<TValue, WaveletTree<TSpec> > & dict)
{
    clear(getFibre(dict, FibreRanks()));
    clear(getFibre(dict, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#empty
..class:Class.RankDictionary
..summary:Returns whether or not the rank dictionary is empty.
..signature:empty(dictionary)
..param.dictionary:The rank dictionary to be checked.
...type:Class.RankDictionary
..returns:$true$ if the dictionary is empty, $false$ otherwise.
..include:seqan/index.h
*/

template <typename TValue, typename TSpec>
inline bool empty(RankDictionary<TValue, WaveletTree<TSpec> > const & dict)
{
    return empty(getFibre(dict, FibreRanks())) && empty(getFibre(dict, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#getValue
..summary:Returns the character of a specified position.
..signature:getCharacter(dictionary, pos)
..class:Class.RankDictionary
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.pos:The position
..include:seqan/index.h
*/

template <typename TValue, typename TSpec, typename TPos>
inline TValue getValue(RankDictionary<TValue, WaveletTree<TSpec> > & dict, TPos pos)
{
    typedef typename Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreTreeStructure>::Type const    TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type                 TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                                       TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                                     TChar;

    unsigned treePos = 0;
    typename Iterator<TWaveletTreeStructure, TopDown<> >::Type iter(dict.waveletTreeStructure, treePos);

    // initialize the return value with the smallest possible value
    TChar character = dict.waveletTreeStructure.minCharValue;

    // while the current node is not a leaf, go right if the bit at the current position is 1
    // go left otherwise
    // note that when going right the return value changes
    while (true)
    {
        TPos rank1 = getRank(dict.ranks[treePos], pos);
        if (getValue(dict.ranks[treePos], pos))
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

    return character;
}

template <typename TValue, typename TSpec, typename TPos>
inline TValue getValue(RankDictionary<TValue, WaveletTree<TSpec> > const & dict, TPos pos)
{
    return getValue(const_cast<RankDictionary<TValue, WaveletTree<TSpec> > &>(dict), pos);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#getRank:
..summary:Returns the rank (number of occurrences) of a specified character up to a specified position. 
..signature:getRank(dictionary, pos, character)
..class:Class.RankDictionary
..cat:Index
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.character:The character of interest.
..param.pos:The position (which is also included in the rank computation).
..returns:The rank (number of occurrences) of a specified character up to a specified position. 
...type:nolink:$unsigned$
..include:seqan/index.h
*/

template <typename TValue, typename TSpec, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, WaveletTree<TSpec> > >::Type
getRank(RankDictionary<TValue, WaveletTree<TSpec> > const & dict, TPos pos, TChar character)
{
    typedef typename Fibre<RankDictionary<TValue, WaveletTree<TSpec> >, FibreTreeStructure>::Type  TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type         TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                               TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                             TChar_;

    TPos sum = pos;
    TPos treePos = 0;

    // determine the leaf containing the character
    // count the number of 1 or 0 up to the computed position
    typename Iterator<TWaveletTreeStructure const, TopDown<> >::Type it(dict.waveletTreeStructure, treePos);
    TChar_ charInTree = dict.waveletTreeStructure.minCharValue;

    while (true)
    {
        TPos addValue = getRank(dict.ranks[treePos], sum);
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

    if (ordEqual(charInTree, character))
        return sum + 1;

    return 0;
}

// ----------------------------------------------------------------------------
// Function _fillStructure()
// ----------------------------------------------------------------------------

// This function is used to fill the bit strings of the wavelet tree.
template <typename TValue, typename TSpec, typename TText>
inline void _fillStructure(RankDictionary<TValue, WaveletTree<TSpec> > & dict, TText const & text)
{
    typedef RankDictionary<TValue, WaveletTree<TSpec> >                 TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type   TWaveletTreeStructure;
    typedef typename Iterator<TWaveletTreeStructure, TopDown<> >::Type  TWaveletTreeIterator;
    typedef typename Size<TRankDictionary>::Type                        TSize;
    typedef typename Iterator<TText const, Standard>::Type              TTextIterator;

    resize(dict.ranks, length(dict.waveletTreeStructure), Exact());

    for (TSize i = 0; i < length(dict.ranks); ++i)
        resize(dict.ranks[i], 0);

    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());

    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
    {
        TWaveletTreeIterator it(dict.waveletTreeStructure, 0);

        while (true)
        {
            // decide whether the character is smaller then the pivot element of the current node 
            if (ordGreater(getCharacter(it), value(textIt)))
            {
                // TODO(esiragusa): use resize() & setValue() instead of appendValue().
                appendValue(dict.ranks[getPosition(it)], false);
                if (!goLeftChild(it))
                    break;
            }
            else
            {
                appendValue(dict.ranks[getPosition(it)], true);
                if (!goRightChild(it))
                    break;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline void updateRanks(RankDictionary<TValue, WaveletTree<TSpec> > & dict)
{
    typedef RankDictionary<TValue, WaveletTree<TSpec> >                 TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    for (TSize i = 0; i < length(getFibre(dict, FibreRanks())); ++i)
        updateRanks(getFibre(dict, FibreRanks())[i]);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()
// ----------------------------------------------------------------------------
/**
.Function.RankDictionary#createRankDictionary
..class:Class.RankDictionary
..summary:This functions creates the dictionary.
..signature:createRankDictionary(dictionary, text)
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.text:A text to be transfered into a wavelet tree.
...type:Class.String
..include:seqan/index.h
*/

template <typename TValue, typename TSpec, typename TText, typename TPrefixSums>
inline void
createRankDictionary(RankDictionary<TValue, WaveletTree<TSpec> > & dict, TText const & text, TPrefixSums const & sums)
{
    createRightArrayBinaryTree(getFibre(dict, FibreTreeStructure()), sums);
//    _resizeStructure(dict, text);
    _fillStructure(dict, text);
    updateRanks(dict);
}

template <typename TValue, typename TSpec, typename TText>
inline void
createRankDictionary(RankDictionary<TValue, WaveletTree<TSpec> > & dict, TText const & text)
{
    typename RankDictionaryBlock_<TValue, WaveletTree<TSpec> >::Type sums;
    prefixSums<TValue>(sums, text);
    createRankDictionary(dict, text, sums);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#open
..class:Class.RankDictionary
..summary:This functions loads a dictionary from disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TValue, typename TSpec>
inline bool open(RankDictionary<TValue, WaveletTree<TSpec> > & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".wtc");
    if (!open(getFibre(dict, FibreRanks()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".wts");
    if (!open(getFibre(dict, FibreTreeStructure()), toCString(name), openMode)) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/**
.Function.RankDictionary#save
..class:Class.RankDictionary
..summary:This functions saves a dictionary to disk.
..signature:open(dictionary, fileName [, openMode])
..param.dictionary:The dictionary.
...type:Class.RankDictionary
..param.fileName:C-style character string containing the file name.
..param.openMode:The combination of flags defining how the file should be opened.
...remarks:To open a file read-only, write-only or to read and write use $OPEN_RDONLY$, $OPEN_WRONLY$, or $OPEN_RDWR$.
...remarks:To create or overwrite a file add $OPEN_CREATE$.
...remarks:To append a file if existing add $OPEN_APPEND$.
...remarks:To circumvent problems, files are always opened in binary mode.
...default:$OPEN_RDWR | OPEN_CREATE | OPEN_APPEND$
..returns:A $bool$ which is $true$ on success.
..include:seqan/index.h
*/

template <typename TValue, typename TSpec>
inline bool save(RankDictionary<TValue, WaveletTree<TSpec> > const & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".wtc");
    if (!save(getFibre(dict, FibreRanks()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".wts");
    if (!save(getFibre(dict, FibreTreeStructure()), toCString(name), openMode)) return false;

    return true;
}

}
#endif  // INDEX_FM_RANKDICTIONARY_WT
