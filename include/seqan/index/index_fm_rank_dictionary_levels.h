// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_RANK_DICTIONARY_LEVELS_H_
#define INDEX_FM_RANK_DICTIONARY_LEVELS_H_
#include <bitset>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryWordSize_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryWordSize_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitsPerBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBitsPerBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValues_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryValues_;

// ----------------------------------------------------------------------------
// Struct RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryEntry_;

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag LevelsRDConfig
// ----------------------------------------------------------------------------

template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS_ = 1>
struct LevelsRDConfig : RDConfig<TSize, TFibre>
{
    static const unsigned LEVELS =  LEVELS_;
};

// ----------------------------------------------------------------------------
// Tag Levels
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = LevelsRDConfig<> >
struct Levels {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitMask_
// ----------------------------------------------------------------------------

// TODO: merge into one method. clean-up! constBlocks
template <typename TWord>
constexpr TWord _bitmask2(unsigned const bitsTotal, unsigned blocks, unsigned const constBlocks, unsigned const blocksize, TWord const value)
{
    return (blocks == constBlocks) ?
           _bitmask2<TWord>(bitsTotal, blocks - 1, constBlocks, blocksize, value << (bitsTotal - blocksize)) :
           (
               (blocks == 0) ?
               value >> (bitsTotal % blocksize) :
               _bitmask2<TWord>(bitsTotal, blocks - 1, constBlocks, blocksize, value | (value >> blocksize))
           );
}

// TODO: RankDictionaryBitMask_ is deprecated!
template <typename TValue>
struct RankDictionaryBitMask_;

template <>
struct RankDictionaryBitMask_<uint32_t>
{
    static const uint32_t VALUE = 0x55555555;
};

template <>
struct RankDictionaryBitMask_<uint64_t>
{
    static const uint64_t VALUE = 0x5555555555555555ull;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryWordSize_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> > :
    BitsPerValue<uint64_t> {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitsPerBlock_
// ----------------------------------------------------------------------------
// The number of bits per block equals the number of bits of the block summary.

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBitsPerBlock_<TValue, Levels<TSpec, TConfig> > :
    BitsPerValue<typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type> {};

// NOTE(esiragusa): This lets a Dna block to have the size of one word - one popcount per block.
//template <typename TSpec, typename TConfig>
//struct RankDictionaryBitsPerBlock_<Dna, Levels<TSpec, TConfig> > :
//    RankDictionaryWordSize_<Dna, Levels<TSpec, TConfig> > {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   TSize_;

    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>                 Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionaryBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef RankDictionary<bool, Levels<TSpec, TConfig> >           TRankDictionary_;

    typedef typename Size<TRankDictionary_>::Type                   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValues_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >                 TRankDictionary_;

    typedef Tuple<TValue, TRankDictionary_::_VALUES_PER_WORD, BitPacked<> > TValues;
    typedef typename TValues::TBitVector                                    TWord;
    typedef Tuple<TValues, TRankDictionary_::_WORDS_PER_BLOCK>              Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreRanks>
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary_;
    typedef RankDictionaryEntry_<TValue, Levels<TSpec, TConfig> >   TEntry_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type TFibreSpec_;

    typedef String<TEntry_, TFibreSpec_>                            Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Struct RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec = Levels<> >
struct RankDictionaryEntry_ {};

// ----------------------------------------------------------------------------
// Struct Levels RankDictionaryEntry_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryEntry_<TValue, Levels<TSpec, TConfig> >
{
    // A bit-compressed block of TValue symbols.
    typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type   values;

    // A summary of counts for each block of TValue symbols.
    typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type    block;
};

// ----------------------------------------------------------------------------
// Class Levels RankDictionary
// ----------------------------------------------------------------------------
// TODO(esiragusa): update doc
/*!
 * @class TwoLevelRankDictionary
 * @extends RankDictionary
 * @headerfile <seqan/index.h>
 *
 * @brief A TwoLevelRankDictionary is a @link RankDictionary @endlink consisting of two levels.
 *
 * @signature template <typename TValue, typename TSpec, typename TConfig>
 *            class RankDictionary<TValue, WaveletTree<TSpec> >;
 *
 * @tparam TValue The alphabet type of the wavelet tree.
 * @tparam TSpec  A tag for specialization purposes. Default: <tt>void</tt>
 *
 * This @link RankDictionary @endlink consists of two levels of rank
 * infromation, in which one stores information of blocks and the other
 * information until a specified block. Combining those two informations
 * leads to constant rank dictionary look ups.
 */

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, Levels<TSpec, TConfig> >
{
    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------

    static const unsigned _BITS_PER_VALUE   = BitsPerValue<TValue>::VALUE;
    static const unsigned _BITS_PER_BLOCK   = RankDictionaryBitsPerBlock_<TValue, Levels<TSpec, TConfig> >::VALUE;
    static const unsigned _BITS_PER_WORD    = Min<RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::VALUE, _BITS_PER_BLOCK>::VALUE;
    static const unsigned _VALUES_PER_WORD  = _BITS_PER_WORD  / _BITS_PER_VALUE;
    static const unsigned _WORDS_PER_BLOCK  = _BITS_PER_BLOCK / _BITS_PER_WORD;
    static const unsigned _VALUES_PER_BLOCK = _VALUES_PER_WORD * _WORDS_PER_BLOCK;

    typedef typename std::conditional<_BITS_PER_WORD == 64, uint64_t, uint32_t>::type TWordType;

    // TODO: rename bitmasks
    static TWordType _BITMASKS[ValueSize<TValue>::VALUE];
    static TWordType _NEWBITMASKS[_VALUES_PER_WORD];

    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreRanks>::Type    ranks;
    typename Size<RankDictionary>::Type                 _length;
    // TODO(esiragusa): open/save _length or remove it.

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    RankDictionary() :
        _length(0)
    {
        // TODO: no code duplication. make it a const-expr?
        auto maxValue = (1 << _BITS_PER_VALUE) - 1;
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            _BITMASKS[i] = _bitmask2<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i);

        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _NEWBITMASKS[i] = _bitmask2<TWordType>(_BITS_PER_WORD, i+1, i+1, _BITS_PER_VALUE, 1);
    }

    template <typename TText>
    RankDictionary(TText const & text) :
        _length(0)
    {
        auto maxValue = (1 << _BITS_PER_VALUE) - 1;
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            _BITMASKS[i] = _bitmask2<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i);
        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _NEWBITMASKS[i] = _bitmask2<TWordType>(_BITS_PER_WORD, i+1, i+1, _BITS_PER_VALUE, 1);

        createRankDictionary(*this, text);
    }
};

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_BITMASKS[ValueSize<TValue>::VALUE];

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_NEWBITMASKS[RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD];

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _toPosInWord()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInWord(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos posInBlock)
{
    return posInBlock % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD;
}

// ----------------------------------------------------------------------------
// Function _toWordPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toWordPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos posInBlock)
{
    return posInBlock / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD;
}

// ----------------------------------------------------------------------------
// Function _toPosInBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _toBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _toPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TBlockPos blockPos)
{
    return blockPos * RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TBlockPos blockPos, TWordPos wordPos)
{
    return dict.ranks[blockPos].values[wordPos];
}

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TBlockPos blockPos, TWordPos wordPos)
{
    return dict.ranks[blockPos].values[wordPos];
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].values;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].values;
}

// ----------------------------------------------------------------------------
// Function _blockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.ranks[_toBlockPos(dict, pos)].block;
}

// ----------------------------------------------------------------------------
// Function _padValues()
// ----------------------------------------------------------------------------
// Set values beyond length(dict) but still within the end of the ranks fibre.

template <typename TValue, typename TSpec, typename TConfig>
inline void _padValues(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;

    TSize beginPos = length(dict);
    TSize endPos   = length(dict.ranks) * TRankDictionary::_VALUES_PER_BLOCK;

    for (TSize pos = beginPos; pos < endPos; ++pos)
        setValue(dict, pos, TValue());
}

// ----------------------------------------------------------------------------
// Function _clearBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline void _clearBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    clear(_blockAt(dict, pos));
}

// ----------------------------------------------------------------------------
// Function _clearBlockAt(bool)
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
inline void _clearBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    _blockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _getBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TBlock const & block, TPos /* pos */, TChar c)
{
    return block[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TConfig, typename TBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TBlock const & block, TPos /* pos */, TChar c, TSmaller smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    // TODO: can _smaller be left out? I.e. can smaller be unequal zero when _getBlockRank is called?
    TSmaller _smaller = block[ordValue(c)-1];
    smaller += _smaller;
    return block[ordValue(c)] - _smaller;
}

// ----------------------------------------------------------------------------
// Function _getBlockRank(bool)
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TBlock const & block, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? block : pos - _toPosInBlock(dict, pos) - block;
}

// ----------------------------------------------------------------------------
// Function _getWordRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */,
             TWord const & word,
             TPosInWord posInWord,
             TValue c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >                TRankDictionary;

    TWord mask = word ^ TRankDictionary::_BITMASKS[ordValue(c)];

    // NOTE: actually it should be: mask & mask >> 1 & mask >> 2 & ... but this is easier and equivalent
    for (TWord i = 1; i < TRankDictionary::_BITS_PER_VALUE; ++i)
        mask &= mask >> 1;

    return popCount(TRankDictionary::_NEWBITMASKS[posInWord] & mask);
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TConfig, typename TWord, typename TPosInWord, typename TPos>
inline typename Size<RankDictionary<Dna, Levels<void, TConfig> > const>::Type
_getWordRank(RankDictionary<Dna, Levels<int, TConfig> > const &,
             TWord const & values,
             TPosInWord posInWord,
             Dna c,
             TPos & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    typedef RankDictionary<Dna, Levels<void, TConfig> >                TRankDictionary;

    TWord word = hiBits(values, TRankDictionary::_BITS_PER_VALUE * (posInWord + 1));

    TWord mA = ~(word | TRankDictionary::_BITMASKS[2]);
    TWord mC = (~(word | TRankDictionary::_BITMASKS[1])) << 1;

    TWord masks[4] = {};
    /*A*/ masks[0] = RankDictionaryBitMask_<TWord>::VALUE & (~word >> 1) & ~word;
    /*C*/ masks[1] = mA;
    /*G*/ masks[2] = mA|mC;
    /*T*/ masks[3] = RankDictionaryBitMask_<TWord>::VALUE;

    TPos rest = TRankDictionary::_VALUES_PER_WORD - (posInWord + 1);

    TPos _smaller = popCount(masks[ordValue(c)-1]);
    smaller += _smaller - rest;

    return popCount(masks[ordValue(c)]) - _smaller;
}

// TODO: original implementation by esirgusa (just for reference and performance comparison)
template <typename TConfig, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<Dna, Levels<bool, TConfig> > const>::Type
_getWordRank(RankDictionary<Dna, Levels<bool, TConfig> > const & /* dict */,
             TWord const & values,
             TPosInWord posInWord,
             Dna c)
{
    typedef RankDictionary<Dna, Levels<bool, TConfig> >                TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    // Clear the last positions.
    TWord word = hiBits(values, TRankDictionary::_BITS_PER_VALUE * (posInWord + 1));

    // And matches when c == G|T.
    TWord odd  = ((ordValue(c) & ordValue(Dna('G'))) ? word : ~word) >> 1;

    // And matches when c == C|T.
    TWord even = ((ordValue(c) & ordValue(Dna('C'))) ? word : ~word);

    // Apply the interleaved mask.
    TWord mask = odd & even & RankDictionaryBitMask_<TWord>::VALUE;

    // The rank is the sum of the bits on.
    TSize valueRank = popCount(mask);

    // If c == A then masked character positions must be subtracted from the count.
    if (c == Dna('A')) valueRank -= TRankDictionary::_VALUES_PER_WORD - (posInWord + 1);

    return valueRank;
}

// TODO: original implementation by esirgusa (just for reference and performance comparison)
template <typename TConfig, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<bool, Levels<bool, TConfig> > const>::Type
_getWordRank(RankDictionary<bool, Levels<bool, TConfig> > const & /* dict */,
             TWord const & values,
             TPosInWord posInWord,
             bool c)
{
    // Negate the values to compute the rank of zero.
    TWord word = c ? values : ~values;

    // Clear the last positions.
    TWord mask = hiBits(word, posInWord + 1);

    // Get the sum of the bits on.
    return popCount(mask);
}

// ----------------------------------------------------------------------------
// Function _getWordRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TWord const & values, TValue c)
{
    return _getWordRank(dict, values, RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD - 1, c);
}

// ----------------------------------------------------------------------------
// Function _getValueRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TValues, typename TPosInBlock>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getValueRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict,
              TValues const & values,
              TPosInBlock posInBlock,
              TValue c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDictionary;
    typedef typename Size<TRankDictionary>::Type            TSize;

    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    TSize valueRank = 0;

    // NOTE(esiragusa): writing the loop in this form prevents the compiler from unrolling it.
//    for (TSize wordPrevPos = 0; wordPrevPos < wordPos; ++wordPrevPos)
//      valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    for (TSize wordPrevPos = 0; wordPrevPos < TRankDictionary::_WORDS_PER_BLOCK; ++wordPrevPos)
        if (wordPrevPos < wordPos) valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    valueRank += _getWordRank(dict, values[wordPos].i, posInWord, c);

    return valueRank;
}

template <typename TValue, typename TConfig, typename TValues, typename TPosInBlock, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<int, TConfig> > const>::Type
_getValueRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict,
              TValues const & values,
              TPosInBlock posInBlock,
              TValue c,
              TSmaller & smaller)
{
    typedef RankDictionary<TValue, Levels<int, TConfig> > TRankDictionary;
    typedef typename Size<TRankDictionary>::Type            TSize;

    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    TSize valueRank = 0;

    // NOTE(esiragusa): writing the loop in this form prevents the compiler from unrolling it.
    //    for (TSize wordPrevPos = 0; wordPrevPos < wordPos; ++wordPrevPos)
    //      valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    for (TSize wordPrevPos = 0; wordPrevPos < TRankDictionary::_WORDS_PER_BLOCK; ++wordPrevPos)
        if (wordPrevPos < wordPos) valueRank += _getWordRank(dict, values[wordPrevPos].i, RankDictionary<TValue, Levels<int, TConfig> >::_VALUES_PER_WORD - 1, c, smaller);

    valueRank += _getWordRank(dict, values[wordPos].i, posInWord, c);

    return valueRank;
}

/*template <typename TValue, typename TConfig, typename TValues, typename TPosInBlock>
inline typename Size<RankDictionary<TValue, Levels<void, TConfig> > const>::Type
_getValueRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict,
              TValues const & values,
              TPosInBlock posInBlock,
              TValue c)
{
    TPosInBlock smaller;
    return _getValueRank(dict, values, posInBlock, c, smaller);
}*/

// ----------------------------------------------------------------------------
// Function _getValuesRanks()
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _getValuesRanks() for Dna.
// TODO: think about why one should do that? probably outdated comment

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type
_getValuesRanks(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    typedef typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                                        TValueSize;

    TBlock blockRank;

    // TODO: only temporary workaround with (uint16_t) cast
    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
        assignValue(blockRank, c, _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue((uint16_t) c)));

    return blockRank;
}

template <typename TValue, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<void, TConfig> >::Type
_getValuesRanks(RankDictionary<TValue, Levels<int, TConfig> > const & dict, TPos pos)
{
    typedef typename RankDictionaryBlock_<TValue, Levels<int, TConfig> >::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                                        TValueSize;

    TBlock blockRank;

    TPos smaller = 0;
    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
    {
        // TODO: only temporary workaround with (uint16_t) cast
        smaller += _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue((uint16_t) c));
        assignValue(blockRank, c, smaller);
    }

    return blockRank;
}

// ----------------------------------------------------------------------------
// Function _getValuesRanks(bool)
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<bool, Levels<TSpec, TConfig> >::Type
_getValuesRanks(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), true);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------
template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TChar c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const           TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Value<TFibreRanks>::Type                               TRankEntry;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);

    TRankEntry const & entry = dict.ranks[blockPos];

    return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c)) +
           _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<void, TConfig> > const>::Type
getRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict, TPos pos, Dna c, TPos & smaller)
{
    typedef RankDictionary<TValue, Levels<int, TConfig> > const           TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Value<TFibreRanks>::Type                               TRankEntry;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);

    TRankEntry const & entry = dict.ranks[blockPos];

    smaller = 0;
    if (c != Dna('A'))
        return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);

    // c == Dna('A')
    return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

// TODO: nicht möglich. auch bei kumulativer version wollen wir bei ordValue(c) == 0
// TODO: die nicht-kumulative version aufrufen, weil andernfalls beim kumulativen wrapper ...[ordValue(c)-1] aufgerufen werden würde
/*template <typename TValue, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Size<RankDictionary<TValue, Levels<void, TConfig> > const>::Type
getRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict, TPos pos, Dna c)
{
    TPos smaller;
    return getRank(dict, pos, c, smaller);
}*/

// ----------------------------------------------------------------------------
// Function getRank(bool)
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type
getRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return getRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------
template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
getValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >             TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    return _valuesAt(dict, blockPos, wordPos)[posInWord];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
getValue(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >             TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    return _valuesAt(dict, blockPos, wordPos)[posInWord];
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline void setValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos, TChar c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >             TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    assignValue(_valuesAt(dict, blockPos, wordPos), posInWord, static_cast<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Better not to have appendValue() - it is not efficient - and thus neither length().

template <typename TValue, typename TSpec, typename TConfig, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TChar c, Tag<TExpand> const tag)
{
    resize(dict, length(dict) + 1, tag);
    setValue(dict, length(dict) - 1, c);
}

// ----------------------------------------------------------------------------
// Function updateRanks()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline void updateRanks(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type       TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type          TFibreRanksIter;

    if (empty(dict)) return;

    TFibreRanksIter ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIter ranksEnd = end(dict.ranks, Standard());

    // Insures the first block ranks start from zero.
    _clearBlockAt(dict, 0u);

    // Clear the uninitialized values.
    _padValues(dict);

    // Iterate through the blocks.
    for (TFibreRanksIter ranksIt = ranksBegin; ranksIt != ranksEnd - 1; ++ranksIt)
    {
        TSize blockPos = ranksIt - ranksBegin;
        TSize curr = _toPos(dict, blockPos);
        TSize next = _toPos(dict, blockPos + 1);

        _blockAt(dict, next) = _blockAt(dict, curr) + _getValuesRanks(dict, next - 1);
    }
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<void, TConfig> > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict, TPos pos, Dna character, TPos & smaller)
{
    smaller = 0;
    if (character == Dna('A'))
        return getRank(dict, pos, character);
    return getRank(dict, pos, character, smaller);
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<void, TConfig> > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<int, TConfig> > const & dict, TPos pos, Dna character)
{
    TPos smaller;
    if (character == Dna('A'))
        return getRank(dict, pos, character);
    return getRank(dict, pos, character, smaller);
}

// TODO: what is this for? wrapper?
template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, Dna character, TPos & /*smaller*/)
{
    // not cumulative!!!!
    return getRank(dict, pos, character);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
length(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict)
{
    return dict._length;
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
reserve(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
    return reserve(dict.ranks, (newCapacity + RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK - 1) /
                               RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK, tag);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
resize(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize newLength, Tag<TExpand> const tag)
{
    dict._length = newLength;
    return resize(dict.ranks, (newLength + RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK - 1) /
                              RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK, tag);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
