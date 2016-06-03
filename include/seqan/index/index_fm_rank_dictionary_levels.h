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

template <typename TValue, typename TSpec>
struct RankDictionarySuperBlock_;

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

template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS_ = 1>
struct LevelsPrefixRDConfig : RDConfig<TSize, TFibre>
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

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >
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

template <typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef RankDictionary<bool, Levels<TSpec, TConfig> >         TRankDictionary_;

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
    typedef Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreRanks> TRanks_;

    typedef String<TEntry_, TFibreSpec_>                            Type;

    typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type superBlockValues;
    typename TRanks_::Type    blocks;
};

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperRanks>
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary_;
    typedef Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreRanks> TRanks_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type TFibreSpec_;

    typedef String<TRanks_, TFibreSpec_>                            Type;
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
// Metafunction RankDictionaryBitMask_
// ----------------------------------------------------------------------------

// TODO: simplify
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

template <typename TWord, typename TValue, typename TSpec, typename TConfig>
constexpr TWord _bitmaskWrapper(RankDictionary<TValue, Levels<TSpec, TConfig> > & /*dict*/, unsigned const bitsTotal, unsigned blocks, unsigned const constBlocks, unsigned const blocksize, TWord const vDefault, TWord const /*vPrefix*/)
{
    return _bitmask2<TWord>(bitsTotal, blocks, constBlocks, blocksize, vDefault);
}

template <typename TWord, typename TValue, typename TSpec, typename TSize, typename TFibre>
constexpr TWord _bitmaskWrapper(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > & /*dict*/, unsigned const bitsTotal, unsigned blocks, unsigned const constBlocks, unsigned const blocksize, TWord const /*vDefault*/, TWord const vPrefix)
{
    return _bitmask2<TWord>(bitsTotal, blocks, constBlocks, blocksize, vPrefix);
}

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

template <typename TValue, typename TConfig>
struct MyBitsPerValue
{
    static const typename BitsPerValue<TValue>::Type VALUE = BitsPerValue<TValue>::VALUE;
};

template <typename TValue, typename TSize, typename TFibre, unsigned LEVELS_>
struct MyBitsPerValue<TValue, LevelsPrefixRDConfig<TSize, TFibre, LEVELS_> >
{
    static const typename BitsPerValue<TValue>::Type VALUE = BitsPerValue<TValue>::VALUE + 1;
};

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, Levels<TSpec, TConfig> >
{
    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------

    static const unsigned _BITS_PER_VALUE   = MyBitsPerValue<TValue, TConfig>::VALUE;
    static const unsigned _BITS_PER_BLOCK   = 64;//RankDictionaryBitsPerBlock_<TValue, Levels<TSpec, TConfig> >::VALUE;
    static const unsigned _BITS_PER_WORD    = Min<RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::VALUE, _BITS_PER_BLOCK>::VALUE;
    static const unsigned _VALUES_PER_WORD  = _BITS_PER_WORD  / _BITS_PER_VALUE;
    static const unsigned _WORDS_PER_BLOCK  = _BITS_PER_BLOCK / _BITS_PER_WORD;
    static const unsigned _VALUES_PER_BLOCK = _VALUES_PER_WORD * _WORDS_PER_BLOCK;
    static const unsigned _VALUES_PER_SUPERBLOCK = _VALUES_PER_BLOCK * 2;

    typedef typename std::conditional<_BITS_PER_WORD == 64, uint64_t, uint32_t>::type TWordType;

    // TODO: rename bitmasks
    static TWordType _BITMASKS[ValueSize<TValue>::VALUE];
    static TWordType _NEWBITMASKS[_VALUES_PER_WORD];

    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreSuperRanks>::Type    superblocks;
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
            _BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, _VALUES_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i, i + (1 << (_BITS_PER_VALUE-1)));

        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _NEWBITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, i+1, i+1, _BITS_PER_VALUE, 1, 1 << (_BITS_PER_VALUE-1)); // 1
    }

    template <typename TText>
    RankDictionary(TText const & text) :
        _length(0)
    {
        auto maxValue = (1 << _BITS_PER_VALUE) - 1;
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            _BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, _VALUES_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i, i + (1 << (_BITS_PER_VALUE-1)));
        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _NEWBITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, i+1, i+1, _BITS_PER_VALUE, 1, 1 << (_BITS_PER_VALUE-1)); // 1

        createRankDictionary(*this, text);
    }
};

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_BITMASKS[ValueSize<TValue>::VALUE];

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_NEWBITMASKS[RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD];








template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperRanks>::Type &
getFibre(RankDictionary<TValue, Levels<TSpec, Levels<TSpec, TConfig> > > & dict, FibreSuperRanks)
{
    return dict.superblocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperRanks>::Type const &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, FibreSuperRanks)
{
    return dict.superblocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline void clear(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict)
{
    for (unsigned i = 0; i < length(dict.superblocks); ++i)
    {
        clear(dict.superblocks[i].blocks);
    }
    clear(dict.superblocks);
}



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
// Function _toSuperBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
 _toSuperBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK;
}

// ----------------------------------------------------------------------------
// Function _toBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return (pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK) / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
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

template <typename TValue, typename TSpec, typename TConfig, typename TSuperBlockPos, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSuperBlockPos superBlockPos, TBlockPos blockPos, TWordPos wordPos)
{
    return dict.superblocks[superBlockPos].blocks[blockPos].values[wordPos];
}

template <typename TValue, typename TSpec, typename TConfig, typename TSuperBlockPos, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TSuperBlockPos superBlockPos, TBlockPos blockPos, TWordPos wordPos)
{
    return dict.superblocks[superBlockPos].blocks[blockPos].values[wordPos];
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)].blocks[_toBlockPos(dict, pos)].values;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)].blocks[_toBlockPos(dict, pos)].values;
}

// ----------------------------------------------------------------------------
// Function _blockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)].blocks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)].blocks[_toBlockPos(dict, pos)].block;
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
    //TSize endPos   = length(dict.ranks) * TRankDictionary::_VALUES_PER_BLOCK
    TSize endPos = 0;
    for (unsigned i = 0; i < length(dict.superblocks); ++i)
        endPos += length(dict.superblocks[i].blocks);
    endPos *= TRankDictionary::_VALUES_PER_BLOCK;

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

// TODO: prototype
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
_getBlockRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & /* dict */, TBlock const & block, TPos /* pos */, TChar c, TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = block[ordValue(c)-1];
    smaller += _smaller; // TODO: _smaller cannot be removed. order of evaluation is not defined!
    return block[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TBlock const & block, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? block : pos - _toPosInBlock(dict, pos) - block;
}

// ----------------------------------------------------------------------------
// Function _getSuperBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSuperBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getSuperBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TSuperBlock const & superblock, TPos /* pos */, TChar c)
{
    return superblock.superBlockValues[ordValue(c)];
}

// TODO: prototype
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TSuperBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
 _getSuperBlockRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & /* dict */, TSuperBlock const & superblock, TPos /* pos */, TChar c, TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = superblock[ordValue(c)-1];
    smaller += _smaller; // TODO: _smaller cannot be removed. order of evaluation is not defined!
    return superblock.superBlockValues[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TSuperBlock, typename TPos>
inline typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type
_getSuperBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TSuperBlock const & superblock, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    // TODO: richtig?
    return c ? superblock.superBlockValues : pos - _toPosInBlock(dict, pos) - superblock.superBlockValues;
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
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const &,
             TWord const & values,
             TPosInWord posInWord,
             TValue c)
{
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >                TRankDictionary;

    return popCount((TRankDictionary::_BITMASKS[ordValue(c)] - values) & TRankDictionary::_NEWBITMASKS[posInWord]);
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TWord, typename TPosInWord, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const &,
             TWord const & values,
             TPosInWord posInWord,
             TValue c,
             TPos & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >                TRankDictionary;

    auto _smaller = popCount((TRankDictionary::_BITMASKS[ordValue(c)-1] - values) & TRankDictionary::_NEWBITMASKS[posInWord]);
    smaller += _smaller;
    return popCount((TRankDictionary::_BITMASKS[ordValue(c)] - values) & TRankDictionary::_NEWBITMASKS[posInWord]) - _smaller;
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

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TValues, typename TPosInBlock, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
_getValueRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict,
              TValues const & values,
              TPosInBlock posInBlock,
              TValue c,
              TSmaller & smaller)
{
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > TRankDictionary;

    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    TSize valueRank = 0;

    // NOTE(esiragusa): writing the loop in this form prevents the compiler from unrolling it.
    //    for (TSize wordPrevPos = 0; wordPrevPos < wordPos; ++wordPrevPos)
    //      valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    for (TSize wordPrevPos = 0; wordPrevPos < TRankDictionary::_WORDS_PER_BLOCK; ++wordPrevPos)
        if (wordPrevPos < wordPos) valueRank += _getWordRank(dict, values[wordPrevPos].i, TRankDictionary::_VALUES_PER_WORD - 1, c, smaller);

    valueRank += _getWordRank(dict, values[wordPos].i, posInWord, c, smaller);

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

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >::Type
_getValuesRanks(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos)
{
    typedef typename RankDictionaryBlock_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                                        TValueSize;

    TBlock blockRank;

    TPos rank0 = _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue(0));
    assignValue(blockRank, 0, rank0);

    for (TValueSize c = 1; c < ValueSize<TValue>::VALUE; ++c)
    {
        // TODO: only temporary workaround with (uint16_t) cast
        TPos smaller = 0;
        TPos rank = _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue((uint16_t) c), smaller);
        assignValue(blockRank, c, rank + smaller);
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

// TODO: needed because otherwise it template specialization would be ambiguous. maybe try to remove all bool-specializations?
template <typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename RankDictionaryBlock_<bool, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >::Type
        _getValuesRanks(RankDictionary<bool, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos)
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
    //typedef typename Fibre<TRankDictionary, FibreSuperRanks>::Type          TFibreSuperRanks;
    //typedef typename Value<FibreSuperRanks>::Type                           TFibreRank;


    //typedef typename Fibre<TRankDictionary, FibreSuperRanks>::Type          TFibreRank;
    //typedef typename Value<TFibreRanks>::Type                               TRankEntry;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos      = _toBlockPos(dict, pos);
    TSize posInBlock    = _toPosInBlock(dict, pos);

    auto const & superblock = dict.superblocks[superBlockPos];
    auto const & entry = superblock.blocks[blockPos];

    return _getSuperBlockRank(dict, superblock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos, TValue c, TPos & smaller)
{
    //typedef RankDictionary<TValue, Levels<TSpec , LevelsPrefixRDConfig<TSize, TFibre> > > const           TRankDictionary;
    //typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    //typedef typename Value<TFibreRanks>::Type                               TRankEntry;

    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);

    auto const & superblock = dict.superblocks[superBlockPos];
    auto const & entry = superblock.blocks[blockPos];

    smaller = 0;
    if (ordValue(c) > 0)
        return _getSuperBlockRank(dict, superblock, pos, static_cast<TValue>(c))
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);

    // c == Dna('A')
    return _getSuperBlockRank(dict, superblock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
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

    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos      = _toBlockPos(dict, pos);
    TSize posInBlock    = _toPosInBlock(dict, pos);
    TSize wordPos       = _toWordPos(dict, posInBlock);
    TSize posInWord     = _toPosInWord(dict, posInBlock);

    return _valuesAt(dict, superBlockPos, blockPos, wordPos)[posInWord];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
getValue(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >             TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    return _valuesAt(dict, superBlockPos, blockPos, wordPos)[posInWord];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > >::Type
getValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > & dict, TPos pos)
{
    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    auto SIZE = Size<typename RankDictionaryValues_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >::TValues>::VALUE;
    unsigned shift = (SIZE - posInWord - 1) * (BitsPerValue<TValue>::VALUE + 1);
    auto value = _valuesAt(dict, superBlockPos, blockPos, wordPos);

    return (value.i >> shift) & value.BIT_MASK2;
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const>::Type
getValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos)
{
    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    auto SIZE = Size<typename RankDictionaryValues_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > >::TValues>::VALUE;
    unsigned shift = (SIZE - posInWord - 1) * (BitsPerValue<TValue>::VALUE + 1);
    auto value = _valuesAt(dict, superBlockPos, blockPos, wordPos);

    return (value.i >> shift) & value.BIT_MASK2;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline void setValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos, TChar c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >             TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    assignValue(_valuesAt(dict, superBlockPos, blockPos, wordPos), posInWord, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos, typename TChar>
inline void setValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > & dict, TPos pos, TChar c)
{
    TSize superBlockPos = _toSuperBlockPos(dict, pos);
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    assignValue2(_valuesAt(dict, superBlockPos, blockPos, wordPos), posInWord, static_cast<TValue>(c));
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
    typedef typename Fibre<TRankDictionary, FibreSuperRanks>::Type       TFibreSuperBlocks;
    typedef typename Iterator<TFibreSuperBlocks, Standard>::Type          TSuperBlockIter;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type       TFibreBlocks;
    typedef typename Iterator<TFibreBlocks, Standard>::Type          TBlockIter;

    typedef typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type TSuperBlock;

    if (empty(dict)) return;

    TSuperBlockIter superBlocksBegin = begin(dict.superblocks, Standard());
    TSuperBlockIter superBlocksEnd = end(dict.superblocks, Standard());

    // Insures the first block ranks start from zero.
    _clearBlockAt(dict, 0u);

    // Clear the uninitialized values.
    _padValues(dict);

    TSuperBlock superBlockSum;// = t0;
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        superBlockSum[i] = 0;
    // Iterate through the blocks.
    std::cout << length(dict.superblocks) << std::endl;
    TSuperBlockIter superBlocksIt = superBlocksBegin;
    for (; superBlocksIt != superBlocksEnd/* - 1*/; ++superBlocksIt) // TODO ?
    {
        //TSize superBlockPos = superBlockIt - superBlocksBegin;
        auto & superBlock = dict.superblocks[superBlocksIt - superBlocksBegin];
        TBlockIter blocksBegin = begin(superBlock.blocks, Standard());
        TBlockIter blocksEnd = end(superBlock.blocks, Standard());

        superBlock.superBlockValues = superBlockSum;

        // TODO: wirklich (un)nötig?
        if (blocksBegin != blocksEnd)
        {
            //Tuple<TSize, ValueSize<TValue>::VALUE> _t0;
            unsigned blocks_per_superblock = TRankDictionary::_VALUES_PER_SUPERBLOCK / TRankDictionary::_VALUES_PER_BLOCK;
            for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
                _blockAt(dict, _toPos(dict, blocks_per_superblock * (superBlocksIt - superBlocksBegin)))[i] = 0; // TRankDictionary::_VALUES_PER_SUPERBLOCK * (superBlocksIt - superBlocksBegin)
        }
        else
            SEQAN_ASSERT(false);

        TSize next, curr;
        for (TBlockIter blocksIt = blocksBegin; blocksIt != blocksEnd - 1; ++blocksIt) // TODO ?
        {
            unsigned blocks_per_superblock = TRankDictionary::_VALUES_PER_SUPERBLOCK / TRankDictionary::_VALUES_PER_BLOCK;
            TSize blockPos = blocksIt - blocksBegin + (blocks_per_superblock * (superBlocksIt - superBlocksBegin)); //  + TRankDictionary::_VALUES_PER_SUPERBLOCK * (superBlocksIt - superBlocksBegin)
            curr = _toPos(dict, blockPos);
            next = _toPos(dict, blockPos + 1);

            _blockAt(dict, next) = _blockAt(dict, curr) + _getValuesRanks(dict, next - 1);
        }
        // TODO kann _blockAt(dict, next) hier undef. sein, wenn er nicht in die schleife reingeht?
        if (blocksBegin != blocksEnd - 1)
            superBlockSum = superBlockSum + (_blockAt(dict, curr) + _getValuesRanks(dict, next - 1));
        //else
        //    SEQAN_ASSERT(false);
    }
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos, TValue c, TPos & smaller)
{
    smaller = 0;
    if (ordValue(c) == 0)
        return getRank(dict, pos, c);
    return getRank(dict, pos, c, smaller);
}

// TODO: prototype for prefix sums for DNA (used by bidirectional FM index)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre> > > const & dict, TPos pos, TValue c)
{
    TPos smaller;
    if (ordValue(c) == 0)
        return getRank(dict, pos, c);
    return getRank(dict, pos, c, smaller);
}

// TODO: what is this for? wrapper?
template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
getCumulativeRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TValue c, TPos & /*smaller*/)
{
    // not cumulative!!!!
    return getRank(dict, pos, c);
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

/*template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
reserve(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
    return reserve(dict.ranks, (newCapacity + RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK - 1) /
                               RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK, tag);
}*/

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
resize(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize newLength, Tag<TExpand> const tag)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDict_;

    dict._length = newLength;
    //auto values = newLength + TRankDict_::_VALUES_PER_BLOCK - 1;
    auto superblocks = (newLength + TRankDict_::_VALUES_PER_SUPERBLOCK - 1) / TRankDict_::_VALUES_PER_SUPERBLOCK; // eq. to ceil(newLength / VALUES_PER_BLOCK)
    auto ret1 = resize(dict.superblocks, superblocks, tag);

    for (unsigned i = 0; i < superblocks-1; ++i)
    {
        resize(dict.superblocks[i].blocks, TRankDict_::_VALUES_PER_SUPERBLOCK / TRankDict_::_VALUES_PER_BLOCK, tag);
    }
    // last superblock might have fewer blocks
    resize(dict.superblocks[superblocks-1].blocks, ((newLength % TRankDict_::_VALUES_PER_SUPERBLOCK) + TRankDict_::_VALUES_PER_BLOCK - 1) / TRankDict_::_VALUES_PER_BLOCK, tag);
    // eq. to ceil((newLength % TRankDict_::_VALUES_PER_SUPERBLOCK) / VALUES_PER_BLOCK)

    return ret1;
    // old
    //return resize(dict.ranks, (newLength + RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK - 1) /
    //                          RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK, tag);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
