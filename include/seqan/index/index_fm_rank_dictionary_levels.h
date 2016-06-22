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
// Metafunction RankDictionarySuperBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionarySuperBlock_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryUltraBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryUltraBlock_;

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
struct LevelsRDConfig : RDConfig<TSize, TFibre> {
    static const unsigned LEVELS = LEVELS_;
};

template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS_ = 1>
struct LevelsPrefixRDConfig : RDConfig<TSize, TFibre> {
    static const unsigned LEVELS = LEVELS_;
};

// ----------------------------------------------------------------------------
// Tag Levels
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = LevelsRDConfig<>>
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
// Metafunction LevelSize_
// ----------------------------------------------------------------------------
//       | Blocks | SuperBlocks  | UltraBlocks
// 1-lvl | 64 bit | ------------ | -----------
// 2-lvl | 32 bit |    64 bit    | -----------
// 3-lvl | 16 bit |    32 bit    |    64 bit

template <unsigned TotalLevels, unsigned Level>
struct LevelSize_
{
    typedef uint64_t Type;
};

template <>
struct LevelSize_<2, 1>
{
    typedef uint32_t Type;
};

template <>
struct LevelSize_<3, 2>
{
    typedef uint32_t Type;
};

template <>
struct LevelSize_<3, 1>
{
    typedef uint16_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >
{
    // biggest datatype, i.e. of topmost level
    typedef typename LevelSize_<TConfig::LEVELS, TConfig::LEVELS>::Type  Type;
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
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 1>::Type   TSize_;
    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>         Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionaryBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 1>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionarySuperBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 2>::Type   TSize_;
    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>         Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 2>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryUltraBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 3>::Type   TSize_;
    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>         Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionaryUltraBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename LevelSize_<TConfig::LEVELS, 3>::Type   Type;
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

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperBlocks>
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary_;
    typedef typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type TSuperBlocks_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type TFibreSpec_;

    typedef String<TSuperBlocks_, TFibreSpec_>                           Type;
};

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreUltraBlocks>
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary_;
    typedef typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type TUltraBlocks_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type TFibreSpec_;

    typedef String<TUltraBlocks_, TFibreSpec_>                           Type;
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
constexpr TWord _bitmaskWrapper(RankDictionary<TValue, Levels<TSpec, TConfig> > & /*dict*/, unsigned const bitsTotal, unsigned blocks, unsigned const blocksize, TWord const vDefault, TWord const /*vPrefix*/)
{
    return _bitmask2<TWord>(bitsTotal, blocks, blocks, blocksize, vDefault);
}

template <typename TWord, typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS>
constexpr TWord _bitmaskWrapper(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > & /*dict*/, unsigned const bitsTotal, unsigned blocks, unsigned const blocksize, TWord const /*vDefault*/, TWord const vPrefix)
{
    return _bitmask2<TWord>(bitsTotal, blocks, blocks, blocksize, vPrefix);
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

// TODO: do not allow prefix levels for bools

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, Levels<TSpec, TConfig> >
{
    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------

    static const unsigned _BITS_PER_VALUE   = MyBitsPerValue<TValue, TConfig>::VALUE;
    static const unsigned _BITS_PER_BLOCK   = RankDictionaryBitsPerBlock_<TValue, Levels<TSpec, TConfig> >::VALUE;
    static const unsigned _BITS_PER_WORD    = Min<RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::VALUE, _BITS_PER_BLOCK>::VALUE;
    static const unsigned _VALUES_PER_WORD  = _BITS_PER_WORD  / _BITS_PER_VALUE;
    static const unsigned _WORDS_PER_BLOCK  = _BITS_PER_BLOCK / _BITS_PER_WORD;
    static const unsigned _VALUES_PER_BLOCK = _VALUES_PER_WORD * _WORDS_PER_BLOCK;
    static const uint64_t _VALUES_PER_SUPERBLOCK = _VALUES_PER_BLOCK * 2; //(((1ull << ((TConfig::LEVELS == 3) ? 16 : 32)) - 1) / _VALUES_PER_BLOCK) * _VALUES_PER_BLOCK; // 2^16 - 1 (3lvl), 2^32 - 1 (2lvl)
    static const uint64_t _VALUES_PER_ULTRABLOCK = _VALUES_PER_BLOCK * 4; //(((1ull << 32) - 1) / _VALUES_PER_SUPERBLOCK) * _VALUES_PER_SUPERBLOCK; // 2^32 - 1

    typedef uint64_t  TWordType; // TODO: how to get underlying data type of RankDictionaryWordSize_?

    static TWordType _CHAR_BITMASKS[ValueSize<TValue>::VALUE]; // filter by character
    static TWordType _TRUNC_BITMASKS[_VALUES_PER_WORD]; // truncate the last values in a word that shell not be counted

    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreUltraBlocks>::Type  ultrablocks;
    typename Fibre<RankDictionary, FibreSuperBlocks>::Type  superblocks;
    typename Fibre<RankDictionary, FibreRanks>::Type        blocks;
    typename Size<RankDictionary>::Type                     _length;
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
            _CHAR_BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i, i + (1 << (_BITS_PER_VALUE-1)));
        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _TRUNC_BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, i+1, _BITS_PER_VALUE, 1, 1 << (_BITS_PER_VALUE-1)); // 1
    }

    template <typename TText>
    RankDictionary(TText const & text) :
        _length(0)
    {
        auto maxValue = (1 << _BITS_PER_VALUE) - 1;
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            _CHAR_BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue-i, i + (1 << (_BITS_PER_VALUE-1)));
        for (unsigned i = 0; i < _VALUES_PER_WORD; ++i)
            _TRUNC_BITMASKS[i] = _bitmaskWrapper<TWordType>(*this, _BITS_PER_WORD, i+1, _BITS_PER_VALUE, 1, 1 << (_BITS_PER_VALUE-1)); // 1

        createRankDictionary(*this, text);
    }
};

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_CHAR_BITMASKS[ValueSize<TValue>::VALUE];

template <typename TValue, typename TSpec, typename TConfig>
typename RankDictionary<TValue, Levels<TSpec, TConfig> >::TWordType RankDictionary<TValue, Levels<TSpec, TConfig> >::_TRUNC_BITMASKS[RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD];

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreRanks>::Type &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, FibreRanks)
{
    return dict.blocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreRanks>::Type const &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, FibreRanks)
{
    return dict.blocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperBlocks>::Type &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, FibreSuperBlocks)
{
    return dict.superblocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreSuperBlocks>::Type const &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, FibreSuperBlocks)
{
    return dict.superblocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreUltraBlocks>::Type &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, FibreUltraBlocks)
{
    return dict.ultrablocks;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreUltraBlocks>::Type const &
getFibre(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, FibreUltraBlocks)
{
    return dict.ultrablocks;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline bool empty(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict)
{
    return empty(getFibre(dict, FibreRanks())) && empty(getFibre(dict, FibreSuperBlocks())) && empty(getFibre(dict, FibreUltraBlocks()));
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline void clear(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict)
{
    clear(dict.blocks);
    clear(dict.superblocks);
    clear(dict.ultrablocks);
}

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
// Function _toPosInSuperBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInSuperBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK;
}

// ----------------------------------------------------------------------------
// Function _toPosInUltraBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInUltraBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_ULTRABLOCK;
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
// Function _toSuperBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toSuperBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK;
}

// ----------------------------------------------------------------------------
// Function _toUltraBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toUltraBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_ULTRABLOCK;
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
    return dict.blocks[blockPos].values[wordPos];
}

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TBlockPos blockPos, TWordPos wordPos)
{
    return dict.blocks[blockPos].values[wordPos];
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].values;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].values;
}

// ----------------------------------------------------------------------------
// Function _blockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].block;
}

// ----------------------------------------------------------------------------
// Function _superBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_superBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_superBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)];
}

// ----------------------------------------------------------------------------
// Function _ultraBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_ultraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    return dict.ultrablocks[_toUltraBlockPos(dict, pos)];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_ultraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos)
{
    return dict.ultrablocks[_toUltraBlockPos(dict, pos)];
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
    TSize endPos   = length(dict.blocks) * TRankDictionary::_VALUES_PER_BLOCK;

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

template <typename TSpec, typename TConfig, typename TPos>
inline void _clearBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    _blockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _clearSuperBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline void _clearSuperBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    clear(_superBlockAt(dict, pos));
}

template <typename TSpec, typename TConfig, typename TPos>
inline void _clearSuperBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    _superBlockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _clearUltraBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline void _clearUltraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    clear(_ultraBlockAt(dict, pos));
}

template <typename TSpec, typename TConfig, typename TPos>
inline void _clearUltraBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos pos)
{
    _ultraBlockAt(dict, pos) = 0u;
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

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getBlockRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & /* dict */, TBlock const & block, TPos /* pos */, TChar c, TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = block[ordValue(c)-1];
    smaller += _smaller; // Note: _smaller cannot be removed. order of evaluation is not defined!
    return block[ordValue(c)] - _smaller;
}

// TODO: könnte falsch sein!
template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline typename std::enable_if<TConfig::LEVELS == 1, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>::type
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TBlock const & block, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? block : pos - _toPosInBlock(dict, pos) - block;// _toPosInSuperBlock(dict, pos) - _toPosInBlock(dict, pos) - block;
}

template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline typename std::enable_if<TConfig::LEVELS >= 2, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>::type
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TBlock const & block, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? block : _toPosInSuperBlock(dict, pos) - _toPosInBlock(dict, pos) - block;
}

// ----------------------------------------------------------------------------
// Function _getSuperBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSuperBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getSuperBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TSuperBlock const & superblock, TPos /* pos */, TChar c)
{
    return superblock[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TSuperBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getSuperBlockRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & /* dict */, TSuperBlock const & superblock, TPos /* pos */, TChar c, TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = superblock[ordValue(c)-1];
    smaller += _smaller; // NOTE: _smaller cannot be removed. order of evaluation is not defined!
    return superblock[ordValue(c)] - _smaller;
}

// TODO: könnte falsch sein!
template <typename TSpec, typename TConfig, typename TSuperBlock, typename TPos>
inline typename std::enable_if<TConfig::LEVELS == 2, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>::type
_getSuperBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TSuperBlock const & superblock, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? superblock : /*_toPosInUltraBlock(dict, pos)*/ pos - _toPosInSuperBlock(dict, pos) - superblock;
}

template <typename TSpec, typename TConfig, typename TSuperBlock, typename TPos>
inline typename std::enable_if<TConfig::LEVELS == 3, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>::type
_getSuperBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TSuperBlock const & superblock, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? superblock : _toPosInUltraBlock(dict, pos) - _toPosInSuperBlock(dict, pos) - superblock;
}

// ----------------------------------------------------------------------------
// Function _getUltraBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TUltraBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getUltraBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TUltraBlock const & ultrablock, TPos /* pos */, TChar c)
{
    return ultrablock[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TUltraBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getUltraBlockRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & /* dict */, TUltraBlock const & ultrablock, TPos /* pos */, TChar c, TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = ultrablock[ordValue(c)-1];
    smaller += _smaller; // NOTE: _smaller cannot be removed. order of evaluation is not defined!
    return ultrablock[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TUltraBlock, typename TPos>
inline typename std::enable_if<TConfig::LEVELS == 3, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>::type
_getUltraBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TUltraBlock const & ultrablock, TPos pos, bool c)
{
    // If c == false then return the complementary rank.
    return c ? ultrablock : pos - _toPosInUltraBlock(dict, pos) - ultrablock;
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

    TWord mask = word ^ TRankDictionary::_CHAR_BITMASKS[ordValue(c)];

    // NOTE: actually it should be: mask & (mask >> 1) & (mask >> 2) & ... but this is shorter and equivalent
    for (TWord i = 1; i < TRankDictionary::_BITS_PER_VALUE; ++i)
        mask &= mask >> 1;

    return popCount(TRankDictionary::_TRUNC_BITMASKS[posInWord] & mask);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const &,
             TWord const & values,
             TPosInWord posInWord,
             TValue c)
{
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >    TRankDictionary;

    return popCount((TRankDictionary::_CHAR_BITMASKS[ordValue(c)] - values) & TRankDictionary::_TRUNC_BITMASKS[posInWord]);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TWord, typename TPosInWord, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const &,
             TWord const & values,
             TPosInWord posInWord,
             TValue c,
             TPos & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >            TRankDictionary;

    auto _smaller = popCount((TRankDictionary::_CHAR_BITMASKS[ordValue(c)-1] - values) & TRankDictionary::_TRUNC_BITMASKS[posInWord]);
    smaller += _smaller;
    return popCount((TRankDictionary::_CHAR_BITMASKS[ordValue(c)] - values) & TRankDictionary::_TRUNC_BITMASKS[posInWord]) - _smaller;
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

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TValues, typename TPosInBlock, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
_getValueRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & dict,
              TValues const & values,
              TPosInBlock posInBlock,
              TValue c,
              TSmaller & smaller)
{
    typedef RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > TRankDictionary;

    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    TSize valueRank = 0;

    // TODO: I guess we can use the original loop since I couldn't measure any performance difference and it is easier to read
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

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >::Type
_getValuesRanks(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & dict, TPos pos)
{
    typedef typename RankDictionaryBlock_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >::Type    TBlock;
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
template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos>
inline typename RankDictionaryBlock_<bool, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >::Type
_getValuesRanks(RankDictionary<bool, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & dict, TPos pos)
{
    return _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), true);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------
template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline typename std::enable_if<TConfig::LEVELS == 3, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>::type //Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TChar c)
{
    //std::cout << "== 3" << std::endl;
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const           TRankDictionary;
    //typedef typename Fibre<TRankDictionary, FibreSuperRanks>::Type          TFibreSuperRanks;
    //typedef typename Value<FibreSuperRanks>::Type                           TFibreRank;
    //typedef typename Fibre<TRankDictionary, FibreSuperRanks>::Type          TFibreRank;
    //typedef typename Value<TFibreRanks>::Type                               TRankEntry;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize posInBlock    = _toPosInBlock(dict, pos);

    auto const & ultraBlock = _ultraBlockAt(dict, pos);
    auto const & superBlock = _superBlockAt(dict, pos);
    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c))
         + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline typename std::enable_if<TConfig::LEVELS == 2, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>::type
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TChar c)
{
    //std::cout << "== 2" << std::endl;
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const           TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize posInBlock    = _toPosInBlock(dict, pos);

    auto const & superBlock = _superBlockAt(dict, pos);
    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline typename std::enable_if<TConfig::LEVELS == 1, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>::type
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TChar c)
{
    //std::cout << "== 1" << std::endl;
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const           TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                            TSize;

    TSize posInBlock    = _toPosInBlock(dict, pos);

    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

// TODO: work around for unidrectional index without smaller
template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos pos, TChar c, TSmaller /**/)
{
    return getRank(dict, pos, c);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & dict, TPos pos, TChar c)
{
    TPos smaller;
    return getRank(dict, pos, static_cast<TValue>(c), smaller);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 3> > > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 3> > > const & dict, TPos pos, TValue c, TPos & smaller)
{
    //std::cout << "==3" << std::endl;
    TSize posInBlock = _toPosInBlock(dict, pos);

    auto const & ultraBlock = _ultraBlockAt(dict, pos);
    auto const & superBlock = _superBlockAt(dict, pos);
    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    smaller = 0;
    if (ordValue(c) > 0)
        return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c), smaller)
             + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c), smaller)
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);

    // c == Dna('A')
    return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c))
         + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 2> > > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 2> > > const & dict, TPos pos, TValue c, TPos & smaller)
{
    //std::cout << "==2" << std::endl;
    TSize posInBlock = _toPosInBlock(dict, pos);

    auto const & superBlock = _superBlockAt(dict, pos);
    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    smaller = 0;
    if (ordValue(c) > 0)
        return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c), smaller)
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);

    // c == Dna('A')
    return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 1> > > const>::Type
getRank(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, 1> > > const & dict, TPos pos, TValue c, TPos & smaller)
{
    //std::cout << "==1" << std::endl;
    TSize posInBlock = _toPosInBlock(dict, pos);

    auto const & entry = dict.blocks[_toBlockPos(dict, pos)];

    smaller = 0;
    if (ordValue(c) > 0)
        return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);

    // c == Dna('A')
    return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

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

    TSize blockPos      = _toBlockPos(dict, pos);
    TSize posInBlock    = _toPosInBlock(dict, pos);
    TSize wordPos       = _toWordPos(dict, posInBlock);
    TSize posInWord     = _toPosInWord(dict, posInBlock);

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

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > >::Type
getValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > & dict, TPos pos)
{
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    auto SIZE = Size<typename RankDictionaryValues_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >::TValues>::VALUE;
    unsigned shift = (SIZE - posInWord - 1) * (BitsPerValue<TValue>::VALUE + 1);
    auto value = _valuesAt(dict, blockPos, wordPos);

    return (value.i >> shift) & value.BIT_MASK2;
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const>::Type
getValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > const & dict, TPos pos)
{
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    auto SIZE = Size<typename RankDictionaryValues_<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > >::TValues>::VALUE;
    unsigned shift = (SIZE - posInWord - 1) * (BitsPerValue<TValue>::VALUE + 1);
    auto value = _valuesAt(dict, blockPos, wordPos);

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

    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    assignValue(_valuesAt(dict, blockPos, wordPos), posInWord, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, typename TPos, typename TChar>
inline void setValue(RankDictionary<TValue, Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS> > > & dict, TPos pos, TChar c)
{
    TSize blockPos   = _toBlockPos(dict, pos);
    TSize posInBlock = _toPosInBlock(dict, pos);
    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    assignValue2(_valuesAt(dict, blockPos, wordPos), posInWord, static_cast<TValue>(c));
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Better not to have appendValue() - it is not efficient - and thus neither length().
// TODO: can we remove this?
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
    typedef typename Iterator<TFibreRanks, Standard>::Type          TBlockIter;

    if (empty(dict)) return;

    // Insures the first superblock ranks start from zero.
    if (TConfig::LEVELS > 2)
        _clearUltraBlockAt(dict, 0u);
    if (TConfig::LEVELS > 1)
        _clearSuperBlockAt(dict, 0u);
    _clearBlockAt(dict, 0u);

    // Clear the uninitialized values.
    _padValues(dict);

    static const unsigned blocks_per_ultrablock = TRankDictionary::_VALUES_PER_ULTRABLOCK / TRankDictionary::_VALUES_PER_BLOCK;
    static const unsigned blocks_per_superblock = TRankDictionary::_VALUES_PER_SUPERBLOCK / TRankDictionary::_VALUES_PER_BLOCK;

    TBlockIter blocksBegin = begin(dict.blocks, Standard());
    TBlockIter blocksEnd = end(dict.blocks, Standard());

    for (TBlockIter blocksIt = blocksBegin; blocksIt != blocksEnd; ++blocksIt)
    {
        TSize blockPos = blocksIt - blocksBegin;
        TSize curr = _toPos(dict, blockPos);
        TSize next = _toPos(dict, blockPos + 1);

        if (blockPos % blocks_per_superblock == 0 && TConfig::LEVELS > 1)
        {
            if (blockPos > 0)
            {
                TSize prevSB = _toPos(dict, blockPos - blocks_per_superblock);
                _superBlockAt(dict, curr) = _superBlockAt(dict, prevSB) + _blockAt(dict, curr);
                _clearBlockAt(dict, curr);
                if (blockPos % blocks_per_ultrablock == 0 && TConfig::LEVELS > 2)
                {
                    TSize prevUB = _toPos(dict, blockPos - blocks_per_ultrablock);
                    _ultraBlockAt(dict, curr) = _ultraBlockAt(dict, prevUB) + _superBlockAt(dict, curr);
                    _clearSuperBlockAt(dict, curr);
                }
            }
        }

        if (blocksIt != blocksEnd - 1)
            _blockAt(dict, next) = _blockAt(dict, curr) + _getValuesRanks(dict, next - 1);
    }
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
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDict_;

    bool res = true;
    res &= reserve(dict.blocks, (newCapacity + TRankDict_::_VALUES_PER_BLOCK - 1) / TRankDict_::_VALUES_PER_BLOCK, tag);
    if (TConfig::LEVELS > 1)
        res &= resize(dict.superblocks, (newCapacity + TRankDict_::_VALUES_PER_SUPERBLOCK - 1) / TRankDict_::_VALUES_PER_SUPERBLOCK, tag);
    if (TConfig::LEVELS > 2)
        res &= resize(dict.ultrablocks, (newCapacity + TRankDict_::_VALUES_PER_ULTRABLOCK - 1) / TRankDict_::_VALUES_PER_ULTRABLOCK, tag);
    return res;
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
resize(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize newLength, Tag<TExpand> const tag)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDict_;

    bool result = true;
    dict._length = newLength;
    result &= resize(dict.blocks, (newLength + TRankDict_::_VALUES_PER_BLOCK - 1) / TRankDict_::_VALUES_PER_BLOCK, tag);
    if (TConfig::LEVELS > 1)
        result &= resize(dict.superblocks, (newLength + TRankDict_::_VALUES_PER_SUPERBLOCK - 1) / TRankDict_::_VALUES_PER_SUPERBLOCK, tag);
    if (TConfig::LEVELS > 2)
        result &= resize(dict.ultrablocks, (newLength + TRankDict_::_VALUES_PER_ULTRABLOCK - 1) / TRankDict_::_VALUES_PER_ULTRABLOCK, tag);
    return result;
}

template <typename TValue, typename TSpec, typename TConfig>
inline bool save(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, const char * fileName, int openMode)
{
    String<char> name;
    bool result = true;
    name = fileName;    append(name, ".bl");    result &= save(getFibre(dict, FibreRanks()), toCString(name), openMode);
    if (TConfig::LEVELS > 1)
        name = fileName;    append(name, ".sbl");   result &= save(getFibre(dict, FibreSuperBlocks()), toCString(name), openMode);
    if (TConfig::LEVELS > 2)
        name = fileName;    append(name, ".ubl");   result &= save(getFibre(dict, FibreUltraBlocks()), toCString(name), openMode);
    return result;
}

template <typename TValue, typename TSpec, typename TConfig>
inline bool open(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, const char * fileName, int openMode)
{
    String<char> name;
    bool result = true;
    name = fileName;    append(name, ".bl");    result &= open(getFibre(dict, FibreRanks()), toCString(name), openMode);
    if (TConfig::LEVELS > 1)
        name = fileName;    append(name, ".sbl");   result &= open(getFibre(dict, FibreSuperBlocks()), toCString(name), openMode);
    if (TConfig::LEVELS > 2)
        name = fileName;    append(name, ".ubl");   result &= open(getFibre(dict, FibreUltraBlocks()), toCString(name), openMode);
    return result;
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
