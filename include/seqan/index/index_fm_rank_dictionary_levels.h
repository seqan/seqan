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
//       & Christopher Pockrandt <christopher.pockrandt@fu-berlin.de>
// ==========================================================================

#ifndef INDEX_FM_RANK_DICTIONARY_LEVELS_H_
#define INDEX_FM_RANK_DICTIONARY_LEVELS_H_

#include <algorithm>

#define TPREFIXLEVELS Levels<TSpec, LevelsPrefixRDConfig<TSize, TFibre, LEVELS, WPB> >

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
// Metafunction RankDictionaryBitString_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBitString_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlockSize_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
struct RankDictionaryBlockSize_;

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryTruncBitmaskSize_
// ----------------------------------------------------------------------------

template <unsigned _VALUES_PER_WORD, typename TValue, typename TSpec>
struct RankDictionaryTruncBitmaskSize_;

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

/*!
 * @class LevelsRDConfig
 * @headerfile <seqan/index.h>
 *
 * @brief LevelsRDConfig allows configuring a @link Levels @endlink RankDictionary.
 *
 * @signature template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS = 1, unsigned WORDS_PER_BLOCK = 0>
 *            struct LevelsRDConfig<TSize, TFibre, LEVELS, WORDS_PER_BLOCK>;
 *
 * @tparam TSize           A data type that can store the length of the input text. Default: <tt>size_t</tt>
 * @tparam TFibre          A tag for specialization purposes of the underlying strings. Default: <tt>Alloc<></tt>
 * @tparam LEVELS          The number of levels (1, 2, or 3). The more levels, the lower the space consumption but possibly slight performance decreases. Default: <tt>1</tt>
 * @tparam WORDS_PER_BLOCK The number of popcount operations per rank query. A lower number implies more space for faster runtime. 0 is a shortcut for the size of the alphabet of the RankDictionary. Default: <tt>0</tt>
 */

template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS = 1, unsigned WORDS_PER_BLOCK = 0>
struct LevelsRDConfig : RDConfig<TSize, TFibre, LEVELS, WORDS_PER_BLOCK> {};

/*!
 * @class LevelsPrefixRDConfig
 * @headerfile <seqan/index.h>
 *
 * @brief LevelsPrefixRDConfig allows configuring a @link Levels @endlink RankDictionary that is recommended for fast searching in bidirectional indices.
 *
 * @signature template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS = 1, unsigned WORDS_PER_BLOCK = 0>
 *            struct LevelsRDConfig<TSize, TFibre, LEVELS, WORDS_PER_BLOCK>;
 *
 * @tparam TSize           A data type that can store the length of the input text. Default: <tt>size_t</tt>
 * @tparam TFibre          A tag for specialization purposes of the underlying strings. Default: <tt>Alloc<></tt>
 * @tparam LEVELS          The number of levels (1, 2, or 3). The more levels, the lower the space consumption but possibly slight performance decreases. Default: <tt>1</tt>
 * @tparam WORDS_PER_BLOCK The number of popcount operations per rank query. A lower number implies more space for faster runtime. 0 is a shortcut for the size of the alphabet of the RankDictionary. Default: <tt>0</tt>
 */
template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS = 1, unsigned WORDS_PER_BLOCK = 0>
struct LevelsPrefixRDConfig : RDConfig<TSize, TFibre, LEVELS, WORDS_PER_BLOCK> {};

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
struct RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> > : BitsPerValue<uint64_t>
{
    typedef uint64_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBitsPerBlock_
// ----------------------------------------------------------------------------
// The number of bits per block equals the number of bits of the block summary.

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBitsPerBlock_<TValue, Levels<TSpec, TConfig> > :
    BitsPerValue<typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type> {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlockSize_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename ValueSize<TValue>::Type TType;
    static constexpr TType VALUE = ValueSize<TValue>::VALUE;
};

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB>
struct RankDictionaryBlockSize_<TValue, TPREFIXLEVELS>
{
    typedef typename ValueSize<TValue>::Type TType;
    static constexpr TType VALUE = ValueSize<TValue>::VALUE - 1;
};

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB>
struct RankDictionaryBlockSize_<bool, TPREFIXLEVELS>
{
    typedef typename ValueSize<bool>::Type TType;
    // actually one bitmask would be sufficient. But optimizing a prefix rank dictionary for bools is not worth it ;-)
    static constexpr TType VALUE = ValueSize<bool>::VALUE;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryTruncBitmaskSize_
// ----------------------------------------------------------------------------

template <unsigned _VALUES_PER_WORD, typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >
{
    static constexpr unsigned VALUE = _VALUES_PER_WORD;
};

template <unsigned _VALUES_PER_WORD, typename TValue, typename TSpec,
          typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB>
struct RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, TPREFIXLEVELS>
{
    static constexpr unsigned VALUE = (_VALUES_PER_WORD + 1)/2 + 1;
};

template <unsigned _VALUES_PER_WORD, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB>
struct RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, bool, TPREFIXLEVELS>
{
    static constexpr unsigned VALUE = _VALUES_PER_WORD;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValues_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >                             TRankDictionary_;

    typedef Tuple<TValue, TRankDictionary_::_VALUES_PER_WORD, BitPacked<> >             TValues;
    typedef typename TValues::TBitVector                                                TWord;
    typedef Tuple<TValues, TRankDictionary_::_WORDS_PER_BLOCK>                          Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryValueWithBits_
// ----------------------------------------------------------------------------

template <unsigned BITS>
struct RankDictionaryValueWithBits_ {};

template <>
struct RankDictionaryValueWithBits_<64>
{
    typedef uint64_t Type;
};

template <>
struct RankDictionaryValueWithBits_<32>
{
    typedef uint32_t Type;
};

template <>
struct RankDictionaryValueWithBits_<16>
{
    typedef uint16_t Type;
};

template <>
struct RankDictionaryValueWithBits_<8>
{
    typedef uint8_t Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlockType_
// ----------------------------------------------------------------------------

template <typename TSize, unsigned LEVELS, unsigned LEVEL>
struct RankDictionaryBlockType_
{
    static constexpr unsigned shift = Max<LEVEL, LEVELS>::VALUE - Min<LEVEL, LEVELS>::VALUE;
    typedef typename RankDictionaryValueWithBits_<BitsPerValue<TSize>::VALUE/(1 << shift)>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type   TSize_;
    typedef Tuple<TSize_, RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE>     Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionaryBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionarySuperBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 2>::Type   TSize_;
    typedef Tuple<TSize_, RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE>     Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionarySuperBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 2>::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryUltraBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 3>::Type   TSize_;
    typedef Tuple<TSize_, RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE>     Type;
};

template <typename TSpec, typename TConfig>
struct RankDictionaryUltraBlock_<bool, Levels<TSpec, TConfig> >
{
    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 3>::Type   Type;
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
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >                             TRankDictionary_;
    typedef typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type   TSuperBlocks_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type                     TFibreSpec_;

    typedef String<TSuperBlocks_, TFibreSpec_>                                          Type;
};

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Levels<TSpec, TConfig> >, FibreUltraBlocks>
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >                             TRankDictionary_;
    typedef typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type   TUltraBlocks_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type                     TFibreSpec_;

    typedef String<TUltraBlocks_, TFibreSpec_>                                          Type;
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

template <typename TWord>
constexpr TWord
_bitmaskRecursive(unsigned const bitsUsed, unsigned const blocksLeft, unsigned const blocksProcessed,
                  unsigned const blocksize, TWord const value1, TWord const value2, TWord const result)
{
    return (blocksLeft == 0)
        ? result
        : _bitmaskRecursive<TWord>(bitsUsed,
                                   blocksLeft - 1,
                                   blocksProcessed + 1,
                                   blocksize,
                                   value1,
                                   value2,
                                   result | (((blocksProcessed % 2) ? value1 : value2)
                                           << (bitsUsed - blocksProcessed * blocksize)));
}

template <typename TWord>
constexpr TWord
_bitmask(unsigned const bitsTotal, unsigned const blocks, unsigned const blocksize, TWord const value1,
         TWord const value2)
{
    return _bitmaskRecursive(bitsTotal - (bitsTotal % blocksize), blocks, 1, blocksize, value1, value2, (TWord) 0);
    // TODO(cpockrandt): Replace recusive version by iterative, easier to understand code once gcc-4.9 support is dropped
    /*unsigned bitsUsed = bitsTotal - (bitsTotal % blocksize);
    TWord bm = 0u;
    for (unsigned i = 1; i <= blocks; ++i)
        bm |= (((i % 2) ? value1 : value2) << (bitsUsed - i * blocksize));
    return bm;*/
}

template <typename T>
struct isLevelsPrefixRD
{
    static const bool VALUE = false;
};

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB>
struct isLevelsPrefixRD<TPREFIXLEVELS>
{
    static const bool VALUE = true;
};

// ----------------------------------------------------------------------------
// Struct Levels RankDictionary
// ----------------------------------------------------------------------------
/*!
 * @class Levels
 * @extends RankDictionary
 * @headerfile <seqan/index.h>
 *
 * @brief Levels is a @link RankDictionary @endlink consisting of up to three levels.
 *
 * @signature template <typename TValue, typename TSpec, typename TConfig>
 *            struct RankDictionary<TValue, Levels<TSpec> >;
 *
 * @tparam TValue The alphabet type of the wavelet tree.
 * @tparam TSpec  A tag for specialization purposes. Default: <tt>void</tt>
 *
 * This @link RankDictionary @endlink consists of up to three levels of rank
 * infromation, in which one stores information of blocks (1st level) and the other
 * information until a specified block. Combining those two information
 * leads to constant rank dictionary look ups. To reduce space one can increase
 * the number of levels up to three.
 */

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, Levels<TSpec, TConfig> >
{
    // ------------------------------------------------------------------------
    // Constants
    // ------------------------------------------------------------------------

    static constexpr unsigned _BITS_PER_VALUE   = BitsPerValue<TValue>::VALUE;
    static constexpr unsigned _BITS_PER_BLOCK   = (TConfig::WORDS_PER_BLOCK == 0 ? BitsPerValue<typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type>::VALUE : RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::VALUE * TConfig::WORDS_PER_BLOCK);
    static constexpr unsigned _BITS_PER_WORD    = Min<RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::VALUE, _BITS_PER_BLOCK>::VALUE;
    static constexpr unsigned _VALUES_PER_WORD  = _BITS_PER_WORD / _BITS_PER_VALUE;
    static constexpr unsigned _WORDS_PER_BLOCK  = _BITS_PER_BLOCK / _BITS_PER_WORD;
    static constexpr unsigned _VALUES_PER_BLOCK = _VALUES_PER_WORD * _WORDS_PER_BLOCK;
    static constexpr uint64_t _VALUES_PER_SUPERBLOCK = Max<1, (((1ull << (BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 2>::Type>::VALUE/2)) - 1) / _VALUES_PER_BLOCK) * _VALUES_PER_BLOCK>::VALUE; // 2^x - 1 values
    static constexpr uint64_t _VALUES_PER_ULTRABLOCK = Max<1, (((1ull << (BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 3>::Type>::VALUE/2)) - 1) / Max<_VALUES_PER_SUPERBLOCK, 1>::VALUE) * _VALUES_PER_SUPERBLOCK>::VALUE; // MAX: workaround for intel & clang: division is not constexpr since divisor could be 0

    typedef typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type   Type;

    static_assert(
        BitsPerValue<typename RankDictionaryBlockType_<typename TConfig::Size, TConfig::LEVELS, 1>::Type>::VALUE >= LogN<_VALUES_PER_BLOCK + 1, 2>::VALUE,
        "The data type of the lowest level has to be larger (or the number of words per block smaller).");

    typedef typename RankDictionaryWordSize_<TValue, Levels<TSpec, TConfig> >::Type TWordType;

    static TWordType _SELECT_BITMASK; // select only every 2nd element

    // filter by character
    static TWordType _CHAR_BITMASKS[RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE];

    // truncate the last values in a word that shall not be counted
    static TWordType
    _TRUNC_BITMASKS[RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >::VALUE];

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

    void _populateBitmasks()
    {
        // TODO(cpockrandt): make all bitmasks a const-expr
        constexpr unsigned maxValue = (1 << _BITS_PER_VALUE) - 1;
        constexpr unsigned padding = _BITS_PER_WORD % _BITS_PER_VALUE;
        // max() workaround for windows C4293 warning
        TWordType paddingWord = (padding > 0) ? (1ull << (_BITS_PER_WORD - std::max(padding, 1u))) : 0ull;
        if (isLevelsPrefixRD<Levels<TSpec, TConfig> >::VALUE && !std::is_same<TValue, bool>::value)
        {
            if (padding > 0)
                _SELECT_BITMASK = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, maxValue, 0u);
            else
                _SELECT_BITMASK = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, 0u, maxValue);
            for (unsigned i = 0; i < RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE; ++i)
            {
                if (padding > 0)
                    _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, i, 1)
                                        | paddingWord;
                else
                    _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, _VALUES_PER_WORD, _BITS_PER_VALUE, 1, i);
            }
            _TRUNC_BITMASKS[0] = 0u;
            for (unsigned i = 0;
                 i < RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >::VALUE - 1;
                 ++i)
            {
                if (padding > 0)
                    _TRUNC_BITMASKS[i+1] = _bitmask<TWordType>(_BITS_PER_WORD,2*i + 1, _BITS_PER_VALUE, 0, 1)
                                           | paddingWord;
                else
                    _TRUNC_BITMASKS[i+1] = _bitmask<TWordType>(_BITS_PER_WORD, 2*i + 1, _BITS_PER_VALUE, 1, 0);
            }
        }
        else
        {
            for (unsigned i = 0; i < RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE; ++i)
                _CHAR_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD,
                                                        _VALUES_PER_WORD,
                                                        _BITS_PER_VALUE,
                                                        maxValue-i,
                                                        maxValue-i);
            for (unsigned i = 0;
                 i < RankDictionaryTruncBitmaskSize_<_VALUES_PER_WORD, TValue, Levels<TSpec, TConfig> >::VALUE;
                 ++i)
                _TRUNC_BITMASKS[i] = _bitmask<TWordType>(_BITS_PER_WORD, i+1, _BITS_PER_VALUE, 1, 1);
        }
    }

    RankDictionary() :
        _length(0)
    {
        _populateBitmasks();
    }

    template <typename TText>
    RankDictionary(TText const & text) :
        _length(0)
    {
        _populateBitmasks();

        createRankDictionary(*this, text);
    }
};

#define TRD RankDictionary<TValue, Levels<TSpec, TConfig> >
template <typename TValue, typename TSpec, typename TConfig>
typename TRD::TWordType TRD::_SELECT_BITMASK;

template <typename TValue, typename TSpec, typename TConfig>
typename TRD::TWordType TRD::_CHAR_BITMASKS[RankDictionaryBlockSize_<TValue, Levels<TSpec, TConfig> >::VALUE];

template <typename TValue, typename TSpec, typename TConfig>
typename TRD::TWordType TRD::_TRUNC_BITMASKS[RankDictionaryTruncBitmaskSize_<TRD::_VALUES_PER_WORD,
                                                                             TValue,
                                                                             Levels<TSpec, TConfig> >::VALUE];
#undef TRD

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
    return empty(getFibre(dict, FibreRanks()))
        && empty(getFibre(dict, FibreSuperBlocks()))
        && empty(getFibre(dict, FibreUltraBlocks()));
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
_toPosInWord(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const posInBlock)
{
    return posInBlock % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD;
}

// ----------------------------------------------------------------------------
// Function _toWordPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toWordPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const posInBlock)
{
    return posInBlock / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD;
}

// ----------------------------------------------------------------------------
// Function _toPosInBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _toPosInSuperBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInSuperBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK;
}

// ----------------------------------------------------------------------------
// Function _toPosInUltraBlock()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPosInUltraBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos % RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_ULTRABLOCK;
}

// ----------------------------------------------------------------------------
// Function _toBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _toSuperBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toSuperBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_SUPERBLOCK;
}

// ----------------------------------------------------------------------------
// Function _toUltraBlockPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toUltraBlockPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TPos const pos)
{
    return pos / RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_ULTRABLOCK;
}

// ----------------------------------------------------------------------------
// Function _toPos()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
_toPos(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */, TBlockPos const blockPos)
{
    return blockPos * RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_BLOCK;
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TBlockPos const blockPos, TWordPos const wordPos)
{
    return dict.blocks[blockPos].values[wordPos];
}

template <typename TValue, typename TSpec, typename TConfig, typename TBlockPos, typename TWordPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::TValues const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict,
          TBlockPos const blockPos,
          TWordPos const wordPos)
{
    return dict.blocks[blockPos].values[wordPos];
}

// ----------------------------------------------------------------------------
// Function _valuesAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].values;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryValues_<TValue, Levels<TSpec, TConfig> >::Type const &
_valuesAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].values;
}

// ----------------------------------------------------------------------------
// Function _blockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].block;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_blockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    return dict.blocks[_toBlockPos(dict, pos)].block;
}

// ----------------------------------------------------------------------------
// Function _superBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_superBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionarySuperBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_superBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    return dict.superblocks[_toSuperBlockPos(dict, pos)];
}

// ----------------------------------------------------------------------------
// Function _ultraBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type &
_ultraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    return dict.ultrablocks[_toUltraBlockPos(dict, pos)];
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryUltraBlock_<TValue, Levels<TSpec, TConfig> >::Type const &
_ultraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
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
inline void _clearBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    clear(_blockAt(dict, pos));
}

template <typename TSpec, typename TConfig, typename TPos>
inline void _clearBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    _blockAt(dict, pos) = 0u;
}

// ----------------------------------------------------------------------------
// Function _clearSuperBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS >= 2, void>
_clearSuperBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    clear(_superBlockAt(dict, pos));
}

template <typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS >= 2, void>
_clearSuperBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    _superBlockAt(dict, pos) = 0u;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS == 1, void>
_clearSuperBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & /* dict */, TPos const /* pos */) {}

// ----------------------------------------------------------------------------
// Function _clearUltraBlockAt()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS == 3, void>
_clearUltraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    clear(_ultraBlockAt(dict, pos));
}

template <typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS == 3, void>
_clearUltraBlockAt(RankDictionary<bool, Levels<TSpec, TConfig> > & dict, TPos const pos)
{
    _ultraBlockAt(dict, pos) = 0u;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS <= 2, void>
_clearUltraBlockAt(RankDictionary<TValue, Levels<TSpec, TConfig> > & /* dict */, TPos const /* pos */) {}

// ----------------------------------------------------------------------------
// Function _getBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */,
              TBlock const & block,
              TPos const /* pos */,
              TChar const c)
{
    return block[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type
_getBlockRank(RankDictionary<TValue, TPREFIXLEVELS> const & /* dict */,
              TBlock const & block,
              TPos const /* pos */,
              TChar const c,
              TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = block[ordValue(c)-1];
    smaller += _smaller; // NOTE: _smaller cannot be removed. order of evaluation is not defined!
    return block[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline
typename std::enable_if_t<TConfig::LEVELS == 1, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict,
              TBlock const & block,
              TPos const pos,
              bool const c)
{
    // If c == false then return the complementary rank.
    return c ? block : pos - _toPosInBlock(dict, pos) - block;
}

template <typename TSpec, typename TConfig, typename TBlock, typename TPos>
inline
typename std::enable_if_t<TConfig::LEVELS >= 2, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>
_getBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict,
              TBlock const & block,
              TPos const pos,
              bool const c)
{
    // If c == false then return the complementary rank.
    return c ? block : _toPosInSuperBlock(dict, pos) - _toPosInBlock(dict, pos) - block;
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TBlock, typename TPos, typename TSmaller>
inline typename Size<RankDictionary<bool, TPREFIXLEVELS> const>::Type
_getBlockRank(RankDictionary<bool, TPREFIXLEVELS> const & dict,
              TBlock const & block,
              TPos const pos,
              bool const /* c */,
              TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    smaller += _getBlockRank(dict, block, pos, false);
    return _getBlockRank(dict, block, pos, true);
}

// ----------------------------------------------------------------------------
// Function _getSuperBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSuperBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getSuperBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */,
                   TSuperBlock const & superblock,
                   TPos const /* pos */,
                   TChar const c)
{
    return superblock[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TSuperBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type
_getSuperBlockRank(RankDictionary<TValue, TPREFIXLEVELS> const & /* dict */,
                   TSuperBlock const & superblock,
                   TPos const /* pos */,
                   TChar const c,
                   TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = superblock[ordValue(c)-1];
    smaller += _smaller; // NOTE: _smaller cannot be removed. order of evaluation is not defined!
    return superblock[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TSuperBlock, typename TPos>
inline
typename std::enable_if_t<TConfig::LEVELS == 2, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>
_getSuperBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict,
                   TSuperBlock const & superblock,
                   TPos const pos,
                   bool const c)
{
    // If c == false then return the complementary rank.
    return c ? superblock : pos - _toPosInSuperBlock(dict, pos) - superblock;
}

template <typename TSpec, typename TConfig, typename TSuperBlock, typename TPos>
inline
typename std::enable_if_t<TConfig::LEVELS == 3, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>
_getSuperBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict,
                   TSuperBlock const & superblock,
                   TPos const pos,
                   bool const c)
{
    // If c == false then return the complementary rank.
    return c ? superblock : _toPosInUltraBlock(dict, pos) - _toPosInSuperBlock(dict, pos) - superblock;
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TSuperBlock, typename TPos, typename TSmaller>
inline typename Size<RankDictionary<bool, TPREFIXLEVELS> const>::Type
_getSuperBlockRank(RankDictionary<bool, TPREFIXLEVELS> const & dict,
                   TSuperBlock const & superblock,
                   TPos const pos,
                   bool const /*c*/,
                   TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    smaller += _getSuperBlockRank(dict, superblock, pos, false);
    return _getSuperBlockRank(dict, superblock, pos, true);
}

// ----------------------------------------------------------------------------
// Function _getUltraBlockRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TUltraBlock, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getUltraBlockRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */,
                   TUltraBlock const & ultrablock,
                   TPos const /* pos */,
                   TChar const c)
{
    return ultrablock[ordValue(c)];
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TUltraBlock, typename TPos, typename TChar, typename TSmaller>
inline typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type
_getUltraBlockRank(RankDictionary<TValue, TPREFIXLEVELS> const & /* dict */,
                   TUltraBlock const & ultrablock,
                   TPos const /* pos */,
                   TChar const c,
                   TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TSmaller _smaller = ultrablock[ordValue(c)-1];
    smaller += _smaller; // NOTE: _smaller cannot be removed. order of evaluation is not defined!
    return ultrablock[ordValue(c)] - _smaller;
}

template <typename TSpec, typename TConfig, typename TUltraBlock, typename TPos>
inline typename std::enable_if_t<TConfig::LEVELS == 3, typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type>
_getUltraBlockRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict,
                   TUltraBlock const & ultrablock,
                   TPos const pos,
                   bool const c)
{
    // If c == false then return the complementary rank.
    return c ? ultrablock : pos - _toPosInUltraBlock(dict, pos) - ultrablock;
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TUltraBlock, typename TPos, typename TSmaller>
inline typename Size<RankDictionary<bool, TPREFIXLEVELS> const>::Type
_getUltraBlockRank(RankDictionary<bool, TPREFIXLEVELS> const & dict,
                   TUltraBlock const & ultrablock,
                   TPos const pos,
                   bool const /*c*/,
                   TSmaller & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    smaller += _getUltraBlockRank(dict, ultrablock, pos, false);
    return _getUltraBlockRank(dict, ultrablock, pos, true);
}

// ----------------------------------------------------------------------------
// Function _getWordRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TWord const & word, TValue const c)
{
    return _getWordRank(dict, word, RankDictionary<TValue, Levels<TSpec, TConfig> >::_VALUES_PER_WORD - 1, c);
}

template <typename TValue, typename TSpec, typename TConfig, typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getWordRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & /* dict */,
             TWord const & word,
             TPosInWord const posInWord,
             TValue const c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDictionary;

    TWord mask = word ^ TRankDictionary::_CHAR_BITMASKS[ordValue(c)];
    // NOTE: actually it should be: mask & (mask >> 1) & (mask >> 2) & ... but this is shorter and equivalent
    for (TWord i = 1; i < TRankDictionary::_BITS_PER_VALUE; ++i)
        mask &= mask >> 1;
    return popCount(TRankDictionary::_TRUNC_BITMASKS[posInWord] & mask);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TWord, typename TPosInWord>
inline typename std::enable_if_t<RankDictionary<TValue, TPREFIXLEVELS>::_BITS_PER_WORD % RankDictionary<TValue, TPREFIXLEVELS>::_BITS_PER_VALUE == 0, typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type>
_getWordRank(RankDictionary<TValue, TPREFIXLEVELS> const &,
             TWord const & word,
             TPosInWord const posInWord,
             TValue const c)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS> TRD;

    TWord const erg1 = (TRD::_CHAR_BITMASKS[ordValue(c)] - (word & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord + 1)/2];
    TWord const erg2 = (TRD::_CHAR_BITMASKS[ordValue(c)] - ((word >> TRD::_BITS_PER_VALUE) & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord/2) + 1];

    return popCount(erg1 | (erg2 << 1));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TWord, typename TPosInWord>
inline
typename std::enable_if_t<RankDictionary<TValue, TPREFIXLEVELS>::_BITS_PER_WORD % RankDictionary<TValue, TPREFIXLEVELS>::_BITS_PER_VALUE != 0, typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type>
_getWordRank(RankDictionary<TValue, TPREFIXLEVELS> const &,
             TWord const & word,
             TPosInWord const posInWord,
             TValue const c)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS>    TRD;

    TWord erg1 = (TRD::_CHAR_BITMASKS[ordValue(c)]
                  - (word & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord/2) + 1];
    TWord erg2 = (TRD::_CHAR_BITMASKS[ordValue(c)]
                  - ((word << TRD::_BITS_PER_VALUE) & TRD::_SELECT_BITMASK)) & TRD::_TRUNC_BITMASKS[(posInWord + 1)/2];

    return popCount(erg1 | (erg2 >> 1));
}

// TODO(cpockrandt): rename functions with smaller value (cause they return a different value than their counterparts)
template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TWord, typename TPosInWord, typename TPos>
inline typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type
_getWordRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
             TWord const & word,
             TPosInWord const posInWord,
             TValue const c,
             TPos & smaller)
{
    // can only be called if ordValue(c) > 0. smaller has to be initialized by the caller!
    TPos _smaller = _getWordRank(dict, word, posInWord, static_cast<TValue>(ordValue(c) - 1));
    smaller += _smaller;
    return _getWordRank(dict, word, posInWord, c) - _smaller;
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TWord, typename TPosInWord>
inline typename Size<RankDictionary<bool, TPREFIXLEVELS> const>::Type
_getWordRank(RankDictionary<bool, TPREFIXLEVELS> const &,
             TWord const & word,
             TPosInWord const posInWord,
             bool const c)
{
    typedef RankDictionary<bool, TPREFIXLEVELS>    TRD;

    return popCount(TRD::_TRUNC_BITMASKS[posInWord] & (word ^ TRD::_CHAR_BITMASKS[ordValue(c)]));
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TWord, typename TPosInWord, typename TPos>
inline typename Size<RankDictionary<bool, TPREFIXLEVELS> const>::Type
_getWordRank(RankDictionary<bool, TPREFIXLEVELS> const & dict,
             TWord const & word,
             TPosInWord const posInWord,
             bool const /* c */,
             TPos & smaller)
{
    TPos _smaller = _getWordRank(dict, word, posInWord, false);
    smaller += _smaller;
    return posInWord - _smaller + 1;
}

// ----------------------------------------------------------------------------
// Function _getValueRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TValues, typename TPosInBlock>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type
_getValueRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict,
              TValues const & values,
              TPosInBlock const posInBlock,
              TValue const c)
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
        if (wordPrevPos < wordPos)
            valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    valueRank += _getWordRank(dict, values[wordPos].i, posInWord, c);

    return valueRank;
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TValues, typename TPosInBlock, typename TSmaller>
inline typename Size<RankDictionary<TValue, TPREFIXLEVELS> const>::Type
_getValueRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
              TValues const & values,
              TPosInBlock const posInBlock,
              TValue const c,
              TSmaller & smaller)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS> TRankDictionary;

    TSize wordPos    = _toWordPos(dict, posInBlock);
    TSize posInWord  = _toPosInWord(dict, posInBlock);

    TSize valueRank = 0;

    // NOTE(esiragusa): writing the loop in this form prevents the compiler from unrolling it.
    //    for (TSize wordPrevPos = 0; wordPrevPos < wordPos; ++wordPrevPos)
    //      valueRank += _getWordRank(dict, values[wordPrevPos].i, c);

    for (TSize wordPrevPos = 0; wordPrevPos < TRankDictionary::_WORDS_PER_BLOCK; ++wordPrevPos)
        if (wordPrevPos < wordPos)
            valueRank += _getWordRank(dict, values[wordPrevPos].i, TRankDictionary::_VALUES_PER_WORD - 1, c, smaller);

    valueRank += _getWordRank(dict, values[wordPos].i, posInWord, c, smaller);

    return valueRank;
}

// ----------------------------------------------------------------------------
// Function _getValuesRanks()
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _getValuesRanks() for Dna.
// TODO(cpockrandt): why could that be beneficial?

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type
_getValuesRanks(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    typedef typename RankDictionaryBlock_<TValue, Levels<TSpec, TConfig> >::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                                        TValueSize;

    TBlock blockRank;

    for (TValueSize c = 0; c < ValueSize<TValue>::VALUE; ++c)
        assignValue(blockRank, c, _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue(c)));

    return blockRank;
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos>
inline typename RankDictionaryBlock_<TValue, TPREFIXLEVELS>::Type
_getValuesRanks(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
                TPos const pos)
{
    typedef typename RankDictionaryBlock_<TValue, TPREFIXLEVELS>::Type    TBlock;
    typedef typename ValueSize<TValue>::Type                                        TValueSize;

    TBlock blockRank;
    TPos rank = _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue(0));
    assignValue(blockRank, 0, rank);

    for (TValueSize c = 1; c < RankDictionaryBlockSize_<TValue, TPREFIXLEVELS>::VALUE; ++c)
    {
        TPos smaller = 0;
        rank = _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), TValue(c), smaller);
        assignValue(blockRank, c, rank + smaller);
    }

    return blockRank;
}

// ----------------------------------------------------------------------------
// Function _getValuesRanks(bool)
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
inline typename RankDictionaryBlock_<bool, Levels<TSpec, TConfig> >::Type
_getValuesRanks(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    return _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), true);
}

template <typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB, typename TPos>
inline typename RankDictionaryBlock_<bool, TPREFIXLEVELS>::Type
_getValuesRanks(RankDictionary<bool, TPREFIXLEVELS> const & dict,
                TPos const pos)
{
    return _getValueRank(dict, _valuesAt(dict, pos), _toPosInBlock(dict, pos), true);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline
typename std::enable_if_t<TConfig::LEVELS == 3, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos, TChar const c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const                           TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreUltraBlocks>::Type>::Type   TFibreUltraBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreSuperBlocks>::Type>::Type   TFibreSuperBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type         TFibreBlocks;
    typedef typename Size<TRankDictionary>::Type                                            TSize;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreUltraBlocks const & ultraBlock(_ultraBlockAt(dict, pos));
    TFibreSuperBlocks const & superBlock(_superBlockAt(dict, pos));
    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c))
         + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline
typename std::enable_if_t<TConfig::LEVELS == 2, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos, TChar const c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const                           TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreSuperBlocks>::Type>::Type   TFibreSuperBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type         TFibreBlocks;
    typedef typename Size<TRankDictionary>::Type                                            TSize;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreSuperBlocks const & superBlock(_superBlockAt(dict, pos));
    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
         + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline
typename std::enable_if_t<TConfig::LEVELS == 1, typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > const>::Type>
getRank(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos, TChar const c)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > const                           TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type         TFibreBlocks;
    typedef typename Size<TRankDictionary>::Type                                            TSize;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
         + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos, typename TChar>
inline TSize getRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict, TPos const pos, TChar const c)
{
    TPos smaller;
    return getRank(dict, pos, static_cast<TValue>(c), smaller);
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos, typename TChar>
inline TSize getRank(RankDictionary<TValue, Levels<TSpec, LevelsRDConfig<TSize, TFibre, LEVELS, WPB> > > const & dict,
                     TPos const pos,
                     TChar const c,
                     TPos & /*smaller*/)
{
    return getRank(dict, pos, static_cast<TValue>(c));
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos, typename TChar>
inline typename std::enable_if_t<LEVELS == 3, TSize>
getRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
        TPos const pos,
        TChar const c,
        TPos & smaller)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS> const   TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreUltraBlocks>::Type>::Type           TFibreUltraBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreSuperBlocks>::Type>::Type           TFibreSuperBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type                 TFibreBlocks;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreUltraBlocks const & ultraBlock(_ultraBlockAt(dict, pos));
    TFibreSuperBlocks const & superBlock(_superBlockAt(dict, pos));
    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    smaller = 0;
    if (SEQAN_LIKELY(ordValue(c) > 0 && ordValue(c) < ValueSize<TValue>::VALUE - 1))
    {
        return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c), smaller)
             + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c), smaller)
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);
    }
    else if (ordValue(c) == 0)
    {
        return _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(c))
             + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
    }
    else
    {
        smaller = _getUltraBlockRank(dict, ultraBlock, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                  + _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                  + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                  + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(ValueSize<TValue>::VALUE - 2));
        return pos + 1 - smaller;
    }
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos, typename TChar>
inline typename std::enable_if_t<LEVELS == 2, TSize>
getRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
        TPos const pos,
        TChar const c,
        TPos & smaller)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS> const   TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreSuperBlocks>::Type>::Type           TFibreSuperBlocks;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type                 TFibreBlocks;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreSuperBlocks const & superBlock(_superBlockAt(dict, pos));
    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    smaller = 0;
    if (SEQAN_LIKELY(ordValue(c) > 0 && ordValue(c) < ValueSize<TValue>::VALUE - 1))
    {
	return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c), smaller)
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);
    }
    else if (ordValue(c) == 0)
    {
        return _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(c))
             + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
    }
    else
    {
        smaller = _getSuperBlockRank(dict, superBlock, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                + _getBlockRank(dict, entry.block, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(ValueSize<TValue>::VALUE - 2));
        return pos + 1 - smaller;
    }
}

template <typename TValue, typename TSpec, typename TSize, typename TFibre, unsigned LEVELS, unsigned WPB,
          typename TPos, typename TChar>
inline typename std::enable_if_t<LEVELS == 1, TSize>
getRank(RankDictionary<TValue, TPREFIXLEVELS> const & dict,
        TPos const pos,
        TChar const c,
        TPos & smaller)
{
    typedef RankDictionary<TValue, TPREFIXLEVELS> const   TRankDictionary;
    typedef typename Value<typename Fibre<TRankDictionary, FibreRanks>::Type>::Type                 TFibreBlocks;

    TSize posInBlock = _toPosInBlock(dict, pos);

    TFibreBlocks const & entry(dict.blocks[_toBlockPos(dict, pos)]);

    smaller = 0;
    if (SEQAN_LIKELY(ordValue(c) > 0 && ordValue(c) < (ValueSize<TValue>::VALUE - 1)))
    {
        return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c), smaller)
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c), smaller);
    }
    else if (ordValue(c) == 0)
    {
        return _getBlockRank(dict, entry.block, pos, static_cast<TValue>(c))
             + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(c));
    }
    else
    {
        smaller = _getBlockRank(dict, entry.block, pos, static_cast<TValue>(ValueSize<TValue>::VALUE - 2))
                + _getValueRank(dict, entry.values, posInBlock, static_cast<TValue>(ValueSize<TValue>::VALUE - 2));
        return pos + 1 - smaller;
    }
}

template <typename TSpec, typename TConfig, typename TPos>
inline typename Size<RankDictionary<bool, Levels<TSpec, TConfig> > const>::Type
getRank(RankDictionary<bool, Levels<TSpec, TConfig> > const & dict, TPos const pos)
{
    return getRank(dict, pos, true);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------
template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline typename Value<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
getValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos)
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
getValue(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict, TPos const pos)
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
inline void setValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TPos const pos, TChar const c)
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
// TODO(cpockrandt): try to remove this (it is currently still used by wavelet trees)
template <typename TValue, typename TSpec, typename TConfig, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TChar const c, Tag<TExpand> const tag)
{
    resize(dict, length(dict) + 1, tag);
    setValue(dict, length(dict) - 1, c);
}

// ----------------------------------------------------------------------------
// Function updateRanks()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize>
inline typename std::enable_if_t<TConfig::LEVELS < 3, void>
_updateUltraBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > & /* dict */, TSize const /* blockPos */) {}

template <typename TValue, typename TSpec, typename TConfig, typename TSize>
inline typename std::enable_if_t<TConfig::LEVELS == 3, void>
_updateUltraBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize const blockPos)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDictionary;
    static const unsigned blocks_per_ub = TRankDictionary::_VALUES_PER_ULTRABLOCK / TRankDictionary::_VALUES_PER_BLOCK;

    if (blockPos % blocks_per_ub == 0)
    {
        TSize curr   = _toPos(dict, blockPos);
        TSize prevUB = _toPos(dict, blockPos - blocks_per_ub);

        _ultraBlockAt(dict, curr) = _ultraBlockAt(dict, prevUB) + _superBlockAt(dict, curr);
        _clearSuperBlockAt(dict, curr);
    }
}

template <typename TValue, typename TSpec, typename TConfig, typename TSize>
inline typename std::enable_if_t<TConfig::LEVELS == 1, void>
_updateSuperBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > & /* dict */, TSize const /* blockPos */) {}

template <typename TValue, typename TSpec, typename TConfig, typename TSize>
inline typename std::enable_if_t<TConfig::LEVELS >= 2, void>
_updateSuperBlock(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize const blockPos)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDictionary;
    static const unsigned blocks_per_sb = TRankDictionary::_VALUES_PER_SUPERBLOCK / TRankDictionary::_VALUES_PER_BLOCK;

    if (blockPos % blocks_per_sb == 0)
    {
        TSize curr   = _toPos(dict, blockPos);
        TSize prevSB = _toPos(dict, blockPos - blocks_per_sb);

        _superBlockAt(dict, curr) = _superBlockAt(dict, prevSB) + _blockAt(dict, curr);
        _clearBlockAt(dict, curr);
    }
}

template <typename TValue, typename TSpec, typename TConfig>
inline void
updateRanks(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> >         TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                    TSize;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type       TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type          TBlockIter;

    if (empty(dict)) return;

    // Insures the first block/superbock/ultrablock ranks start from zero.
    _clearUltraBlockAt(dict, 0u);
    _clearSuperBlockAt(dict, 0u);
    _clearBlockAt(dict, 0u);

    // Clear the uninitialized values.
    _padValues(dict);

    TBlockIter blocksBegin = begin(dict.blocks, Standard());
    TBlockIter blocksEnd = end(dict.blocks, Standard());

    for (TBlockIter blocksIt = blocksBegin; blocksIt != blocksEnd; ++blocksIt)
    {
        TSize blockPos = blocksIt - blocksBegin;
        TSize curr = _toPos(dict, blockPos);
        TSize next = _toPos(dict, blockPos + 1);

        if (blockPos > 0)
        {
            _updateSuperBlock(dict, blockPos);
            _updateUltraBlock(dict, blockPos);
        }

        if (blocksIt != blocksEnd - 1)
        {
            _blockAt(dict, next) = _blockAt(dict, curr) + _getValuesRanks(dict, next - 1);
        }
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
reserve(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize const newCapacity, Tag<TExpand> const tag)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDict_;

    if (TConfig::LEVELS > 1)
    {
        reserve(dict.superblocks,
                (newCapacity + TRankDict_::_VALUES_PER_SUPERBLOCK - 1) / TRankDict_::_VALUES_PER_SUPERBLOCK,
                tag);
    }
    if (TConfig::LEVELS > 2)
    {
        reserve(dict.ultrablocks,
                (newCapacity + TRankDict_::_VALUES_PER_ULTRABLOCK - 1) / TRankDict_::_VALUES_PER_ULTRABLOCK,
                tag);
    }
    return reserve(dict.blocks, (newCapacity + TRankDict_::_VALUES_PER_BLOCK - 1) / TRankDict_::_VALUES_PER_BLOCK, tag);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Levels<TSpec, TConfig> > >::Type
resize(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, TSize const newLength, Tag<TExpand> const tag)
{
    typedef RankDictionary<TValue, Levels<TSpec, TConfig> > TRankDict_;

    dict._length = newLength;
    if (TConfig::LEVELS > 1)
    {
        resize(dict.superblocks,
               (newLength + TRankDict_::_VALUES_PER_SUPERBLOCK - 1) / TRankDict_::_VALUES_PER_SUPERBLOCK,
               tag);
    }
    if (TConfig::LEVELS > 2)
    {
        resize(dict.ultrablocks,
               (newLength + TRankDict_::_VALUES_PER_ULTRABLOCK - 1) / TRankDict_::_VALUES_PER_ULTRABLOCK,
               tag);
    }
    return resize(dict.blocks, (newLength + TRankDict_::_VALUES_PER_BLOCK - 1) / TRankDict_::_VALUES_PER_BLOCK, tag);
}

template <typename TValue, typename TSpec, typename TConfig>
inline bool save(RankDictionary<TValue, Levels<TSpec, TConfig> > const & dict,
                 const char * fileName,
                 int const openMode)
{
    String<char> name;
    bool result = true;

    result &= save(getFibre(dict, FibreRanks()), toCString(fileName), openMode);

    if (TConfig::LEVELS > 1)
    {
        name = fileName;
        append(name, ".sbl");
        result &= save(getFibre(dict, FibreSuperBlocks()), toCString(name), openMode);
    }

    if (TConfig::LEVELS > 2)
    {
        name = fileName;
        append(name, ".ubl");
        result &= save(getFibre(dict, FibreUltraBlocks()), toCString(name), openMode);
    }

    return result;
}

template <typename TValue, typename TSpec, typename TConfig>
inline bool open(RankDictionary<TValue, Levels<TSpec, TConfig> > & dict, const char * fileName, int const openMode)
{
    String<char> name;
    bool result = true;

    name = fileName;
    append(name, "");
    result &= open(getFibre(dict, FibreRanks()), toCString(name), openMode);

    if (TConfig::LEVELS > 1)
    {
        name = fileName;
        append(name, ".sbl");
        result &= open(getFibre(dict, FibreSuperBlocks()), toCString(name), openMode);
    }

    if (TConfig::LEVELS > 2)
    {
        name = fileName;
        append(name, ".ubl");
        result &= open(getFibre(dict, FibreUltraBlocks()), toCString(name), openMode);
    }

    return result;
}

}

#undef TPREFIXLEVELS

#endif  // INDEX_FM_RANK_DICTIONARY_LEVELS_H_
