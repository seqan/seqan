// ==========================================================================
//                                  FMIndex
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
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

#ifndef INDEX_FM_RANK_SUPPORT_BIT_STRING_H_
#define INDEX_FM_RANK_SUPPORT_BIT_STRING_H_

#include <seqan/index_fm.h>
#include <seqan/misc/misc_bit_twiddling.h>

namespace seqan {

//Rank Support Bit String
template <typename TSpec = void>
struct RankSupportBitString;

// ==========================================================================
// Tags
// ==========================================================================

// FM index fibres

/**
.Tag.RankSupportBitString Fibres
..summary:Tag to select a specific fibre (e.g. table, object, ...) of a @Class.RankSupportBitString@.
..remarks:These tags can be used to get @Metafunction.Fibre.Fibres@ of a rank support bit string.
..cat:Index

..tag.FibreBits:The bit string.
..tag.FibreBlocks:The block string.
..tag.FibreSuperBlocks:The super block string.

..see:Metafunction.Fibre
..see:Function.getFibre
..include:seqan/index.h
*/

struct FibreBits_;
struct FibreBlocks_;
struct FibreSuperBlocks_;

typedef Tag<FibreBits_> const           FibreBits;         
typedef Tag<FibreBlocks_> const         FibreBlocks;      
typedef Tag<FibreSuperBlocks_> const    FibreSuperBlocks; 

typedef FibreBits           RankSupportBitStringBits;
typedef FibreBlocks         RankSupportBitStringBlocks;
typedef FibreSuperBlocks    RankSupportBitStringSuperBlocks;

// ==========================================================================
// Metafunctions
// ==========================================================================

template <typename TSpec>
struct DefaultOverflowImplicit<RankSupportBitString<TSpec> >
{
    typedef Generous Type;
};

//The limiting factor of the size is the underlying data type of the super block string
template <typename TSpec>
struct Size<RankSupportBitString<TSpec> >
{
    typedef unsigned long long Type;
};

template <typename TSpec>
struct Position<RankSupportBitString<TSpec> > :
	Size<RankSupportBitString<TSpec> > {};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBits>
{
    typedef String<unsigned long long> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreBlocks>
{
    typedef String<unsigned short> Type;
};

template <typename TSpec>
struct Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>
{
    typedef String<typename Size<RankSupportBitString<TSpec> >::Type> Type;
};

// ==========================================================================
// Classes
// ==========================================================================

/**
.Class.RankSupportBitString:
..summary:A bit string supporting rank queries in constant time.
..cat:Index
..signature:RankSupportBitString<TSpec>
..param.TSpec:Specialisation tag.
...default:void
..remarks:The constant rank query time is achieved by evaluating precomputed subsolutions. In order to do so, the bit string is divided into blocks of length l. A super block string stores for each block of l blocks the number of bits set from the beginning. In addition a block string stores the number of bits set in each block from the start of the last super block block. Therefore it is possible to compute the result of a rank query in constant time by adding information from the bit, block and super block string.
..include:seqan/index.h
*/

// Forward declaration because setBit is used in constructor
template <typename TSpec, typename TPos, typename TBit>
inline void setBit(RankSupportBitString<TSpec> & bitString, TPos const pos, TBit const setBit);

template <typename TSpec>
struct RankSupportBitString
{
    typedef typename Fibre<RankSupportBitString, FibreBits>::Type         TBitString;
    typedef typename Fibre<RankSupportBitString, FibreBlocks>::Type       TBlockString;
    typedef typename Fibre<RankSupportBitString, FibreSuperBlocks>::Type  TSuperBlockString;

    TBitString                                  bits;
    TBlockString                                blocks;
    TSuperBlockString                           superBlocks;
    typename Size<RankSupportBitString>::Type   length_;

    RankSupportBitString() :
        length_(0)
    {}

    template <typename TString>
    RankSupportBitString(TString const & input) :
        length_(length(input))
    {
        typedef typename Iterator<TString const>::Type TIter_;

        resize(*this, length_);
        typename Position<RankSupportBitString>::Type i = 0;
        for (TIter_ it = begin(input, Standard()); it != end(input, Standard()); ++it, ++i)
            setBit(*this, i, getValue(it));
        updateRanks_(*this);
    }

    inline bool operator==(const RankSupportBitString & other) const
    {
        return length_ == other.length_ &&
               bits == other.bits &&
               blocks == other.blocks &&
               superBlocks == other.superBlocks;
    }
};

// ==========================================================================
// Functions
// ==========================================================================

/**
.Function.appendValue:
..param.target:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TBit>
inline void appendValue(RankSupportBitString<TSpec> & bitString, TBit const bit)
{
	typename Size<RankSupportBitString<TSpec> >::Type len = length(bitString);
    resize(bitString, len + 1);
    setBit(bitString, len, bit);
    updateLastRank_(bitString); 
}

// ==========================================================================
/**
.Function.clear
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline void clear(RankSupportBitString<TSpec> & bitString)
{
    clear(bitString.bits);
    clear(bitString.blocks);
    clear(bitString.superBlocks);
    bitString.length_ = 0;
}

// ==========================================================================
// This function returns the number of bits set in a block until a specified position
template <typename TValue>
inline unsigned getRankInBlock_(TValue const value, False)
{
    return popCount(value);
}

template <typename TValue>
inline unsigned getRankInBlock_(TValue const value, True)
{
    return popCount(value);
}

template <typename TValue>
inline unsigned getRankInBlock_(TValue const value)
{
    return getRankInBlock_(value, typename Eval<(BitsPerValue<TValue>::VALUE > 32)>::Type());
}

template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
getRankInBlock_(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;

    TFibreBitsValue const mask = ((TFibreBitsValue)2u << getPosInBlock_(bitString, pos)) - 1;
    return getRankInBlock_(bitString.bits[getBlockPos_(bitString, pos)] & mask);
}

// ==========================================================================
// This function returns the number of bits set in a block until a specified position
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
getSuperBlockPos_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                 TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type      TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type    TFibrelocks;
    typedef typename Value<TFibreBits>::Type                            TFibreBitsValue;
    typedef typename Value<TFibrelocks>::Type                     		TFibreBlocksValue;

    TFibreBlocksValue const bitsPerValue_ = BitsPerValue<TFibreBitsValue>::VALUE;
    return pos / (bitsPerValue_ * bitsPerValue_);
}

// ==========================================================================
/**
.Function.getRank
..summary:Returns the rank (the number of bits set from the start of the bit string) of a specified position.
..signature:getRank(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of a bit.
..returns:Value type of the super block fibre (default unsigned long).
..include:seqan/index.h
..example.code:
*/
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
getRank(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    return bitString.superBlocks[getSuperBlockPos_(bitString, pos)]
         + bitString.blocks[getBlockPos_(bitString, pos)]
         + getRankInBlock_(bitString, pos);
}


/**
.Function.empty
..param.object:
...type:Class.RankSupportBitString
*/
// template <typename TSpec, typename TPos>
// inline bool
// empty(RankSupportBitString<TSpec> const & bitString, TPos const pos)
// {
//     return 
// }


// ==========================================================================
/**
.Function.getBit
..summary:Returns whether a specified bit is set or not.
..signature:getBit(bitString, pos)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..returns:Returns whether a specified bit is set or not.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));
...
mark all 'a's
...

for (unsigned i = 0; i < length(bitString); ++i)
    if(getBit(bitString, i))
        std::cout << "a found at: " << i << std::endl;
*/

template <typename TSpec, typename TPos>
inline bool getBit(RankSupportBitString<TSpec> const & bitString, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
   
    return (bitString.bits[getBlockPos_(bitString, pos)] & ((TFibreBitsValue)1u << getPosInBlock_(bitString, pos))) != (TFibreBitsValue)0;
}

// ==========================================================================
// This function returns the position in the block string of the corresponding block.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
getBlockPos_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type  TFibreBits;
    typedef typename Value<TFibreBits>::Type                        TFibreBitsValue;

    return pos / BitsPerValue<TFibreBitsValue>::VALUE;
}

// ==========================================================================
/**
.Function.getFibre
..param.container:
...type:Class.RankSupportBitString
..param.fibreTag:
...type:Tag.RankSupportBitString Fibres
*/
template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBits)
{
    return string.bits;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBits)
{
    return string.bits;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreBlocks)
{
    return string.blocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type const &
getFibre(RankSupportBitString<TSpec> const & string, const FibreBlocks)
{
    return string.blocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type &
getFibre(RankSupportBitString<TSpec>&string, const FibreSuperBlocks)
{
    return string.superBlocks;
}

template <typename TSpec>
inline typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type const &
getFibre(const RankSupportBitString<TSpec>&string, const FibreSuperBlocks)
{
    return string.superBlocks;
}

// ==========================================================================
// This function returns the position of a specified bit within a block.
template <typename TSpec, typename TPos>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type>::Type
getPosInBlock_(RankSupportBitString<TSpec> const & /*bitString*/, TPos const pos)
{
    typedef RankSupportBitString<TSpec>                             TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type  TFibreBits;
    typedef typename Value<TFibreBits>::Type                        TFibreBitsValue;

    return pos % BitsPerValue<TFibreBitsValue>::VALUE;
}

// ==========================================================================
/**
.Function.length
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec>
inline typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type
length(RankSupportBitString<TSpec> const & bitString)
{
    return bitString.length_;
}


// ==========================================================================
/*
.Function.updateRanks_
..summary:Adds the block and super block information to the bit string.
..signature:updateRanks_(bitString)
..param.bitString:The bit string to be completed.
...type:Class.RankSupportBitString
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBit(bitString, 1);

updateRanks_(bitString);
*/

template <typename TSpec, typename TPos>
inline void updateRanksImpl_(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    if (!empty(bitString))
    {
        typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
        typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
        typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
        typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
        typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
        typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
        typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

        TFibreSuperBlocksValue i = getBlockPos_(bitString, pos);

        TFibreBlocksValue blockSum_;
        if (i == 0)
            blockSum_ = 0;
        else
            blockSum_ = bitString.blocks[i - 1];

        TFibreSuperBlocksValue sBlockSum_;
        TFibreSuperBlocksValue superBlockPos = getSuperBlockPos_(bitString, pos);
        if (superBlockPos != 0)
        {
            sBlockSum_ = bitString.superBlocks[superBlockPos];
        }
        else
            sBlockSum_ = 0;
        
        if (i == 0)
            ++i;
        for (; i < length(bitString.bits); ++i)
        {
        	blockSum_ += getRankInBlock_(bitString.bits[i - 1]);
            if ((i % BitsPerValue<TFibreBitsValue>::VALUE) == 0u)
            {
                sBlockSum_ += blockSum_;
                bitString.superBlocks[++superBlockPos] = sBlockSum_;
                blockSum_ = 0;
            }
            bitString.blocks[i] = blockSum_;
       }
    }
}

template <typename TSpec, typename TPos>
inline void updateRanks_(RankSupportBitString<TSpec> & bitString, TPos pos)
{
    updateRanksImpl_(bitString, pos);
}

template <typename TSpec>
inline void updateRanks_(RankSupportBitString<TSpec> & bitString)
{
    updateRanksImpl_(bitString, 0u);
}

// This function update the rank information of the last block.
// ==========================================================================
template <typename TSpec>
inline void updateLastRank_(RankSupportBitString<TSpec> & bitString)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    // It is only necessary to update the rank at the start of a new block
    TFibreSuperBlocksValue pos = length(bitString) - 1;
    if (getPosInBlock_(bitString, pos) != 0u) return;

    // Here we compute the new rank if the last bit does not reside in the 
    // first block and the last as well as the second last bit are in the 
    // same block.
    TFibreSuperBlocksValue superBlockPos = getSuperBlockPos_(bitString, pos);
    TFibreSuperBlocksValue blockPos = getBlockPos_(bitString, pos);
    if (blockPos > 0u && superBlockPos == getSuperBlockPos_(bitString, pos - 1))
        getFibre(bitString, FibreBlocks())[blockPos] = getFibre(bitString, FibreBlocks())[blockPos - 1]
                + getRankInBlock_(bitString.bits[blockPos - 1]);

    // Here we compute the new rank if the last and second last bit are in
    // different blocks.
    if ((superBlockPos > getSuperBlockPos_(bitString, pos - 1)) && superBlockPos > 0u)
        getFibre(bitString, FibreSuperBlocks())[superBlockPos] = getFibre(bitString, FibreSuperBlocks())[superBlockPos - 1]
                + getFibre(bitString, FibreBlocks())[blockPos - 1]
                + getRankInBlock_(bitString.bits[blockPos - 1]);
}


// ==========================================================================
template <typename TSpec, typename TSize, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
reserve(RankSupportBitString<TSpec> & bitString, TSize const size, Tag<TExpand> const tag)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type          TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    TFibreBlocksValue bitsPerBlock_ = BitsPerValue<TFibreBitsValue>::VALUE;
    TFibreSuperBlocksValue numberOfBlocks_ = (size + bitsPerBlock_ - 1) / bitsPerBlock_;

    resserve(bitString.blocks, numberOfBlocks_, tag);
    reserve(bitString.superBlocks, (numberOfBlocks_ + bitsPerBlock_ - 1) / bitsPerBlock_, tag);
    return reserve(bitString.bits, numberOfBlocks_, tag) * bitsPerBlock_;
    //reserve(bitString, size);
}

// ==========================================================================
/**
.Function.resize
..param.object:
...type:Class.RankSupportBitString
*/
template <typename TSpec, typename TLength, typename TValue, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
resize(RankSupportBitString<TSpec> & bitString, TLength const length_, TValue const value_, Tag<TExpand> const tag)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<TRankSupportBitString, FibreBits>::Type          TFibreBits;
    typedef typename Fibre<TRankSupportBitString, FibreBlocks>::Type        TFibreBlocks;
    typedef typename Fibre<TRankSupportBitString, FibreSuperBlocks>::Type   TFibreSuperBlocks;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;
    typedef typename Value<TFibreBlocks>::Type                              TFibreBlocksValue;
    typedef typename Value<TFibreSuperBlocks>::Type                         TFibreSuperBlocksValue;

    TFibreBlocksValue bitsPerBlock_ = BitsPerValue<TFibreBitsValue>::VALUE;
    TFibreSuperBlocksValue numberOfBlocks_ = (length_ + bitsPerBlock_ - 1) / bitsPerBlock_;

    TLength currentLength = length(bitString);

    resize(bitString.bits, numberOfBlocks_, 0, tag);
    resize(bitString.blocks, numberOfBlocks_, 0, tag);
    resize(bitString.superBlocks, (numberOfBlocks_ + bitsPerBlock_ - 1) / bitsPerBlock_, 0, tag);

    if (value_ != 0 && currentLength < length_)
    {
        for (unsigned i = currentLength; i < (typename MakeUnsigned<TLength const>::Type)length_; ++i)
            setBit(bitString, i, 1);
        updateRanks_(bitString, currentLength);
    }

    return bitString.length_ = length_;
}

template <typename TSpec, typename TLength, typename TExpand>
inline typename Size<RankSupportBitString<TSpec> >::Type
resize(RankSupportBitString<TSpec> & bitString, TLength const length_, Tag<TExpand> const tag)
{
    return resize(bitString, length_, 0, tag);
}

// ==========================================================================
/**
.Function.setBit
..summary:Set a specified bit to true or false.
..signature:setBit(bitString, pos, bit)
..param.bitString:The bit string.
...type:Class.RankSupportBitString
..param.pos:Position of the bit.
..param.bit:The value of the bit.
...remarks:Note that values different from 0 are interpreted as 1.
..include:seqan/index.h
..example.code:
String<Dna5> genome = "ACGTACGT";

RankSupportBitString<> bitString;
resize(bitString, length(genome));

for (unsigned i = 0; i < length(genome); ++i)
    if(genome[i] < Dna5('c'))
        setBit(bitString, 1);

updateRanks_(bitString);
*/
template <typename TSpec, typename TPos, typename TBit>
inline void setBit(RankSupportBitString<TSpec> & bitString, TPos const pos, TBit const newBit)
{
    typedef RankSupportBitString<TSpec>                                     TRankSupportBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type    TFibreBits;
    typedef typename Value<TFibreBits>::Type                                TFibreBitsValue;

    TFibreBitsValue const shiftValue = (TFibreBitsValue)1u << getPosInBlock_(bitString, pos);
    if (newBit != (TBit)0)
        bitString.bits[getBlockPos_(bitString, pos)] |= shiftValue;
    else
        bitString.bits[getBlockPos_(bitString, pos)] &= ~shiftValue;
}

// ==========================================================================
template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);

    String<char> name;
    name = fileName;    append(name, ".bit");
    if (!open(getFibre(string, FibreBits()), toCString(name), openMode))
        return false;

    name = fileName;    append(name, ".bl");    open(getFibre(string, FibreBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".sbl");   open(getFibre(string, FibreSuperBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".len");    open(lengthString, toCString(name), openMode);
    string.length_ = lengthString[0];
    return true;
}

template <typename TSpec>
inline bool open(
    RankSupportBitString<TSpec> & string,
    const char * fileName)
{
    return open(string, fileName, OPEN_RDONLY);
}

template <typename TSpec, typename TSetSpec>
inline bool open(
    StringSet<RankSupportBitString<TSpec>, TSetSpec> & strings,
    const char * fileName)
{
    return open(strings, fileName, OPEN_RDONLY);
}

// ==========================================================================
template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName,
    int openMode)
{
    typedef typename Value<typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type>::Type TValue;
    String<TValue> lengthString;
    resize(lengthString, 1);
    lengthString[0] = length(string);
    String<char> name;
    name = fileName;    append(name, ".len");    save(lengthString, toCString(name), openMode);
    name = fileName;    append(name, ".bit");       save(getFibre(string, FibreBits()), toCString(name), openMode);
    name = fileName;    append(name, ".bl");    save(getFibre(string, FibreBlocks()), toCString(name), openMode);
    name = fileName;    append(name, ".sbl");   save(getFibre(string, FibreSuperBlocks()), toCString(name), openMode);
    return true;
}

template <typename TSpec>
inline bool save(
    RankSupportBitString<TSpec> const & string,
    const char * fileName)
{
    return save(string, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

template <typename TSpec, typename TSetSpec>
inline bool save(
    StringSet<RankSupportBitString<TSpec>, TSetSpec> const & strings,
    const char * fileName)
{
    return save(strings, fileName, DefaultOpenMode<RankSupportBitString<TSpec> >::VALUE);
}

// ==========================================================================
template <typename TValue>
inline void printBits(TValue entrie)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    TValue one = 1;
    std::cout << "entrie: " << entrie << std::endl;
    std::cout << bitsPerValue << std::endl;
    for (TValue i = 0; i < bitsPerValue; ++i)
    {
        std::cout << ((entrie >> i) & one);
    }
    std::cout << std::endl;
}

template <typename TValue, typename TSize>
inline std::ostream & printBits(std::ostream & stream, TValue entrie, TSize blockSize)
{
    unsigned bitsPerValue = BitsPerValue<TValue>::VALUE;
    bool temp;
    for (int i = bitsPerValue - 1; i >= 0; --i)
    {
        temp = (entrie >> i) & 1;
        stream << temp;
        if ((bitsPerValue - i) % blockSize == 0)
            stream << " ";
    }
    return stream;
}

template <typename TSpec>
inline std::ostream & operator<<(std::ostream & stream, const RankSupportBitString<TSpec> & rankSupportBitString)
{
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBits>::Type                TBitString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreBlocks>::Type             TBlockString;
    typedef typename Fibre<RankSupportBitString<TSpec>, FibreSuperBlocks>::Type        TSuperBlockString;

    typedef typename Value<TBitString>::Type          TBitStringValue;
    typedef typename Value<TBlockString>::Type       TBlockStringValue;
    typedef typename Value<TSuperBlockString>::Type  TSuperBlockStringValue;

    unsigned bitsPerBlock = BitsPerValue<typename Value<TBitString>::Type>::VALUE;

    TBitString const & bits = rankSupportBitString.bits;
    TBlockString const & blockString = rankSupportBitString.blocks;
    TSuperBlockString const & superBlocks = rankSupportBitString.superBlocks;

    stream << "  ";
    for (TBitStringValue i = 0; i < length(bits); i++)
    {
        printBits(stream, bits[i], bitsPerBlock);
    }
    stream << std::endl;

    for (TBlockStringValue i = 0; i < length(blockString); i++)
    {
        stream << blockString[i] << " ";
    }
    stream << std::endl;

    for (TSuperBlockStringValue i = 0; i < length(superBlocks); i++)
    {
        stream << superBlocks[i] << " ";
    }
    return stream;
}

}


#endif  // INDEX_FM_BRANK_SUPPORT_BIT_STRING_H_
